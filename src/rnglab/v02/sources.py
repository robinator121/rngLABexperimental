"""Source schedules, actual-input auditing, combiners and finite file-backed bits."""
from pathlib import Path
from functools import lru_cache
import hashlib,hmac
import numpy as np
from rnglab.providers import SyntheticSource,PCGSource,DirectProvider,HybridProvider,ChaChaSource,bits_to_uniform
from rnglab.v02.monitor import BitMonitor


class ScheduledSource:
    def __init__(self,seed,faults=()):
        self.base=SyntheticSource(seed)
        self.position=0
        self.faults=[]
        for i,f in enumerate(faults):
            f=dict(f)
            if f.get('kind') not in {'bias','freeze','periodic','autocorrelation','partial_loss'}:raise ValueError('unknown scheduled fault')
            if not 0<=f.get('strength',1)<=1:raise ValueError('strength outside [0,1]')
            for key in ['onset_bit','ramp_bits','recovery_bits','active_bits','period_bits','end_bit','drift_period_bits']:
                if key in f and (not isinstance(f[key],int) or f[key]<0):raise ValueError(f'invalid {key}')
            if f.get('end_bit',2**63)<=f.get('onset_bit',0):raise ValueError('end_bit must exceed onset')
            if 'period_bits' in f and (f['period_bits']<1 or not 0<f.get('active_bits',0)<=f['period_bits']):raise ValueError('invalid intermittent duty cycle')
            if not f.get('pattern','0011') or set(f.get('pattern','0011'))-{'0','1'}:raise ValueError('invalid periodic pattern')
            if f.get('frozen_bit',0) not in [0,1]:raise ValueError('invalid frozen bit')
            if not 0<=f.get('known_bits',0)<=f.get('block_bits',32) or f.get('block_bits',32)<1:raise ValueError('invalid partial loss')
            self.faults.append((f,np.random.Generator(np.random.PCG64([*seed,1011,i])),None))

    def read_bits(self,count):
        out=self.base.read_bits(count);pos=np.arange(self.position,self.position+count)
        for index,(f,rng,previous) in enumerate(self.faults):
            onset=f.get('onset_bit',0);ramp=f.get('ramp_bits',0);end=f.get('end_bit')
            envelope=np.clip((pos-onset)/ramp,0,1) if ramp else (pos>=onset).astype(float)
            if end is not None:
                recovery=f.get('recovery_bits',0)
                envelope*=np.clip((end+recovery-pos)/recovery,0,1) if recovery else (pos<end)
            if f.get('period_bits'):
                envelope*=((pos-onset)%f['period_bits']<f['active_bits'])
            if f.get('drift_period_bits'):
                envelope*=.5*(1-np.cos(2*np.pi*(pos-onset)/f['drift_period_bits']))
            mask=rng.random(count)<f.get('strength',1)*envelope
            kind=f['kind']
            if kind=='partial_loss':mask&=pos%f.get('block_bits',32)<f.get('known_bits',0)
            if kind in {'bias','freeze','partial_loss'}:
                out[mask]=1 if kind=='bias' else f.get('frozen_bit',0)
            elif kind=='periodic':
                pattern=np.array(list(f.get('pattern','0011')),dtype=np.uint8)
                out[mask]=pattern[pos[mask]%len(pattern)]
            elif count:
                if previous is None:mask[0]=False
                origin=np.maximum.accumulate(np.where(~mask,np.arange(1,count+1),0))
                out=np.concatenate(([0 if previous is None else previous],out))[origin].astype(np.uint8)
            self.faults[index]=(f,rng,int(out[-1]) if count else previous)
        self.position+=count
        return out


class CombinedSource:
    def __init__(self,sources,mode='xor',relationship='independent'):
        if len(sources)<2 or mode not in {'xor','hash'} or relationship not in {'independent','identical','attacker_copy','attacker_complement'}:raise ValueError('invalid combiner')
        self.sources=sources;self.mode=mode;self.relationship=relationship
        self.position=0;self.pending=np.empty(0,dtype=np.uint8);self.raw_bits_consumed=0

    def _inputs(self,count):
        streams=[s.read_bits(count) for s in self.sources]
        if self.relationship in {'identical','attacker_copy'}:streams[1]=streams[0].copy()
        if self.relationship=='attacker_complement':streams[1]=1-streams[0]
        self.raw_bits_consumed+=count*len(streams)
        return streams

    def read_bits(self,count):
        if count<0:raise ValueError('count must be nonnegative')
        if self.mode=='xor':
            out=np.bitwise_xor.reduce(self._inputs(count))
        else:
            while len(self.pending)<count:
                # Fixed-length, domain-separated hash conditioning; not an entropy proof.
                data=b'rnglab-v02/combine/'+b''.join(np.packbits(x).tobytes() for x in self._inputs(256))
                digest=np.unpackbits(np.frombuffer(hashlib.sha256(data).digest(),dtype=np.uint8))
                self.pending=np.concatenate((self.pending,digest))
            out,self.pending=self.pending[:count].copy(),self.pending[count:]
        self.position+=count
        return out


@lru_cache(maxsize=32)
def file_digest(path,size,mtime):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class FileSource:
    """Read finite measured/synthetic data, with explicit extraction policy; never cycles."""
    def __init__(self,path,encoding,provenance,offset_bits=0,limit_bits=None):
        if not provenance or encoding not in {'packed-msb','ascii01'}:raise ValueError('explicit provenance and encoding required')
        path=Path(path)
        if not path.is_absolute():path=Path(__file__).parents[3]/path
        stat=path.stat()
        self.provenance={'path':str(path.resolve()),'sha256':file_digest(str(path),stat.st_size,stat.st_mtime_ns),'encoding':encoding,'description':provenance,'exhaustion':'raise EOFError','offset_bits':offset_bits}
        self.encoding=encoding
        if encoding=='packed-msb':
            self.raw=np.memmap(path,dtype=np.uint8,mode='r');self.total_bits=len(self.raw)*8
        else:
            text=''.join(path.read_text(encoding='ascii').split())
            if set(text)-{'0','1'}:raise ValueError('nonbinary input')
            self.bits=np.array(list(text),dtype=np.uint8)
            self.total_bits=len(self.bits)
        if offset_bits<0 or offset_bits>self.total_bits:raise ValueError('invalid file bit offset')
        self.start=offset_bits
        self.limit=self.total_bits-offset_bits if limit_bits is None else limit_bits
        if self.limit<0 or self.start+self.limit>self.total_bits:raise ValueError('invalid file slice')
        self.position=0

    def read_bits(self,count):
        if count<0:raise ValueError('negative count')
        if self.position+count>self.limit:raise EOFError('source exhausted; no wraparound')
        start=self.start+self.position
        if self.encoding=='packed-msb':out=np.unpackbits(self.raw[start//8:(start+count+7)//8])[start%8:start%8+count]
        else:out=self.bits[start:start+count].copy()
        self.position+=count
        return out


def source_model(spec,seed):
    if 'sources' in spec:
        sources=[source_model(s,[*seed,i+700]) for i,s in enumerate(spec['sources'])]
        return CombinedSource(sources,spec.get('combine','xor'),spec.get('relationship','independent'))
    if 'file' in spec:return FileSource(spec['file'],spec['encoding'],spec['provenance'],spec.get('offset_bits',0),spec.get('limit_bits'))
    return ScheduledSource(seed,spec.get('faults',[]))


class AuditedSource:
    def __init__(self,source,label,monitor,clock):
        self.source,self.label,self.monitor,self.clock=source,label,monitor,clock
        self.position=0;self.events=[]

    def read_bits(self,count):
        bits=self.source.read_bits(count)
        offset=self.position;self.position+=count
        alarm=self.monitor.update(bits,self.clock())
        self.events.append({'stage':self.label,'time':self.clock(),'source_bit_offset':offset,'bits':count,
                            'sha256':hashlib.sha256(np.packbits(bits).tobytes()).hexdigest(),'alarm':alarm})
        return bits


class ResearchHybrid(HybridProvider):
    def __init__(self,initial,reseed,interval,entropy_bits=256):
        if entropy_bits<64 or entropy_bits%8:raise ValueError('entropy_bits must be a byte multiple >=64')
        self.entropy_bits=entropy_bits
        super().__init__(initial,reseed,interval)

    def _reseed(self):
        offset=self.reseed_source.position
        fresh=np.packbits(self.reseed_source.read_bits(self.entropy_bits)).tobytes()
        self.epoch+=1
        self.key=hmac.new(self.key,b'rnglab-v0.2/reseed/'+self.epoch.to_bytes(8,'big')+fresh,hashlib.sha256).digest()
        self.stream=ChaChaSource(self.key,self.epoch)
        self.events.append({'sample_index':self.position,'epoch':self.epoch,'source_bit_offset':offset,'entropy_bits':self.entropy_bits})


class OutputFaultProvider:
    """Post-CSPRNG faults, explicitly distinct from faults in seed material."""
    def __init__(self,provider,seed,faults):
        self.provider=provider
        # Reuse scheduled faults by replacing only the underlying bit supplier.
        self.fault_source=ScheduledSource(seed,faults)
    @property
    def position(self):return self.provider.position
    @property
    def events(self):return self.provider.events
    def sample(self,count,dimension):
        _,bits=self.provider.sample(count,dimension)
        class Supplied:
            def read_bits(self,n):
                if n!=len(bits):raise ValueError('shape mismatch')
                return bits.copy()
        self.fault_source.base=Supplied()
        bits=self.fault_source.read_bits(len(bits))
        return bits_to_uniform(bits,count,dimension),bits


def research_provider(spec,seed,steps):
    kind=spec.get('provider','synthetic');audits=[];holder={}
    if kind=='prng':provider=DirectProvider(PCGSource(seed))
    elif kind=='synthetic':provider=DirectProvider(source_model(spec,seed))
    elif kind in {'csprng','hybrid','triggered'}:
        interval=spec.get('reseed_interval',2048) if kind=='hybrid' else None
        entropy_bits=spec.get('entropy_bits',256)
        max_reseeds=steps//interval+1 if interval else steps//spec.get('monitor_stride',1024)+1
        monitor=BitMonitor(alpha=.005,max_looks=max_reseeds+1,window_bits=256)
        clock=lambda:holder['provider'].position if holder else 0
        initial=AuditedSource(source_model(spec.get('initial_source',{}),seed),'initial',monitor,clock)
        reseed=AuditedSource(source_model(spec.get('reseed_source',{}),[*seed,9271]),'reseed',monitor,clock)
        provider=ResearchHybrid(initial,reseed,interval,entropy_bits)
        holder['provider']=provider;audits=[initial,reseed]
    else:raise ValueError('unknown research provider')
    base=provider
    if spec.get('output_faults'):
        provider=OutputFaultProvider(provider,[*seed,819],spec['output_faults'])
    return provider,base,audits
