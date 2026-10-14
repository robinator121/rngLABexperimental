"""Predeclared finite-horizon Bonferroni monitor for an IID fair-bit null."""
import numpy as np
from scipy.stats import binomtest


def binomial_p(k,n,p=.5):
    return float(binomtest(int(k),int(n),p).pvalue) if n else 1.


def test_bits(bits, uniforms=None):
    b=np.asarray(bits,dtype=np.uint8);n=len(b)
    if n<64: raise ValueError('Monitoring requires at least 64 bits')
    result={'frequency':binomial_p(b.sum(),n),
            'runs':binomial_p(np.count_nonzero(b[1:]!=b[:-1]),n-1)}
    # Disjoint pairs: each category count is exactly binomial under the null.
    pair=b[:n//2*2].reshape(-1,2);codes=2*pair[:,0]+pair[:,1]
    result['serial']=min(1.,4*min(binomial_p(np.count_nonzero(codes==k),len(codes),.25) for k in range(4)))
    for lag in [2,8,32]:
        blocks=b[:n//(2*lag)*(2*lag)].reshape(-1,2*lag)
        equal=(blocks[:,:lag]==blocks[:,lag:]).ravel()
        result[f'lag{lag}']=binomial_p(equal.sum(),len(equal))
    result['period4']=min(1.,4*min(binomial_p(b[j::4].sum(),len(b[j::4])) for j in range(4)))
    if uniforms is not None and len(uniforms):
        counts,_=np.histogram(uniforms,bins=16,range=(0,1))
        result['uniform']=min(1.,16*min(binomial_p(k,len(uniforms),1/16) for k in counts))
    return result


class BitMonitor:
    def __init__(self,alpha=.01,max_looks=32,window_bits=65536):
        if not 0<alpha<1 or max_looks<1 or window_bits<64:raise ValueError('invalid monitor settings')
        self.alpha,self.max_looks,self.window_bits=alpha,max_looks,window_bits
        self.looks=0;self.buffer=np.empty(0,dtype=np.uint8);self.records=[]
        self.first_alarm=None;self.consumed=0

    def update(self,bits,time,uniforms=None):
        self.looks+=1
        if self.looks>self.max_looks:raise ValueError('predeclared monitor horizon exceeded')
        self.consumed+=len(bits)
        self.buffer=np.concatenate((self.buffer,bits))[-self.window_bits:]
        pvals=test_bits(self.buffer,uniforms)
        # 8 tests reserved even when no uniform-variate test is available.
        cutoff=self.alpha/(self.max_looks*8)
        alarm=any(p<=cutoff for p in pvals.values())
        if alarm and self.first_alarm is None:self.first_alarm=int(time)
        record={'time':int(time),'consumed_bits':self.consumed,'window_bits':len(self.buffer),
                'bit_bias':float(self.buffer.mean()-.5),'alarm':bool(alarm),'cutoff':cutoff,
                'min_p':float(min(pvals.values())),**{f'p_{k}':v for k,v in pvals.items()}}
        self.records.append(record)
        return alarm
