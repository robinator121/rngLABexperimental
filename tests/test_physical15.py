from pathlib import Path
import hashlib
import json
import math
import tempfile
import unittest
from unittest.mock import patch
from scipy.stats import binom
from rnglab15.integrity import capture,digest,lines,snapshot,inventory
from rnglab15.ingest import waveforms,sensors,events,seconds_per_tick
from rnglab15.assessment import Health,Timing,screen
from rnglab15.conditioning import treatments,toeplitz_reference,pack_bits
from rnglab15.applications import FiniteSource,defects,paired_summary,model_run,finite_bead_energy,ising_exact_energy
from rnglab15.pipeline import run


class Physical15Tests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
    def tearDown(self):
        self.tmp.cleanup()
    def file(self,name,data):
        p=self.root/name;p.write_bytes(data);return p
    def test_capture_fragmentation_and_no_overwrite(self):
        p=self.root/'capture'
        r=capture(iter([b'abc',b'\r',b'\ndef']),p,10)
        self.assertEqual(p.read_bytes(),b'abc\r\ndef')
        self.assertEqual(r['sha256'],hashlib.sha256(p.read_bytes()).hexdigest())
        with self.assertRaises(FileExistsError):capture([b'x'],p,2)
    def test_interrupted_capture_retains_partial(self):
        def broken():
            yield b'raw';raise RuntimeError('transport disconnected')
        with self.assertRaises(RuntimeError):capture(broken(),self.root/'run',20)
        self.assertEqual((self.root/'run.partial').read_bytes(),b'raw')
        self.assertFalse((self.root/'run').exists())
    def test_capture_budget_does_not_truncate_silently(self):
        with self.assertRaises(ValueError):capture([b'a',b'bc'],self.root/'run',2)
        self.assertEqual((self.root/'run.partial').read_bytes(),b'a')
    def test_snapshot_roundtrip_and_corruption(self):
        source=self.file('raw',b'a'*2000000)
        target,h=snapshot(source,self.root/'vault')
        self.assertEqual(digest(target),digest(source));self.assertEqual(target.stat().st_mode&0o777,0o444)
        self.assertEqual(snapshot(source,self.root/'vault')[1],h)
        target.chmod(0o644);target.write_bytes(b'corrupt')
        with self.assertRaises(ValueError):snapshot(source,self.root/'vault')
    def test_symlinks_rejected(self):
        p=self.file('source',b'abc');link=self.root/'link';link.symlink_to(p)
        with self.assertRaises(ValueError):digest(link)
    def test_oversized_line_bounded(self):
        p=self.file('long',b'x'*1025)
        with self.assertRaises(ValueError):list(lines(p,max_line=1024))
    def test_waveform_explicit_schema_and_validation(self):
        p=self.file('wave',b'0 80 2 3 4\r\n0 160 4 3 2\n')
        records=list(waveforms(p,samples=3));self.assertEqual(records[0].samples,(2,3,4))
        with self.assertRaises(ValueError):list(waveforms(p,samples=4))
        for data in (b'0 80 2 3 x',b'0 80 2 3 65536',b'0 80 2 3 1_0'):
            p.write_bytes(data)
            with self.assertRaises(ValueError):list(waveforms(p,samples=3))
    def test_sensor_na_and_reset_not_repaired(self):
        p=self.file('sensor.csv',b'time_us,ax\n100,2\n50,NA\n60,3\n')
        rows=list(sensors(p,['ax']));self.assertIsNone(rows[1][2]['ax']);self.assertEqual(rows[1][1],50)
    def test_bad_csv_header_and_shape(self):
        for data in (b'time_us,ax,ax\n1,2,3\n',b'time_us,ax\n1,2,3\n',b'time_us,ax\n1\n'):
            p=self.file('x',data)
            with self.assertRaises(ValueError):list(sensors(p,['ax']))
    def test_event_channels_never_merge(self):
        p=self.file('events',b'channel,timestamp\n0,80\n1,10\n0,160\n')
        self.assertEqual([e[1] for e in events(p)],['0','1','0'])
    def test_unknown_units_and_clock_evidence(self):
        self.assertIsNone(seconds_per_tick({}))
        with self.assertRaises(ValueError):seconds_per_tick({'seconds_per_tick':1})
        for value in ('NaN',-1,0):
            with self.assertRaises(ValueError):seconds_per_tick({'seconds_per_tick':value,'clock_evidence':'x'})
    def test_timing_reset_and_lattice(self):
        t=Timing()
        for x in (80,160,160,40,120,280):t.add(x)
        r=t.report();self.assertEqual(r['zero_gaps'],1);self.assertEqual(r['backward_gaps'],1)
        self.assertEqual(r['positive_gap_gcd_raw_units'],80);self.assertIsNone(r['positive_gap_rate_hz'])
    def test_rct_chunk_state_exact_boundary(self):
        h=Health(h=1,alpha=2**-4,window=512)
        for _ in range(h.rct_cutoff-1):self.assertTrue(h.step(0))
        self.assertFalse(h.step(0));self.assertTrue(h.report()['latched_failure'])
    def test_apt_exact_tail_and_minimal_cutoff(self):
        h=Health(h=0.5,alpha=2**-20,window=512)
        p=2**-0.5;c=h.apt_cutoff
        self.assertLessEqual(binom.sf(c-2,511,p),h.alpha)
        self.assertGreater(binom.sf(c-3,511,p),h.alpha)
    def test_health_reject_invalid_parameters(self):
        for h in (0,-1,2,float('nan')):
            with self.assertRaises(ValueError):Health(h=h)
    def test_predictable_balanced_stream_no_entropy_credit(self):
        r=screen(bytes([0,1])*1000)
        self.assertGreater(r['iid_only_99pct_marginal_min_entropy_lower_bits'],0.8)
        self.assertEqual(r['heldout_predictors_context_lengths']['1']['accuracy'],1)
        self.assertEqual(r['validated_entropy_credit_bits'],0)
    def test_stuck_stream_zero_iid_estimate(self):
        self.assertEqual(screen(bytes(2000))['iid_only_99pct_marginal_min_entropy_lower_bits'],0)
    def test_vn_known_vector_and_tails(self):
        p=self.file('bits',bytes([0,1,1,0])*8+bytes([1]))
        r=treatments(p,self.root/'treated')
        self.assertEqual((self.root/'treated/vn.bin').read_bytes(),b'UU')
        self.assertEqual(r['vn_unpaired_input_bits'],1)
        self.assertEqual(r['validated_entropy_credit_bits'],0)
    def test_conditioning_fixed_blocks_and_chunk_boundary(self):
        bits=bytes([0,1])*513
        p=self.file('bits',bits);r=treatments(p,self.root/'treated')
        self.assertEqual((self.root/'treated/raw.bin').read_bytes(),bytes([0x55])*128)
        expected=hashlib.sha256(b'RNGLab15 research\0'+bytes([0x55])*64).digest()*2
        self.assertEqual((self.root/'treated/sha256.bin').read_bytes(),expected)
        self.assertEqual(r['raw_tail_bits_discarded'],2)
    def test_toeplitz_budget_and_identity(self):
        # Diagonal seed gives the identity matrix for 4x4.
        bits=[1,0,1,1];seed=[0,0,0,1,0,0]
        self.assertEqual(toeplitz_reference(bits,seed,3,conditional_min_entropy=4,epsilon=0.8),bits[:3])
        with self.assertRaises(ValueError):toeplitz_reference(bits,seed,4,conditional_min_entropy=3,epsilon=0.01)
    def test_finite_source_no_padding_no_replay(self):
        s=FiniteSource(b'\x00\xff');self.assertEqual(s.random_float(),0.5/256)
        self.assertEqual(s.random_float(),255.5/256)
        with self.assertRaises(EOFError):s.read_bytes(1)
        self.assertEqual(s.position,2)
        self.assertEqual(s.stream_digest,hashlib.sha256(b'\x00\xff').hexdigest())
    def test_defects_reproducible_and_targeted(self):
        data=bytes(range(256))*4
        for name in ('stuck','bias','repeat','quantized','drift','periodic_dropout'):
            self.assertEqual(defects(data,name,1),defects(data,name,1))
            self.assertEqual(len(defects(data,name,1)),len(data))
        self.assertEqual(set(defects(data,'stuck',1)),{0})
    def test_exact_model_references(self):
        self.assertAlmostEqual(finite_bead_energy(128,2),0.5/math.tanh(1),delta=0.001)
        self.assertTrue(-2<ising_exact_energy(4,3)<0)
    def test_all_models_finite_slot(self):
        import random
        data=random.Random(44).randbytes(256)
        for model in ('pi','vmc','ising','pimc'):
            r=model_run(model,data,31)
            self.assertLessEqual(r['proposal_bytes_consumed'],256)
            self.assertTrue(math.isfinite(r['outcome']))
    def test_holm_family_correctly_adjusted(self):
        rows=[]
        for i in range(20):
            for name,fail in [('os',False),('stuck',True)]:
                rows.append(dict(identity=i,model='pi',treatment=name,cohort='confirmatory',tolerance_failure=fail,outcome=int(fail)))
        r=paired_summary(rows)[0];self.assertTrue(r['reject_fwer_0_05']);self.assertEqual(r['pairs'],20)
    def test_end_to_end_registry_and_privacy(self):
        source=self.file('data.csv',b'time_us,ax\n80,2\n160,3\n240,2\n320,3\n')
        cfg=self.root/'config.json'
        cfg.write_text(json.dumps({'schema_version':1,'sources':[{'path':str(source),'kind':'sensor_csv','columns':['ax'],'privacy_review':'non_sensitive_physical_measurement'}]}))
        r=run(cfg,self.root/'run',self.root/'vault');self.assertEqual(r['status'],'complete')
        self.assertEqual(r['validated_entropy_credit_bits'],0)
        for name,sha in r['output_hashes'].items():self.assertEqual(digest(self.root/'run'/name),sha)
        with self.assertRaises(FileExistsError):run(cfg,self.root/'run',self.root/'vault')
        value=json.loads(cfg.read_text());del value['sources'][0]['privacy_review'];cfg.write_text(json.dumps(value))
        with self.assertRaises(ValueError):run(cfg,self.root/'failed',self.root/'vault')
        self.assertEqual(json.loads((self.root/'failed/registry.json').read_text())['status'],'failed')
    def test_pipeline_rejects_input_mutation(self):
        source=self.file('data.csv',b'time_us,ax\n80,2\n160,3\n')
        cfg=self.file('config.json',json.dumps({'schema_version':1,'sources':[{'path':str(source),'kind':'sensor_csv','columns':['ax'],'privacy_review':'non_sensitive_physical_measurement'}]}).encode())
        from rnglab15 import pipeline
        original=pipeline.describe_sensors
        def mutate(*args):
            r=original(*args);source.write_bytes(b'changed');return r
        with patch.object(pipeline,'describe_sensors',mutate):
            with self.assertRaises(ValueError):run(cfg,self.root/'run',self.root/'vault')
        self.assertEqual(json.loads((self.root/'run/registry.json').read_text())['status'],'failed')

class ExternalAdapterTests(unittest.TestCase):
    def test_nist_small_input_is_explicit_blocker(self):
        from rnglab15.nist import assess
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=root/'bits';p.write_bytes(bytes([0,1])*100)
            r=assess(p,root/'nist');self.assertEqual(r['status'],'blocked')
            self.assertEqual(r['validated_entropy_credit_bits'],0)
    def test_nist_rejects_wrong_symbol_width(self):
        from rnglab15.nist import assess
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=root/'bits';p.write_bytes(bytes([255])*1_000_000)
            with patch('rnglab15.nist.shutil.which',return_value='/bin/true'):
                with self.assertRaises(ValueError):assess(p,root/'nist')
    def test_nist_subprocess_is_bounded_and_logged(self):
        from rnglab15.nist import assess
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=root/'bits';p.write_bytes(bytes([0,1])*500_000)
            with patch('rnglab15.nist.shutil.which',return_value='/bin/true'),patch('rnglab15.nist.subprocess.run') as runner:
                runner.return_value.returncode=0
                r=assess(p,root/'nist',timeout=30)
                self.assertEqual(r['status'],'executed')
                self.assertEqual(runner.call_args.kwargs['timeout'],30)
                self.assertEqual(r['command'][1:3],['-i','-a'])
    def test_campaign_replays_exact_archived_os_fixture(self):
        from rnglab15.applications import simulation_campaign
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            first=simulation_campaign(root/'a',trials=2)
            second=simulation_campaign(root/'b',trials=2,reference_path=root/'a/os-reference.bin')
            self.assertEqual(first,second)
    def test_short_waveform_schema_runs_without_long_window_assumptions(self):
        from rnglab15.pipeline import describe_waveform
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=root/'wave';p.write_text('0 80 '+' '.join(['3300']*34+['65535']*8+['1800']*22)+'\n')
            out=root/'out';out.mkdir()
            r=describe_waveform(p,{'samples':64,'exploratory_gap_divisor':80},out)
            self.assertEqual(r['channels']['0']['waveform_samples'],64)
            self.assertEqual(r['channels']['0']['clipped_event_fraction'],1)
    def test_health_alarm_stays_latched_after_recovery(self):
        h=Health(h=1,alpha=0.25)
        for x in [0,0,0]:h.step(x)
        self.assertFalse(h.step(1));self.assertFalse(h.step(0))
    def test_finite_bit_reads_have_exact_consumption_digest(self):
        s=FiniteSource(b'\xa5');self.assertEqual(s.read_bits(4),[1,0,1,0])
        self.assertEqual(s.stream_digest,hashlib.sha256(b'\xa5').hexdigest())
    def test_code_provenance_without_repository_layout(self):
        from rnglab15.pipeline import code_manifest,git_parent
        with tempfile.TemporaryDirectory() as tmp:
            result=code_manifest(Path(tmp))
            self.assertIn('rnglab15/pipeline.py',result)
            self.assertIsNone(git_parent(Path(tmp)))
