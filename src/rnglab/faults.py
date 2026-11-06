"""RNGLab v0.18 namespace: controlled fault adapters.

Re-exports DegradedFaultSource from rnglab14 unchanged (bias, Markov rho on
contiguous emitted bits, stuck bit, periodic coupling). Its rho is defined on
the bit stream, not on consumed uniforms.

The Phase II channel faults (acceptance-uniform shift, proposal/acceptance
coupling) are implemented only in the frozen standalone study
research/phase2-20261030/experiment/fault_study.py and are not importable from
here. The Phase 3 severity adapters will be added to this module as new v0.18
code and validated against those definitions.
"""
from rnglab14.sources.faults import DegradedFaultSource

__all__ = ["DegradedFaultSource"]
