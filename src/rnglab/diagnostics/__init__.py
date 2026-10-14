"""Descriptive checks and p-values calibrated ONLY under the IID null."""
import numpy as np
from scipy.stats import binomtest, chisquare, norm


def diagnose(bits, uniforms, lags=(1, 2, 8, 32), bins=16):
    bits = np.asarray(bits, dtype=np.uint8)
    n = len(bits)
    if n < 2 or not len(uniforms):
        raise ValueError("diagnostics require at least two bits and one uniform")
    ones = int(bits.sum())
    result = {"bit_count": n, "ones": ones, "bit_bias": ones/n - 0.5,
              "monobit_z": (2*ones-n)/np.sqrt(n),
              "monobit_p_iid": float(binomtest(ones, n, 0.5).pvalue)}
    for lag in lags:
        a, b = bits[:-lag].astype(float), bits[lag:].astype(float)
        if len(a) < 2 or np.std(a) == 0 or np.std(b) == 0:
            rho, p = None, None
        else:
            rho = float(np.corrcoef(a, b)[0, 1])
            p = float(2 * norm.sf(abs(rho) * np.sqrt(len(a))))
        result[f"lag{lag}_rho"] = rho
        result[f"lag{lag}_p_iid_approx"] = p
    hist, _ = np.histogram(uniforms, bins=bins, range=(0, 1))
    result["uniform_chi2_p_iid"] = float(chisquare(hist).pvalue) if len(uniforms)/bins >= 5 else None
    return result, hist
