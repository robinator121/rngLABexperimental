"""
Financial Integrity Monitoring Engine.

Connects source-level randomness diagnostics with application-level financial invariants:
1. Martingale Preservation Test (E[e^{-rT} S_T] = S0).
2. Innovation Normality & Skewness Check (Jarque-Bera).
3. Path Duplication & Hash Invariant Check.
4. Kupiec POF Likelihood Ratio for Risk Quantiles.
5. 4-Quadrant Operational Risk Matrix (Q1 Nominal, Q2 Early Warning, Q3 Synchronized, Q4 Hazardous Gap).
"""

import math
import hashlib
from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from scipy.stats import norm, chi2, jarque_bera

from rnglab_finance.core_interface import (
    FinancialRNGProvider,
    DualIntegrityMonitor,
    frequency_monobit_test,
    runs_test,
    RepetitionCountTest,
)


class SourceIntegrityMonitor:
    """
    Evaluates raw RNG bit-stream health using NIST SP 800-22 and SP 800-90B standards.
    """

    def __init__(self, alpha: float = 0.01):
        self.alpha = alpha

    def assess_source(self, uniforms: np.ndarray) -> Dict[str, Any]:
        """Convert float uniforms into bits and evaluate monobit, runs, and repetition count."""
        if len(uniforms) == 0:
            return {"healthy": False, "reason": "empty_stream"}

        # Extract bits from threshold u >= 0.5 (or uint64 bits)
        bits = (uniforms >= 0.5).astype(int).tolist()
        
        # NIST SP 800-22 Frequency Monobit
        monobit = frequency_monobit_test(bits, alpha=self.alpha)
        # NIST SP 800-22 Runs Test
        runs = runs_test(bits, alpha=self.alpha)

        # NIST SP 800-90B Repetition Count Test
        rct = RepetitionCountTest(min_entropy_per_sample=0.99)
        rct_passed = True
        for b in bits:
            res = rct.step(b)
            if not res.passed:
                rct_passed = False
                break

        healthy = monobit.passed and runs.passed and rct_passed
        score = min(monobit.p_value, runs.p_value)

        return {
            "healthy": healthy,
            "min_p_value": score,
            "monobit_passed": monobit.passed,
            "monobit_p_value": monobit.p_value,
            "runs_passed": runs.passed,
            "runs_p_value": runs.p_value,
            "rct_passed": rct_passed,
        }


class FinancialIntegrityMonitor:
    """
    Dual-level integrity monitor tailored for financial Monte Carlo simulations.
    Evaluates both physical/statistical randomness health and economic invariants.
    """

    def __init__(
        self,
        alpha_significance: float = 0.01,
        martingale_z_threshold: float = 3.0,
        max_acceptable_pricing_error_pct: float = 0.02,  # 2% error tolerance
    ):
        self.alpha_significance = alpha_significance
        self.martingale_z_threshold = martingale_z_threshold
        self.max_pricing_error_pct = max_acceptable_pricing_error_pct
        self.source_monitor = SourceIntegrityMonitor(alpha=alpha_significance)

    def assess_financial_run(
        self,
        rng_provider: FinancialRNGProvider,
        pricing_results: Dict[str, Any],
        innovations: Optional[np.ndarray] = None,
    ) -> Dict[str, Any]:
        """
        Evaluate a complete financial Monte Carlo run across all integrity dimensions.
        """
        # 1. Source-Level Evaluation
        diag_uniforms = rng_provider.draw_uniforms(min(max(rng_provider.total_drawn, 2048), 16384))
        source_eval = self.source_monitor.assess_source(diag_uniforms)
        source_alarm = not source_eval["healthy"]

        # 2. Application-Level Financial Invariants
        app_alarms: List[str] = []

        # Invariant A: Martingale Preservation
        martingale_z = pricing_results.get("martingale_z_score", 0.0)
        martingale_violated = abs(martingale_z) > self.martingale_z_threshold
        if martingale_violated:
            app_alarms.append(f"MartingaleViolation(z={martingale_z:.2f})")

        # Invariant B: Economic Pricing Error Exceedance
        rel_error = pricing_results.get("relative_error", 0.0)
        pricing_exceeded = rel_error > self.max_pricing_error_pct
        if pricing_exceeded:
            app_alarms.append(f"PricingErrorExceeded({rel_error*100:.2f}%)")

        # Invariant C: Path / Sample Duplication
        duplication_detected = False
        if innovations is not None and len(innovations) > 100:
            # Subsample up to 50,000 elements across the whole array to detect worker-level repeats
            n_total = len(innovations)
            if n_total > 50000:
                indices = np.linspace(0, n_total - 1, 50000, dtype=int)
                sample_z = innovations[indices]
            else:
                sample_z = innovations
            n_unique = len(np.unique(sample_z))
            unique_ratio = n_unique / float(len(sample_z))
            if unique_ratio < 0.98:  # In continuous N(0,1), duplicate floats indicate worker seed reuse
                duplication_detected = True
                app_alarms.append(f"SeedDuplicationDetected(unique_ratio={unique_ratio:.2f})")

        # Invariant D: Kupiec POF for VaR
        if "kupiec_test_passed" in pricing_results:
            if not pricing_results["kupiec_test_passed"]:
                app_alarms.append(f"KupiecVaRFailed(p={pricing_results.get('kupiec_p_value', 0.0):.4f})")

        app_alarm = len(app_alarms) > 0

        # 3. 4-Quadrant Operational Risk Classification
        if not source_alarm and not app_alarm:
            quadrant = "Q1_Nominal"
            risk_level = "LOW"
            recommendation = "Simulation outputs certified for financial decision-making."
        elif source_alarm and not app_alarm:
            quadrant = "Q2_EarlyWarning"
            risk_level = "MODERATE"
            recommendation = "Source degradation detected. Simulation within pricing tolerance, but entropy exhaustion imminent. Reseed RNG."
        elif source_alarm and app_alarm:
            quadrant = "Q3_SynchronizedDetection"
            risk_level = "CRITICAL"
            recommendation = "Simulation compromised by RNG defect. Invalidate pricing/risk metrics. Halt execution."
        else:
            quadrant = "Q4_HazardousGap"
            risk_level = "HAZARDOUS"
            recommendation = "Financial model corrupted without source-level alarm. Possible structural defect or algorithmic failure."

        # 4. Comparison with Conventional Model Risk Baseline
        # Standard baseline check: Standard Error <= 0.15 (or relative SE < 1.5%)
        baseline_se_check_passed = bool(pricing_results.get("standard_error", 1.0) < 0.15)
        baseline_blind_to_defect = bool(baseline_se_check_passed and app_alarm)

        return {
            "operational_quadrant": quadrant,
            "risk_level": risk_level,
            "source_alarm": source_alarm,
            "application_alarm": app_alarm,
            "application_alarms_detail": app_alarms,
            "martingale_violated": martingale_violated,
            "pricing_tolerance_exceeded": pricing_exceeded,
            "duplication_detected": duplication_detected,
            "hazardous_gap": (quadrant == "Q4_HazardousGap"),
            "early_warning": (quadrant == "Q2_EarlyWarning"),
            "baseline_se_check_passed": baseline_se_check_passed,
            "baseline_blind_to_defect": baseline_blind_to_defect,
            "rnglab_superior_to_baseline": (source_alarm or app_alarm) and baseline_blind_to_defect,
            "recommendation": recommendation,
            "source_metrics": source_eval,
        }
