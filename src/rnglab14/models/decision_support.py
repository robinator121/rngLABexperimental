"""RNGLab v0.14: Monte Carlo Decision-Support System (DSS) Simulator.

Studies how randomness-source degradation alters stochastic decision-making:
- Option ranking inversion
- Value-at-Risk (VaR) and Expected Shortfall misestimation
- Decision confidence / regret under corrupted distributions
- Probability of making an incorrect operational decision (P_error)
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Dict, List, Optional, Sequence, Union

import numpy as np

from rnglab14.sources.base import RandomSource, StochasticChannels


@dataclass(frozen=True)
class DecisionOption:
    name: str
    true_mean_payoff: float
    true_volatility: float
    skewness: float = 0.0


@dataclass(frozen=True)
class DecisionEvaluationResult:
    total_scenarios: int
    option_names: List[str]
    estimated_means: Dict[str, float]
    estimated_var_95: Dict[str, float]
    selected_option: str
    optimal_option: str
    ranking_correct: bool
    decision_regret: float
    confidence_spread: float
    provenance_hash: str


class MonteCarloDecisionSupport:
    """Stochastic Decision Support System evaluating alternative strategic investments or policies under uncertainty."""

    def __init__(self, options: Optional[Sequence[DecisionOption]] = None):
        if options is None:
            # Canonical 3-option choice problem:
            # Option A: Conservative low risk, moderate return (mean=100, vol=15)
            # Option B: High return, high volatility (mean=120, vol=40) -> Optimal choice under risk neutrality
            # Option C: High risk, moderate return (mean=105, vol=35)
            self.options = [
                DecisionOption("Conservative_A", true_mean_payoff=100.0, true_volatility=15.0),
                DecisionOption("Growth_B", true_mean_payoff=120.0, true_volatility=40.0),
                DecisionOption("Speculative_C", true_mean_payoff=105.0, true_volatility=35.0),
            ]
        else:
            self.options = list(options)

        # True best option under risk neutrality
        self.optimal_option = max(self.options, key=lambda opt: opt.true_mean_payoff).name

    def evaluate_decisions(
        self,
        source_or_channels: Union[RandomSource, StochasticChannels],
        scenarios: int = 5000,
    ) -> DecisionEvaluationResult:
        if isinstance(source_or_channels, StochasticChannels):
            src = source_or_channels.proposal
        else:
            src = source_or_channels

        option_payoffs: Dict[str, List[float]] = {opt.name: [] for opt in self.options}

        for _ in range(scenarios):
            # Box-Muller normal generation using source
            for opt in self.options:
                u1 = max(1e-12, src.random_float())
                u2 = src.random_float()
                z = math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)
                payoff = opt.true_mean_payoff + z * opt.true_volatility
                option_payoffs[opt.name].append(payoff)

        estimated_means = {}
        estimated_var95 = {}

        for opt in self.options:
            arr = np.array(option_payoffs[opt.name])
            m = float(np.mean(arr))
            v95 = float(np.percentile(arr, 5))  # 5th percentile (95% Value at Risk)
            estimated_means[opt.name] = round(m, 4)
            estimated_var95[opt.name] = round(v95, 4)

        # Decision rule: pick option with highest estimated mean payoff
        selected = max(estimated_means.keys(), key=lambda k: estimated_means[k])
        true_best_val = next(o.true_mean_payoff for o in self.options if o.name == self.optimal_option)
        true_selected_val = next(o.true_mean_payoff for o in self.options if o.name == selected)
        regret = max(0.0, true_best_val - true_selected_val)

        # Confidence spread: gap between 1st and 2nd option
        sorted_means = sorted(estimated_means.values(), reverse=True)
        conf_spread = sorted_means[0] - sorted_means[1] if len(sorted_means) > 1 else 0.0

        return DecisionEvaluationResult(
            total_scenarios=scenarios,
            option_names=[o.name for o in self.options],
            estimated_means=estimated_means,
            estimated_var_95=estimated_var95,
            selected_option=selected,
            optimal_option=self.optimal_option,
            ranking_correct=(selected == self.optimal_option),
            decision_regret=round(regret, 4),
            confidence_spread=round(conf_spread, 4),
            provenance_hash=src.stream_digest,
        )
