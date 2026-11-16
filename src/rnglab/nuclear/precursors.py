"""RNGLab v0.20 nuclear module, milestone B2: synthetic delayed-neutron precursor groups.

A precursor group is a synthetic nuclide that decays by beta-delayed neutron
emission with probability P_n and by ordinary beta- decay otherwise, both to
stable daughters. Only the beta-delayed-neutron channel emits a neutron. The
groups are abstract: half-lives and P_n are chosen by the user, and Z starts at
170, outside the chart of known nuclides, so no group stands for real
delayed-neutron data. The same groups can be simulated nucleus by nucleus
(Layer A) or as populations (B2); the aggregate precursor kinetics of Layer B
will be checked against them.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .decay import Channel, DecayScheme, Mode, Nuclide


@dataclass(frozen=True)
class PrecursorGroup:
    name: str
    half_life: float
    neutron_probability: float

    def __post_init__(self):
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("group name must be a nonempty string")
        if not (isinstance(self.half_life, (int, float)) and math.isfinite(self.half_life) and self.half_life > 0):
            raise ValueError("half-life must be finite and positive")
        if not (isinstance(self.neutron_probability, (int, float)) and 0.0 < self.neutron_probability <= 1.0):
            raise ValueError("neutron emission probability must lie in (0, 1]")


def precursor_scheme(groups: Sequence[PrecursorGroup], *, name: str = "precursor-groups",
                     base_z: int = 170) -> DecayScheme:
    """A decay scheme with one precursor and its stable daughters per group."""
    if not groups:
        raise ValueError("need at least one group")
    nuclides = []
    for index, group in enumerate(groups):
        z = base_z + 3 * index
        a = 2 * z + 100
        channels = [Channel(f"{group.name}-n", Mode.BETA_DELAYED_NEUTRON, group.neutron_probability)]
        if group.neutron_probability < 1.0:
            channels.append(Channel(f"{group.name}-b", Mode.BETA_MINUS, 1.0 - group.neutron_probability))
        nuclides.append(Nuclide(group.name, z, a, group.half_life, tuple(channels)))
        nuclides.append(Nuclide(f"{group.name}-n", z + 1, a - 1, math.inf))
        if group.neutron_probability < 1.0:
            nuclides.append(Nuclide(f"{group.name}-b", z + 1, a, math.inf))
    return DecayScheme(name, nuclides)
