"""Add workloads by registering a dimension, known target, and contributions."""
from dataclasses import dataclass
from typing import Callable
import numpy as np


@dataclass(frozen=True)
class Workload:
    name: str
    dimension: int
    target: float
    contributions: Callable


WORKLOADS = {
    "pi": Workload("pi", 2, float(np.pi), lambda u: 4.0 * (np.sum(u*u, axis=1) <= 1)),
    "integral": Workload("integral", 1, 1/3, lambda u: u[:, 0]**2),
}
