"""RNGLab v0.14: Path Integral Monte Carlo (PIMC) Quantum Harmonic Oscillator.

Discretizes imaginary time beta = 1/T into M time-slices ('beads') on a ring polymer.
Isolates:
- Proposal Channel: Local bead displacement proposal x_j -> x_j + delta x
- Acceptance Channel: Metropolis acceptance ratio comparison u in [0, 1)
- Initialization Channel: Initial bead positions

Calculates:
- Thermodynamic energy estimator
- Virial energy estimator (low-variance quantum estimator)
- Quantum bead spread (Radius of gyration R_g^2)
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
import struct
from typing import Any, Mapping, Optional, Union

import numpy as np
from scipy.stats import t as student_t

from rnglab14.sources.base import RandomSource, StochasticChannels
from rnglab12.uncertainty import summarize_series


def exact_finite_bead_energy(beads: int, beta: float) -> float:
    """Primitive-action harmonic reference; m=hbar=omega=1.

    Z_M = [2 sinh(M asinh(beta/(2M)))]**-1. Its logarithmic beta
    derivative gives this reference, including finite-bead discretization.
    """
    if type(beads) is not int or beads < 1 or not math.isfinite(beta) or beta <= 0:
        raise ValueError("positive integer beads and finite positive beta required")
    z = beta / (2.0 * beads)
    return 0.5 / (math.tanh(beads * math.asinh(z)) * math.sqrt(1.0 + z*z))


def harmonic_energy_estimators(path: np.ndarray, beta: float) -> tuple[np.ndarray, np.ndarray]:
    """Return thermodynamic and centroid-virial energies, last axis=beads.

    Differentiating S=M*sum(dx**2)/(2*beta)+beta*mean(V) gives
    E_T=M/(2*beta)-M*sum(dx**2)/(2*beta**2)+mean(V).
    The beta normalization prefactor contributes M/(2*beta).
    """
    x = np.asarray(path, dtype=float)
    if x.ndim < 1 or x.shape[-1] < 1 or not np.isfinite(x).all():
        raise ValueError("finite path with a nonempty bead axis required")
    if not math.isfinite(beta) or beta <= 0:
        raise ValueError("finite positive beta required")
    m = x.shape[-1]
    springs = np.sum((np.roll(x, -1, axis=-1) - x)**2, axis=-1)
    x2 = np.mean(x*x, axis=-1)
    centroid = np.mean(x, axis=-1)
    thermo = m/(2.0*beta) - m*springs/(2.0*beta**2) + 0.5*x2
    virial = 1.0/(2.0*beta) + x2 - 0.5*centroid**2
    return thermo, virial


def harmonic_centroid_heatbath(path: np.ndarray, beta: float,
                              source: RandomSource) -> np.ndarray:
    """Return a path with a refreshed harmonic centroid; do not mutate input.

    Writing x_j = c + y_j, sum(y_j)=0, the primitive harmonic action factors
    as S(y) + beta*c**2/2. Thus c ~ Normal(0, 1/beta) is the exact conditional
    law at every bead count. Internal coordinates and springs are preserved
    up to floating-point rounding; no Metropolis acceptance draw is needed.

    Consume exactly 16 bytes from the supplied source: two big-endian 64-bit
    words map to 53-bit uniforms u,v in [0,1), then use the cosine Box-Muller
    variate sqrt(-2*log1p(-u))*cos(2*pi*v)/sqrt(beta). The sine variate is
    discarded. Correct conditional sampling assumes IID uniform source words
    independent of the retained internal coordinates. This interface does not
    verify those assumptions, certify entropy, or establish full-path mixing.
    """
    x = np.asarray(path, dtype=float)
    if x.ndim != 1 or x.size < 1 or not np.isfinite(x).all():
        raise ValueError("finite one-dimensional nonempty path required")
    if not math.isfinite(beta) or beta <= 0:
        raise ValueError("finite positive beta required")
    if not isinstance(source, RandomSource):
        raise TypeError("explicit RandomSource required")
    # Scale before summing to avoid overflow in the mean of large finite paths.
    centroid = math.fsum(float(value)/x.size for value in x)
    with np.errstate(over="ignore", invalid="ignore"):
        internal = x - centroid
    if not np.isfinite(internal).all():
        raise ValueError("centered path must remain finite")
    a, b = struct.unpack(">QQ", source.read_bytes(16))
    u, v = (a >> 11)*2.0**-53, (b >> 11)*2.0**-53
    new_centroid = (math.sqrt(-2.0*math.log1p(-u))
                    * math.cos(2.0*math.pi*v)/math.sqrt(beta))
    with np.errstate(over="ignore", invalid="ignore"):
        updated = internal + new_centroid
    if not np.isfinite(updated).all():
        raise ValueError("updated path must remain finite")
    return updated


def sample_exact_finite_bead_paths(source: RandomSource, beads: int, beta: float,
                                 independent_paths: int) -> np.ndarray:
    """Exact Gaussian finite-bead reference conditional on IID source words.

    Bulk bytes -> 53-bit uniforms -> Box-Muller -> harmonic precision solve.
    This is a distinct reference sampler, not a repair of Metropolis mixing.
    It consumes only the explicitly supplied source; independent channels are
    neither required nor inferred. A corrupted source can invalidate the law.
    """
    exact_finite_bead_energy(beads, beta)  # validate the mathematical domain
    if type(independent_paths) is not int or independent_paths < 2:
        raise ValueError("at least two reference paths required")
    identity = np.eye(beads)
    difference = np.roll(identity, -1, axis=0) - identity
    precision = beta/beads*identity + beads/beta*(difference.T @ difference)
    covariance = np.linalg.inv(precision)
    cholesky = np.linalg.cholesky(covariance)
    count = independent_paths*beads
    pairs = (count+1)//2
    words = np.frombuffer(source.read_bytes(pairs*16), dtype=">u8")
    uniforms = (words >> np.uint64(11)).astype(float).reshape(pairs,2)*2.0**-53
    radius = np.sqrt(-2.0*np.log1p(-uniforms[:,0]))
    angles = 2.0*math.pi*uniforms[:,1]
    normals = np.column_stack((radius*np.cos(angles),radius*np.sin(angles))).ravel()[:count]
    return normals.reshape(independent_paths,beads) @ cholesky.T


def summarize_independent_reference(values: np.ndarray, reference: float,
                                    confidence: float = .95) -> dict:
    """Student-t approximation for independently sampled path estimators.

    Gaussian paths do not make quadratic energy observables Gaussian. The t
    interval is asymptotic, so held-out finite-sample calibration is required.
    """
    a = np.asarray(values,dtype=float)
    if a.ndim != 1 or len(a)<2 or not np.isfinite(a).all():
        raise ValueError("finite one-dimensional independent estimator sample required")
    if not math.isfinite(confidence) or not 0<confidence<1 or not math.isfinite(reference):
        raise ValueError("finite reference and confidence strictly between zero and one required")
    mean=float(a.mean()); std=float(a.std(ddof=1)); se=std/math.sqrt(len(a))
    critical=float(student_t.ppf((1+confidence)/2,len(a)-1))
    bounds=[mean-critical*se,mean+critical*se]
    if std==0:
        return {"mean":mean,"sample_std":std,"n":len(a),"status":"ABSTAIN_CONSTANT_REFERENCE_SAMPLE",
                "mcse":None,"interval":None,"covered":None,"reference":reference}
    return {"mean":mean,"sample_std":std,"n":len(a),"status":"IID_REFERENCE_ASYMPTOTIC_T",
        "confidence":confidence,"mcse":se,"interval":bounds,
        "covered":bool(bounds[0]<=reference<=bounds[1]),"reference":reference,
        "assumptions":"IID healthy source words and correct Gaussian transform; finite-sample t coverage not guaranteed"}


@dataclass(frozen=True)
class PIMCResult:
    beads_M: int
    temperature_T: float
    beta: float
    total_sweeps: int
    burn_in: int
    acceptance_rate: float
    virial_energy_mean: float
    virial_energy_std: float
    thermodynamic_energy_mean: float
    exact_quantum_energy: float
    energy_bias: float
    radius_of_gyration: float
    provenance_hash: str
    exact_finite_bead_energy: float
    discretization_bias: float
    virial_uncertainty: dict
    thermodynamic_uncertainty: dict
    sampler_protocol: str = "LOCAL_METROPOLIS_V1"
    centroid_update: Optional[dict] = None
    initialization_protocol: Optional[dict] = None
    source_provenance: Optional[dict] = None
    trace: Optional[dict] = None


class PathIntegralHarmonicOscillator:
    """Path Integral Monte Carlo (PIMC) for 1D Harmonic Oscillator (m=1, hbar=1, omega=1)."""

    def __init__(
        self,
        beads: int = 20,
        temperature: float = 0.5,
        displacement_delta: float = 0.5,
    ):
        if type(beads) is not int or beads < 2:
            raise ValueError("sampler requires integer beads >= 2")
        if not math.isfinite(temperature) or temperature <= 0:
            raise ValueError("temperature must be finite and positive")
        if not math.isfinite(displacement_delta) or displacement_delta <= 0:
            raise ValueError("displacement_delta must be finite and positive")
        self.M = beads
        self.T = temperature
        self.beta = 1.0 / temperature
        self.dtau = self.beta / self.M
        self.delta = displacement_delta
        # Exact quantum thermodynamic energy: E_exact = 1/2 * coth(beta / 2)
        self.exact_energy = 0.5 * (1.0 / math.tanh(self.beta / 2.0))
        self.finite_bead_energy = exact_finite_bead_energy(self.M, self.beta)
        self.path = np.zeros(self.M, dtype=float)

    def run_exact_reference(self, source: RandomSource, independent_paths: int = 2048,
                            confidence: float = .95) -> dict:
        """Separate exact harmonic reference protocol with independent path draws."""
        paths=sample_exact_finite_bead_paths(source,self.M,self.beta,independent_paths)
        td,cv=harmonic_energy_estimators(paths,self.beta)
        return {"protocol":"EXACT_FINITE_BEAD_GAUSSIAN_REFERENCE_V1",
            "beads":self.M,"beta":self.beta,"independent_paths":independent_paths,
            "finite_bead_reference":self.finite_bead_energy,"continuum_reference":self.exact_energy,
            "discretization_bias":self.finite_bead_energy-self.exact_energy,
            "thermodynamic":summarize_independent_reference(td,self.finite_bead_energy,confidence),
            "centroid_virial":summarize_independent_reference(cv,self.finite_bead_energy,confidence),
            "consumed_source_digest":source.stream_digest,"source_manifest":source.manifest.to_dict(),
            "channels_used":["supplied source only"],"independence_verified":False,
            "scope":"Known harmonic Gaussian reference under IID uniform input assumption; not evidence Metropolis has mixed"}

    def run_simulation(
        self,
        source_or_channels: Union[RandomSource, StochasticChannels],
        sweeps: int = 2000,
        burn_in: int = 500,
        *,
        centroid_source: Optional[RandomSource] = None,
        initialization: str = "uniform",
        initial_path: Optional[np.ndarray] = None,
        initial_path_provenance: Optional[Mapping[str, Any]] = None,
        uncertainty_batch_sizes: tuple[int, ...] = (32, 128),
        record_trace: bool = False,
    ) -> PIMCResult:
        """Run local Metropolis, optionally followed by one centroid Gibbs update.

        The default preserves the existing transition and source consumption.
        A supplied centroid_source is used after each local sweep, including
        burn-in. Its separate argument does not prove independent origin; it
        may alias a local channel, and corrupted/dependent input can invalidate
        the conditional Gaussian update. This option is harmonic-model specific.

        An explicit initial_path is copied and takes precedence over the selected
        initialization mode, without consuming initialization-source bytes. Exact
        Gaussian initialization uses the existing two-path reference sampler and
        discards its second path. All optional arguments are validated before
        source consumption. Retained full paths/observables are stored only when
        record_trace=True; legacy scalar result fields retain their rounding.

        Every mode records the SHA-256 of the starting path (little-endian
        float64 bytes, bead order). initial_path_provenance, allowed only with
        initial_path, is a JSON-serialisable mapping stored as its JSON round
        trip; it is recorded as given, not verified.
        """
        if (type(sweeps) is not int or type(burn_in) is not int or burn_in < 0
                or sweeps - burn_in < 32):
            raise ValueError("integer sweeps/burn_in and at least 32 retained sweeps required")
        if centroid_source is not None and not isinstance(centroid_source, RandomSource):
            raise TypeError("centroid_source must be an explicit RandomSource or None")
        if not isinstance(initialization, str) or initialization not in ("uniform", "exact_gaussian"):
            raise ValueError("initialization must be 'uniform' or 'exact_gaussian'")
        if type(record_trace) is not bool:
            raise ValueError("record_trace must be a bool")
        try:
            batch_sizes = tuple(uncertainty_batch_sizes)
        except TypeError as error:
            raise ValueError("unique integer uncertainty batch sizes >= 2 required") from error
        if (not batch_sizes or any(type(size) is not int or size < 2 for size in batch_sizes)
                or len(set(batch_sizes)) != len(batch_sizes)):
            raise ValueError("unique integer uncertainty batch sizes >= 2 required")
        supplied_path = None
        if initial_path is not None:
            supplied_path = np.array(initial_path, dtype=float, copy=True)
            if supplied_path.shape != (self.M,) or not np.isfinite(supplied_path).all():
                raise ValueError("initial_path must be a finite one-dimensional path with M beads")
        path_provenance = None
        if initial_path_provenance is not None:
            if supplied_path is None:
                raise ValueError("initial_path_provenance requires initial_path")
            if not isinstance(initial_path_provenance, Mapping) or not all(
                    isinstance(key, str) for key in initial_path_provenance):
                raise ValueError("initial_path_provenance must be a mapping with string keys")
            try:
                path_provenance = json.loads(json.dumps(dict(initial_path_provenance), allow_nan=False))
            except (TypeError, ValueError) as error:
                raise ValueError("initial_path_provenance must be JSON-serialisable without NaN or infinity") from error
        if isinstance(source_or_channels, StochasticChannels):
            prop_src = source_or_channels.proposal
            acc_src = source_or_channels.acceptance
            init_src = source_or_channels.initialization
        else:
            prop_src = source_or_channels
            acc_src = source_or_channels
            init_src = source_or_channels

        role_sources = {"proposal": prop_src, "acceptance": acc_src,
                        "initialization": init_src,
                        "resampling": (source_or_channels.resampling
                                       if isinstance(source_or_channels, StochasticChannels) else None),
                        "centroid": centroid_source}
        if any(source is not None and not isinstance(source, RandomSource)
               for source in role_sources.values()) or any(
                   role_sources[role] is None for role in ("proposal", "acceptance", "initialization")):
            raise TypeError("local channels and optional role sources must be RandomSource instances")
        provenance_before = {role: {"consumed_bytes_before": source.consumed_bytes,
            "stream_digest_before": source.stream_digest,
            "manifest": {key: value for key, value in source.manifest.to_dict().items()
                         if key != "created_at"}}
            for role, source in role_sources.items() if source is not None}
        centroid_before = centroid_source.consumed_bytes if centroid_source is not None else None
        centroid_digest_before = centroid_source.stream_digest if centroid_source is not None else None
        M = self.M
        dtau = self.dtau
        delta = self.delta

        # Spring constant between adjacent beads: k_spring = M / (beta * hbar^2) = 1 / dtau
        k_spring = 1.0 / dtau

        init_before = init_src.consumed_bytes
        if supplied_path is not None:
            self.path[:] = supplied_path
            initialization_mode = "explicit_path"
        elif initialization == "exact_gaussian":
            self.path[:] = sample_exact_finite_bead_paths(init_src, M, self.beta, 2)[0]
            initialization_mode = "exact_gaussian"
        else:
            for j in range(M):
                self.path[j] = init_src.read_uniform(-0.5, 0.5)
            initialization_mode = "uniform"
        initialization_protocol = {"mode": initialization_mode,
            "requested_initialization": initialization,
            "source_bytes_consumed": init_src.consumed_bytes-init_before,
            "reference_paths_drawn": 2 if initialization_mode == "exact_gaussian" else 0,
            "reference_paths_discarded": 1 if initialization_mode == "exact_gaussian" else 0,
            "initial_path_sha256": hashlib.sha256(self.path.astype("<f8").tobytes()).hexdigest(),
            "initial_path_digest_scheme": "SHA-256 of little-endian float64 bytes in bead order",
            "initial_path_provenance": path_provenance,
            "independence_verified": False,
            "assumptions": ("Exact Gaussian initialization assumes IID uniform source words; "
                            "explicit paths carry no stationarity claim; "
                            "supplied provenance is recorded as given, not verified")}
        retained = sweeps-burn_in
        trace_initial = self.path.copy() if record_trace else None
        trace_paths = np.empty((retained, M)) if record_trace else None
        trace_centroid = np.empty(retained) if record_trace else None
        trace_rg2 = np.empty(retained) if record_trace else None

        virial_energies = []
        thermo_energies = []
        rg_vals = []
        accepted = 0
        total_proposals = 0

        for sweep in range(sweeps):
            # Sweep through all beads
            for j in range(M):
                j_prev = (j - 1) % M
                j_next = (j + 1) % M
                x_old = self.path[j]

                # 1. Coordinate proposal
                u_prop = prop_src.random_float()
                dx = (u_prop * 2.0 - 1.0) * delta
                x_new = x_old + dx
                total_proposals += 1

                # Change in kinetic spring energy:
                # E_kin = 1/2 * k_spring * [(x_j - x_prev)^2 + (x_next - x_j)^2]
                dE_spring = 0.5 * k_spring * (
                    (x_new - self.path[j_prev]) ** 2 + (self.path[j_next] - x_new) ** 2
                    - (x_old - self.path[j_prev]) ** 2 - (self.path[j_next] - x_old) ** 2
                )
                # Change in potential energy: V(x) = 1/2 * x^2
                # Action contribution: S_pot = dtau * V(x)
                dS_pot = dtau * 0.5 * (x_new ** 2 - x_old ** 2)

                dS_total = dE_spring + dS_pot

                # 2. Acceptance decision
                if dS_total <= 0.0:
                    self.path[j] = x_new
                    accepted += 1
                else:
                    u_acc = acc_src.random_float()
                    if u_acc < math.exp(-dS_total):
                        self.path[j] = x_new
                        accepted += 1

            if centroid_source is not None:
                self.path[:] = harmonic_centroid_heatbath(self.path, self.beta, centroid_source)

            if sweep >= burn_in:
                e_thermo, e_virial = harmonic_energy_estimators(self.path, self.beta)
                x_cm = float(np.mean(self.path))
                rg2 = float(np.mean((self.path - x_cm) ** 2))

                virial_energies.append(e_virial)
                thermo_energies.append(e_thermo)
                rg_vals.append(math.sqrt(max(0.0, rg2)))
                if record_trace:
                    index = sweep-burn_in
                    trace_paths[index] = self.path
                    trace_centroid[index] = x_cm
                    trace_rg2[index] = rg2

        v_arr = np.array(virial_energies)
        t_arr = np.array(thermo_energies)

        v_mean = float(np.mean(v_arr))
        v_std = float(np.std(v_arr, ddof=1)) if len(v_arr) > 1 else 0.0
        t_mean = float(np.mean(t_arr))
        rg_mean = float(np.mean(rg_vals))

        return PIMCResult(
            beads_M=M,
            temperature_T=self.T,
            beta=round(self.beta, 4),
            total_sweeps=sweeps,
            burn_in=burn_in,
            acceptance_rate=round(accepted / max(1, total_proposals), 4),
            virial_energy_mean=round(v_mean, 6),
            virial_energy_std=round(v_std, 6),
            thermodynamic_energy_mean=round(t_mean, 6),
            exact_quantum_energy=round(self.exact_energy, 6),
            energy_bias=round(v_mean - self.exact_energy, 6),
            radius_of_gyration=round(rg_mean, 6),
            provenance_hash=prop_src.stream_digest,
            exact_finite_bead_energy=self.finite_bead_energy,
            discretization_bias=self.finite_bead_energy - self.exact_energy,
            virial_uncertainty=summarize_series(v_arr, reference=self.finite_bead_energy, batch_sizes=batch_sizes),
            thermodynamic_uncertainty=summarize_series(t_arr, reference=self.finite_bead_energy, batch_sizes=batch_sizes),
            sampler_protocol=("LOCAL_METROPOLIS_CENTROID_HEATBATH_V1"
                              if centroid_source is not None else "LOCAL_METROPOLIS_V1"),
            centroid_update=({
                "enabled": True,
                "schedule": "once after every local sweep, including burn-in",
                "total_updates": sweeps,
                "draw_bytes_per_update": 16,
                "centroid_draw_bytes": 16*sweeps,
                "source_manifest": centroid_source.manifest.to_dict(),
                "source_digest_before": centroid_digest_before,
                "source_digest_after": centroid_source.stream_digest,
                "source_byte_counter_delta": centroid_source.consumed_bytes-centroid_before,
                "source_shared_with_local_channels": any(
                    centroid_source is s for s in (prop_src, acc_src, init_src)),
                "independence_verified": False,
                "assumptions": ("IID uniform source words independent of retained internal coordinates; "
                                "no entropy or full-path mixing certification"),
            } if centroid_source is not None else {"enabled": False}),
            initialization_protocol=initialization_protocol,
            source_provenance={"roles": {role: ({**provenance_before[role],
                "consumed_bytes_after": source.consumed_bytes,
                "consumed_bytes_delta": source.consumed_bytes-provenance_before[role]["consumed_bytes_before"],
                "stream_digest_after": source.stream_digest,
                "aliases": [other for other, candidate in role_sources.items()
                            if candidate is source and other != role]}
                if source is not None else None) for role, source in role_sources.items()},
                "independence_verified": False,
                "manifest_scope": "Source identity fields; construction-time created_at omitted for deterministic replay",
                "accounting": ("Counters and digests identify consumed output, not entropy. "
                               "Aliased roles repeat the same whole-source delta; do not sum their deltas.")},
            trace=({"initial_path": trace_initial, "paths": trace_paths,
                    "TD": t_arr.copy(), "CV": v_arr.copy(),
                    "centroid": trace_centroid, "rg2": trace_rg2} if record_trace else None),
        )
