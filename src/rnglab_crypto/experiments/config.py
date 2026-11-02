from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class CampaignConfig:
    namespace: str = "blue-confirm-v1"
    repetitions: int = 128
    calibration_repetitions: int = 199
    operations: int = 256
    fault_onset: int = 64
    diagnostic_bytes: int = 4096
    weak_seed_bits: tuple[int, ...] = (8, 12)
    nonce_bits: tuple[int, ...] = (8, 12, 16)
    synthetic_only: bool = True

    def __post_init__(self):
        for name, lo, hi in (("repetitions", 2, 512), ("calibration_repetitions", 19, 999),
                             ("operations", 32, 4096), ("fault_onset", 1, 4095),
                             ("diagnostic_bytes", 128, 65536)):
            value = getattr(self, name)
            if type(value) is not int or not lo <= value <= hi:
                raise ValueError(f"invalid {name}")
        if not isinstance(self.namespace, str) or not self.namespace or len(self.namespace) > 100:
            raise ValueError("bounded nonempty seed namespace required")
        if self.synthetic_only is not True or self.fault_onset >= self.operations:
            raise ValueError("synthetic_only and fault within observation horizon required")
        for name, maximum in (("weak_seed_bits", 16), ("nonce_bits", 24)):
            values = tuple(getattr(self, name))
            if not values or len(values) > 4 or len(set(values)) != len(values) or any(type(v) is not int or not 1 <= v <= maximum for v in values):
                raise ValueError(f"invalid bounded {name}")
            object.__setattr__(self, name, values)
        if self.repetitions * sum(2**b for b in self.weak_seed_bits) > 2_000_000:
            raise ValueError("bounded enumeration budget exceeded")
        if self.repetitions*self.operations*(len(self.nonce_bits)+6) > 4_000_000:
            raise ValueError("trajectory budget exceeded")

    def to_dict(self):
        return asdict(self)

    @classmethod
    def smoke(cls):
        return cls(namespace="blue-smoke-v1", repetitions=4, calibration_repetitions=19,
                   operations=64, fault_onset=16, diagnostic_bytes=512,
                   weak_seed_bits=(4,), nonce_bits=(4, 8))
