"""RNGLab v0.14: Physical Randomness Acquisition Adapters & Emulators.

Interfaces for physical entropy sources:
1. Geiger-Müller Radioactive Decay Event Stream:
   - Models Poisson inter-arrival intervals: Delta t ~ Exp(lambda)
   - Detector dead time: Paralyzable vs Non-paralyzable modes
   - Timer quantization & clock jitter
2. Photodiode & Optical Noise:
   - ADC sampling (8, 12, 16-bit)
   - Thermal drift and electronic shot noise
3. Microcontroller Serial Line Adapter:
   - Packet framing, CRC-16 checks, dropped packet detection
"""
from __future__ import annotations

import math
import random
import struct
import time
from typing import List, Optional

from .base import RandomSource, SourceManifest


class GeigerMullerPulseSource(RandomSource):
    """Geiger-Müller radiation detector event timing stream.
    
    Can operate in:
    - Real hardware mode (reading timestamped pulse events)
    - Synthetic calibration mode (simulating Poisson decay with detector dead time)
    """

    def __init__(
        self,
        mean_cps: float = 120.0,
        dead_time_seconds: float = 90e-6,  # 90 microseconds typical GM dead time
        paralyzable: bool = False,
        timer_clock_hz: float = 1e6,       # 1 MHz timer resolution (1 microsecond)
        seed: Optional[int] = 42,
        is_synthetic: bool = True,
        max_bytes: Optional[int] = None,
    ):
        if not is_synthetic:
            raise ValueError("Live hardware acquisition is not implemented; use the explicit emulator")
        if not all(math.isfinite(v) for v in [mean_cps,timer_clock_hz,dead_time_seconds]) or mean_cps <= 0 or timer_clock_hz <= 0 or dead_time_seconds < 0:
            raise ValueError("invalid emulator parameters")
        manifest = SourceManifest(
            source_id="gm_detector_sim" if is_synthetic else "gm_detector_live",
            family="physical_gm",
            is_deterministic=is_synthetic,
            is_physical=not is_synthetic,
            seed=seed,
            metadata={
                "mean_cps": mean_cps,
                "dead_time_seconds": dead_time_seconds,
                "paralyzable": paralyzable,
                "timer_clock_hz": timer_clock_hz,
                "is_synthetic": is_synthetic,
                "physics": "Poisson decay with detector dead time",
            },
        )
        super().__init__(manifest, max_bytes=max_bytes)
        self.mean_cps = mean_cps
        self.dead_time = dead_time_seconds
        self.paralyzable = paralyzable
        self.clock_hz = timer_clock_hz
        self._rng = random.Random(seed)
        self._current_time = 0.0
        self._buffer = bytearray()

    def _next_event_interval(self) -> float:
        """Observed intervals for Poisson incident pulses, not hardware acquisition.

        Nonparalyzable: d+Exp(lambda). Paralyzable: sum gaps until a gap>=d;
        every lost pulse extends the exclusion interval. A finite loop budget
        rejects pathological dead-time/rate combinations rather than hanging.
        """
        if not self.paralyzable:
            return self.dead_time+self._rng.expovariate(self.mean_cps)
        elapsed=0.
        for _ in range(100000):
            gap=self._rng.expovariate(self.mean_cps);elapsed+=gap
            if gap>=self.dead_time:return elapsed
        raise ValueError("paralyzable emulator event budget exhausted")

    def _read_bytes(self, n: int) -> bytes:
        while len(self._buffer) < n:
            dt = self._next_event_interval()
            # Quantize interval by timer clock ticks
            ticks = int(round(dt * self.clock_hz))
            # Extract least significant byte of quantized arrival time
            self._buffer.append(ticks & 0xFF)
        out = bytes(self._buffer[:n])
        del self._buffer[:n]
        return out


class PhotodiodeNoiseSource(RandomSource):
    """Photodiode / Optical noise digitized via an Analog-to-Digital Converter (ADC).
    
    Models:
    - Thermal Gaussian shot noise
    - Configurable ADC bit depth (e.g. 8, 12, or 16-bit)
    - Slow thermal baseline drift
    """

    def __init__(
        self,
        adc_bits: int = 12,
        noise_std_dev: float = 15.0,
        thermal_drift_rate: float = 0.001,
        seed: Optional[int] = 12345,
        max_bytes: Optional[int] = None,
    ):
        manifest = SourceManifest(
            source_id=f"photodiode_adc{adc_bits}b",
            family="physical_optical",
            is_deterministic=True,
            is_physical=False,
            seed=seed,
            metadata={
                "adc_bits": adc_bits,
                "noise_std_dev": noise_std_dev,
                "thermal_drift_rate": thermal_drift_rate,
            },
        )
        super().__init__(manifest, max_bytes=max_bytes)
        self.adc_bits = adc_bits
        self.max_adc = (1 << adc_bits) - 1
        self.midpoint = self.max_adc / 2.0
        self.noise_std = noise_std_dev
        self.drift_rate = thermal_drift_rate
        self.step_idx = 0
        self._rng = random.Random(seed)
        self._buffer = bytearray()

    def _sample_adc(self) -> int:
        # Baseline drift
        drift = math.sin(self.step_idx * self.drift_rate) * (self.max_adc * 0.1)
        mean_level = self.midpoint + drift
        # Gaussian noise
        raw_val = self._rng.gauss(mean_level, self.noise_std)
        # ADC quantization and clipping
        quantized = max(0, min(self.max_adc, int(round(raw_val))))
        self.step_idx += 1
        return quantized

    def _read_bytes(self, n: int) -> bytes:
        while len(self._buffer) < n:
            val = self._sample_adc()
            # Extract lower 8 bits of ADC sample
            self._buffer.append(val & 0xFF)
        out = bytes(self._buffer[:n])
        del self._buffer[:n]
        return out


class MicrocontrollerSerialSource(RandomSource):
    """Serial stream parser for microcontrollers (Raspberry Pi, ESP32, Arduino).
    
    Frames incoming byte packets with header [0xAA, 0x55], length, payload, and CRC-16.
    Detects and accounts for dropped or corrupt packets.
    """

    def __init__(
        self,
        packet_payload_size: int = 32,
        simulated_loss_rate: float = 0.0,
        seed: Optional[int] = 999,
        max_bytes: Optional[int] = None,
    ):
        manifest = SourceManifest(
            source_id="mcu_serial_adapter",
            family="physical_serial",
            is_deterministic=True,
            is_physical=False,
            seed=seed,
            metadata={
                "framing": "SYNC [0xAA, 0x55] + LEN + PAYLOAD + CRC16",
                "simulated_loss_rate": simulated_loss_rate,
            },
        )
        super().__init__(manifest, max_bytes=max_bytes)
        self.payload_size = packet_payload_size
        self.loss_rate = simulated_loss_rate
        self._rng = random.Random(seed)
        self._buffer = bytearray()
        self.dropped_packets = 0
        self.valid_packets = 0

    @staticmethod
    def _crc16(data: bytes) -> int:
        crc = 0xFFFF
        for byte in data:
            crc ^= byte
            for _ in range(8):
                if crc & 1:
                    crc = (crc >> 1) ^ 0xA001
                else:
                    crc >>= 1
        return crc & 0xFFFF

    def _generate_framed_packet(self) -> bytes:
        payload = self._rng.randbytes(self.payload_size)
        crc = self._crc16(payload)
        frame = b"\xAA\x55" + struct.pack(">H", len(payload)) + payload + struct.pack(">H", crc)
        return frame

    def _read_bytes(self, n: int) -> bytes:
        while len(self._buffer) < n:
            frame = self._generate_framed_packet()
            # Check for simulated packet loss
            if self.loss_rate > 0.0 and self._rng.random() < self.loss_rate:
                self.dropped_packets += 1
                continue

            # Verify sync and CRC
            sync = frame[:2]
            length = struct.unpack(">H", frame[2:4])[0]
            payload = frame[4 : 4 + length]
            received_crc = struct.unpack(">H", frame[4 + length : 6 + length])[0]

            if sync == b"\xAA\x55" and self._crc16(payload) == received_crc:
                self.valid_packets += 1
                self._buffer.extend(payload)
            else:
                self.dropped_packets += 1

        out = bytes(self._buffer[:n])
        del self._buffer[:n]
        return out
