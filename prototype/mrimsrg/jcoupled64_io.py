"""Dependency-light reader and writer for the lossless J-coupled checkpoint."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import struct
from typing import Any

import numpy as np


J64_MAGIC = b"mrimsrg_j64_v1\0\0"


@dataclass(frozen=True)
class JCoupled64:
    hw: float
    emax: int
    orbits: np.ndarray
    zero_body: float
    one_body: np.ndarray
    records: tuple[tuple[int, int, int, int, int, float], ...]


def _read_exact(stream: Any, size: int, field: str) -> bytes:
    value = stream.read(size)
    if len(value) != size:
        raise ValueError(f"truncated jcoupled64 while reading {field}")
    return value


def read_jcoupled64(path: Path) -> JCoupled64:
    with path.open("rb") as stream:
        if _read_exact(stream, 16, "magic") != J64_MAGIC:
            raise ValueError("unsupported jcoupled64 payload")
        hw = struct.unpack("<d", _read_exact(stream, 8, "hw"))[0]
        emax = struct.unpack("<i", _read_exact(stream, 4, "emax"))[0]
        norb, nobme, ntbme = struct.unpack(
            "<QQQ", _read_exact(stream, 24, "counts")
        )
        orbits = np.asarray(
            [
                struct.unpack("<iiii", _read_exact(stream, 16, "orbit"))
                for _ in range(norb)
            ],
            dtype=np.int64,
        )
        zero_body = struct.unpack("<d", _read_exact(stream, 8, "zero body"))[0]
        if nobme != norb * (norb + 1) // 2:
            raise ValueError("invalid jcoupled64 one-body count")
        one_body = np.zeros((norb, norb), dtype=np.float64)
        for a in range(norb):
            for b in range(a, norb):
                value = struct.unpack("<d", _read_exact(stream, 8, "OBME"))[0]
                one_body[a, b] = value
                one_body[b, a] = value
        records = []
        for _ in range(ntbme):
            a, b, c, d, j = struct.unpack(
                "<iiiii", _read_exact(stream, 20, "TBME indices")
            )
            value = struct.unpack("<d", _read_exact(stream, 8, "TBME"))[0]
            records.append((a, b, c, d, j, value))
        if stream.read(1):
            raise ValueError("jcoupled64 payload has trailing bytes")
    return JCoupled64(hw, emax, orbits, zero_body, one_body, tuple(records))


def write_jcoupled64(path: Path, payload: JCoupled64, operator: Any) -> None:
    """Materialize an imsrg++ scalar operator in the acceptance-only format."""
    with path.open("xb") as stream:
        stream.write(J64_MAGIC)
        stream.write(
            struct.pack(
                "<diQQQ",
                payload.hw,
                payload.emax,
                len(payload.orbits),
                len(payload.orbits) * (len(payload.orbits) + 1) // 2,
                len(payload.records),
            )
        )
        for orbit in payload.orbits:
            stream.write(struct.pack("<iiii", *(int(value) for value in orbit)))
        stream.write(struct.pack("<d", float(operator.ZeroBody)))
        for a in range(len(payload.orbits)):
            for b in range(a, len(payload.orbits)):
                stream.write(struct.pack("<d", float(operator.GetOneBody(a, b))))
        for a, b, c, d, j, _ in payload.records:
            value = float(operator.TwoBody.GetTBME_J_norm(j, j, a, b, c, d))
            stream.write(struct.pack("<iiiiid", a, b, c, d, j, value))
