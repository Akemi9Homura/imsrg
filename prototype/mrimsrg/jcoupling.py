"""Independent scalar two-body J/m coupling oracle.

The production ``imsrg++`` matrices store normalized antisymmetrized pair
matrix elements.  This module implements the same convention directly from
Clebsch--Gordan coefficients so that MR density and commutator code can be
checked without calling the production implementation.

For canonical J-orbits ``a <= b`` and ``c <= d`` we store

    O_J(ab,cd) = 1 / [N_ab N_cd (2J+1)]
                 sum_{M,m_a,...} C(ab;JM) C(cd;JM) O_m(ab,cd),

where ``N_ab=sqrt(1+delta_ab)``.  This is the convention used by the
``TwoBodyME`` matrices and by the validated ``jcoupled64`` bridge.  A rank-zero
FCIQMC/shell-model-obs TBTD is a reduced density and must therefore be divided
by ``sqrt(2J+1)`` before being stored in these blocks.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import sqrt

import numpy as np
from sympy import Rational
from sympy.physics.wigner import clebsch_gordan


@dataclass(frozen=True)
class JOrbit:
    index: int
    n: int
    l: int
    j2: int
    tz2: int
    substates: tuple[int, ...]


@dataclass(frozen=True)
class CoupledBlock:
    J: int
    parity: int
    Tz: int
    pairs: tuple[tuple[int, int], ...]
    matrix: np.ndarray


@lru_cache(maxsize=None)
def _cg(j2a: int, j2b: int, J: int, m2a: int, m2b: int) -> float:
    M2 = m2a + m2b
    if abs(M2) > 2 * J:
        return 0.0
    return float(
        clebsch_gordan(
            Rational(j2a, 2),
            Rational(j2b, 2),
            J,
            Rational(m2a, 2),
            Rational(m2b, 2),
            Rational(M2, 2),
        )
    )


def extract_j_orbits(orbits: np.ndarray) -> tuple[JOrbit, ...]:
    """Validate the prototype orbit table and return its spherical multiplets."""

    if orbits.ndim != 2 or orbits.shape[1] != 6:
        raise ValueError("orbits must have columns (jindex,n,l,j2,m2,tz2)")
    result: list[JOrbit] = []
    for index in sorted(int(value) for value in np.unique(orbits[:, 0])):
        rows = np.flatnonzero(orbits[:, 0] == index)
        first = orbits[rows[0]]
        n, l, j2, tz2 = map(int, (first[1], first[2], first[3], first[5]))
        if np.any(orbits[rows, 1] != n) or np.any(orbits[rows, 2] != l):
            raise ValueError(f"J-orbit {index} has inconsistent n/l labels")
        if np.any(orbits[rows, 3] != j2) or np.any(orbits[rows, 5] != tz2):
            raise ValueError(f"J-orbit {index} has inconsistent j/tz labels")
        expected_m2 = np.arange(-j2, j2 + 1, 2, dtype=orbits.dtype)
        actual_m2 = np.sort(orbits[rows, 4])
        if not np.array_equal(actual_m2, expected_m2):
            raise ValueError(f"J-orbit {index} does not contain a complete j multiplet")
        ordered_rows = tuple(int(rows[np.where(orbits[rows, 4] == m2)[0][0]]) for m2 in expected_m2)
        result.append(JOrbit(index, n, l, j2, tz2, ordered_rows))
    if tuple(orbit.index for orbit in result) != tuple(range(len(result))):
        raise ValueError("J-orbit indices must be contiguous and zero based")
    return tuple(result)


def _allowed_pair_J(a: JOrbit, b: JOrbit, J: int) -> bool:
    if 2 * J < abs(a.j2 - b.j2) or 2 * J > a.j2 + b.j2:
        return False
    return a.index != b.index or J % 2 == 0


def _channel_pairs(j_orbits: tuple[JOrbit, ...], J: int, parity: int, Tz: int) -> tuple[tuple[int, int], ...]:
    pairs: list[tuple[int, int]] = []
    for a in j_orbits:
        for b in j_orbits[a.index :]:
            if (a.l + b.l) % 2 != parity or (a.tz2 + b.tz2) // 2 != Tz:
                continue
            if _allowed_pair_J(a, b, J):
                pairs.append((a.index, b.index))
    return tuple(pairs)


def _all_channels(j_orbits: tuple[JOrbit, ...]) -> tuple[tuple[int, int, int], ...]:
    labels: set[tuple[int, int, int]] = set()
    for a in j_orbits:
        for b in j_orbits[a.index :]:
            for J in range(abs(a.j2 - b.j2) // 2, (a.j2 + b.j2) // 2 + 1):
                if _allowed_pair_J(a, b, J):
                    labels.add((J, (a.l + b.l) % 2, (a.tz2 + b.tz2) // 2))
    return tuple(sorted(labels))


def couple_scalar_two_body(tensor: np.ndarray, orbits: np.ndarray) -> tuple[CoupledBlock, ...]:
    """Project an antisymmetrized scalar m-scheme tensor to normalized J blocks."""

    j_orbits = extract_j_orbits(orbits)
    norb = orbits.shape[0]
    if tensor.shape != (norb, norb, norb, norb):
        raise ValueError("two-body tensor shape does not match the orbit table")
    blocks: list[CoupledBlock] = []
    for J, parity, Tz in _all_channels(j_orbits):
        pairs = _channel_pairs(j_orbits, J, parity, Tz)
        matrix = np.zeros((len(pairs), len(pairs)), dtype=np.result_type(tensor, np.float64))
        for ibra, (a, b) in enumerate(pairs):
            oa, ob = j_orbits[a], j_orbits[b]
            norm_ab = sqrt(2.0 if a == b else 1.0)
            for iket, (c, d) in enumerate(pairs):
                oc, od = j_orbits[c], j_orbits[d]
                norm_cd = sqrt(2.0 if c == d else 1.0)
                value = 0.0
                for p in oa.substates:
                    m2p = int(orbits[p, 4])
                    for q in ob.substates:
                        m2q = int(orbits[q, 4])
                        cg_ab = _cg(oa.j2, ob.j2, J, m2p, m2q)
                        if cg_ab == 0.0:
                            continue
                        M2 = m2p + m2q
                        for r in oc.substates:
                            m2r = int(orbits[r, 4])
                            for s in od.substates:
                                m2s = int(orbits[s, 4])
                                if m2r + m2s != M2:
                                    continue
                                cg_cd = _cg(oc.j2, od.j2, J, m2r, m2s)
                                if cg_cd != 0.0:
                                    value += cg_ab * cg_cd * tensor[p, q, r, s]
                matrix[ibra, iket] = value / (norm_ab * norm_cd * (2 * J + 1))
        blocks.append(CoupledBlock(J, parity, Tz, pairs, matrix))
    return tuple(blocks)


def reconstruct_scalar_two_body(
    blocks: tuple[CoupledBlock, ...], orbits: np.ndarray
) -> np.ndarray:
    """Expand normalized J blocks to the antisymmetrized m-scheme tensor."""

    j_orbits = extract_j_orbits(orbits)
    by_label = {(block.J, block.parity, block.Tz): block for block in blocks}
    pair_positions = {
        label: {pair: position for position, pair in enumerate(block.pairs)}
        for label, block in by_label.items()
    }
    norb = orbits.shape[0]
    tensor = np.zeros((norb, norb, norb, norb), dtype=np.float64)
    for p in range(norb):
        op = j_orbits[int(orbits[p, 0])]
        for q in range(norb):
            if p == q:
                continue
            oq = j_orbits[int(orbits[q, 0])]
            for r in range(norb):
                o_r = j_orbits[int(orbits[r, 0])]
                for s in range(norb):
                    if r == s or int(orbits[p, 4] + orbits[q, 4]) != int(orbits[r, 4] + orbits[s, 4]):
                        continue
                    os = j_orbits[int(orbits[s, 0])]
                    if (op.l + oq.l - o_r.l - os.l) % 2 != 0 or op.tz2 + oq.tz2 != o_r.tz2 + os.tz2:
                        continue
                    a, b = op.index, oq.index
                    c, d = o_r.index, os.index
                    swap_ab = a > b
                    swap_cd = c > d
                    if swap_ab:
                        a, b = b, a
                    if swap_cd:
                        c, d = d, c
                    Jmin = max(abs(op.j2 - oq.j2), abs(o_r.j2 - os.j2)) // 2
                    Jmax = min(op.j2 + oq.j2, o_r.j2 + os.j2) // 2
                    value = 0.0
                    for J in range(Jmin, Jmax + 1):
                        label = (J, (op.l + oq.l) % 2, (op.tz2 + oq.tz2) // 2)
                        block = by_label.get(label)
                        if block is None:
                            continue
                        ibra = pair_positions[label].get((a, b))
                        iket = pair_positions[label].get((c, d))
                        if ibra is None or iket is None:
                            continue
                        phase = sqrt(2.0 if a == b else 1.0) * sqrt(2.0 if c == d else 1.0)
                        if swap_ab:
                            phase *= -1.0 if (op.j2 + oq.j2) // 2 - J + 1 & 1 else 1.0
                        if swap_cd:
                            phase *= -1.0 if (o_r.j2 + os.j2) // 2 - J + 1 & 1 else 1.0
                        value += (
                            phase
                            * block.matrix[ibra, iket]
                            * _cg(op.j2, oq.j2, J, int(orbits[p, 4]), int(orbits[q, 4]))
                            * _cg(o_r.j2, os.j2, J, int(orbits[r, 4]), int(orbits[s, 4]))
                        )
                    tensor[p, q, r, s] = value
    return tensor


def coupled_scalar_contraction(
    left: tuple[CoupledBlock, ...], right: tuple[CoupledBlock, ...]
) -> float:
    """Return ``1/4 sum_m left_m right_m`` from normalized coupled blocks."""

    right_by_label = {(block.J, block.parity, block.Tz): block for block in right}
    value = 0.0
    for block in left:
        label = (block.J, block.parity, block.Tz)
        other = right_by_label.get(label)
        if other is None or other.pairs != block.pairs:
            raise ValueError(f"incompatible coupled block {label}")
        value += (2 * block.J + 1) * np.einsum(
            "ij,ij->", block.matrix, other.matrix, optimize=True
        )
    return float(value)
