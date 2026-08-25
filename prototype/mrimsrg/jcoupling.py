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
from sympy.physics.wigner import clebsch_gordan, wigner_6j


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


@dataclass(frozen=True)
class MR1BContributions:
    """Gebrerufael Eq. (4.89b) IV--VI in the spherical one-body basis."""

    iv: np.ndarray
    v: np.ndarray
    vi: np.ndarray

    @property
    def total(self) -> np.ndarray:
        return self.iv + self.v + self.vi


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


def _block_map(blocks: tuple[CoupledBlock, ...]) -> dict[tuple[int, int, int], CoupledBlock]:
    return {(block.J, block.parity, block.Tz): block for block in blocks}


def _pair_phase(a: JOrbit, b: JOrbit, J: int) -> float:
    return -1.0 if ((a.j2 + b.j2) // 2 - J + 1) % 2 else 1.0


def _unnormalized_element(
    blocks: tuple[CoupledBlock, ...],
    block_map: dict[tuple[int, int, int], CoupledBlock],
    j_orbits: tuple[JOrbit, ...],
    J: int,
    a: int,
    b: int,
    c: int,
    d: int,
) -> float:
    oa, ob, oc, od = j_orbits[a], j_orbits[b], j_orbits[c], j_orbits[d]
    parity_bra = (oa.l + ob.l) % 2
    if parity_bra != (oc.l + od.l) % 2 or oa.tz2 + ob.tz2 != oc.tz2 + od.tz2:
        return 0.0
    block = block_map.get((J, parity_bra, (oa.tz2 + ob.tz2) // 2))
    if block is None:
        return 0.0
    phase = 1.0
    if a > b:
        phase *= _pair_phase(oa, ob, J)
        a, b = b, a
    if c > d:
        phase *= _pair_phase(oc, od, J)
        c, d = d, c
    try:
        ibra = block.pairs.index((a, b))
        iket = block.pairs.index((c, d))
    except ValueError:
        return 0.0
    pair_normalization = sqrt(2.0 if a == b else 1.0) * sqrt(2.0 if c == d else 1.0)
    return float(phase * pair_normalization * block.matrix[ibra, iket])


@lru_cache(maxsize=None)
def _six_j(j2a: int, j2b: int, J: int, j2c: int, j2d: int, Jp: int) -> float:
    try:
        return float(
            wigner_6j(
                Rational(j2a, 2),
                Rational(j2b, 2),
                J,
                Rational(j2c, 2),
                Rational(j2d, 2),
                Jp,
            )
        )
    except ValueError:
        return 0.0


def _pandya_element(
    blocks: tuple[CoupledBlock, ...],
    block_map: dict[tuple[int, int, int], CoupledBlock],
    j_orbits: tuple[JOrbit, ...],
    J: int,
    a: int,
    b: int,
    c: int,
    d: int,
) -> float:
    """Gebrerufael Eq. (4.49): O^J_(a bar b,c bar d)."""

    oa, ob, oc, od = j_orbits[a], j_orbits[b], j_orbits[c], j_orbits[d]
    Jp_min = max(abs(oa.j2 - od.j2), abs(oc.j2 - ob.j2)) // 2
    Jp_max = min(oa.j2 + od.j2, oc.j2 + ob.j2) // 2
    value = 0.0
    for Jp in range(Jp_min, Jp_max + 1):
        value -= (
            (2 * Jp + 1)
            * _six_j(oa.j2, ob.j2, J, oc.j2, od.j2, Jp)
            * _unnormalized_element(blocks, block_map, j_orbits, Jp, a, d, c, b)
        )
    return value


def mr_lambda2_one_body_coupled(
    x: tuple[CoupledBlock, ...],
    y: tuple[CoupledBlock, ...],
    lambda2: tuple[CoupledBlock, ...],
    orbits: np.ndarray,
    x_hermiticity: int,
    y_hermiticity: int,
) -> MR1BContributions:
    """Slow J-scheme oracle for the MR 2B--2B--lambda2 one-body term.

    ``x_hermiticity`` and ``y_hermiticity`` are +1 for Hermitian and -1 for
    anti-Hermitian inputs.  The implementation follows Gebrerufael Eq. (4.89b)
    IV--VI with explicit spherical-orbit sums.  IV and V carry an additional
    factor 1/2 relative to the thesis' unnormalized unrestricted-pair notation
    when evaluated through ``TwoBodyME`` normalized-pair blocks; the independent
    CG/m-scheme test below fixes this conversion.
    """

    if x_hermiticity not in (-1, 1) or y_hermiticity not in (-1, 1):
        raise ValueError("input hermiticity must be +1 or -1")
    j_orbits = extract_j_orbits(orbits)
    norb = len(j_orbits)
    maps = [_block_map(blocks) for blocks in (x, y, lambda2)]
    max_pair_J = max((oa.j2 + ob.j2) // 2 for oa in j_orbits for ob in j_orbits)

    def element(which: int, J: int, a: int, b: int, c: int, d: int) -> float:
        blocks = (x, y, lambda2)[which]
        return _unnormalized_element(blocks, maps[which], j_orbits, J, a, b, c, d)

    def pandya(which: int, J: int, a: int, b: int, c: int, d: int) -> float:
        blocks = (x, y, lambda2)[which]
        return _pandya_element(blocks, maps[which], j_orbits, J, a, b, c, d)

    def raw(one: int, two: int, x_index: int, y_index: int) -> np.ndarray:
        term = np.zeros(3)
        jhat1_sq = j_orbits[one].j2 + 1.0
        for J in range(max_pair_J + 1):
            Jhat_sq = 2 * J + 1.0
            for t in range(norb):
                for s in range(norb):
                    for w in range(norb):
                        for r in range(norb):
                            for v in range(norb):
                                # The 1/8 and 1/2 coefficients are the normalized-
                                # pair equivalents of thesis IV and V before the
                                # final one-body Hermiticity permutation.
                                term[0] += (
                                    Jhat_sq
                                    / (8.0 * jhat1_sq)
                                    * element(x_index, J, one, t, s, w)
                                    * element(y_index, J, r, v, two, t)
                                    * element(2, J, r, v, s, w)
                                )
                                term[1] += (
                                    Jhat_sq
                                    / (2.0 * jhat1_sq)
                                    * pandya(x_index, J, one, t, s, r)
                                    * pandya(2, J, s, r, v, w)
                                    * pandya(y_index, J, v, w, two, t)
                                )
        for J1 in range(max_pair_J + 1):
            for J2 in range(max_pair_J + 1):
                angular_weight = (2 * J1 + 1.0) * (2 * J2 + 1.0)
                for t in range(norb):
                    for s in range(norb):
                        if j_orbits[t].j2 != j_orbits[s].j2:
                            continue
                        for r in range(norb):
                            for v in range(norb):
                                for w in range(norb):
                                    term[2] -= (
                                        angular_weight
                                        / (2.0 * jhat1_sq * (j_orbits[t].j2 + 1.0))
                                        * element(x_index, J1, one, t, two, s)
                                        * element(2, J2, s, w, r, v)
                                        * element(y_index, J2, r, v, t, w)
                                    )
        return term

    raw_terms = np.zeros((norb, norb, 3))
    for one in range(norb):
        for two in range(norb):
            if (
                j_orbits[one].l != j_orbits[two].l
                or j_orbits[one].j2 != j_orbits[two].j2
                or j_orbits[one].tz2 != j_orbits[two].tz2
            ):
                continue
            raw_terms[one, two] = raw(one, two, 0, 1) - raw(one, two, 1, 0)

    output_hermiticity = -x_hermiticity * y_hermiticity
    completed = raw_terms + output_hermiticity * raw_terms.transpose(1, 0, 2)
    return MR1BContributions(completed[:, :, 0], completed[:, :, 1], completed[:, :, 2])
