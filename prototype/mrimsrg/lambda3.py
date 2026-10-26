"""Rank-0 three-body cumulant lambda3 of a J = 0 reference, J coupled.

The m-scheme cumulant uses the gamma2 convention of densities.py,

    gamma3[abc, def] = <Psi| a^dagger_a a^dagger_b a^dagger_c a_f a_e a_d |Psi>
    lambda3          = gamma3 - A(gamma1 gamma1 gamma1) - A(gamma1 lambda2)

and is evaluated in the natural-orbit basis, where gamma1 is diagonal, so
A(ggg) only lives on the triple diagonal and A(g lambda2) only connects
triples that share an orbital.  It vanishes unless every index is a
fractionally occupied orbit, so only those "active" orbits are kept.

The coupled block follows imsrg++ ThreeBodyME::GetME_pn: plain Clebsch-Gordan
coupling of the antisymmetric m-scheme tensor, first (ab) to J12, then c,

    Lambda^J[(ab)J12 c, (de)J45 f]
        = sum_m <ja ma jb mb|J12 M12> <J12 M12 jc mc|J M> (same for def)
          * lambda3[(a ma)(b mb)(c mc), (d md)(e me)(f mf)],

averaged over M (rank-0 projection, exact for a J = 0 state).  States are
stored with a <= b; (ba) follows from (ab) by the phase (-1)^(ja+jb-J12+1).
The full m-scheme contraction is then
    sum_m v lambda3 = sum_J (2J+1) sum_{unrestricted states} V^J Lambda^J.
"""

from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from itertools import combinations, product

import numpy as np
from scipy import sparse
from sympy import Rational
from sympy.physics.wigner import clebsch_gordan


@lru_cache(maxsize=None)
def cg2(j1: int, m1: int, j2: int, m2: int, J: int, M: int) -> float:
    """<j1 m1 j2 m2 | J M> with all arguments doubled."""
    if m1 + m2 != M or abs(m1) > j1 or abs(m2) > j2 or abs(M) > J:
        return 0.0
    if J < abs(j1 - j2) or J > j1 + j2 or (j1 + j2 + J) % 2:
        return 0.0
    return float(clebsch_gordan(Rational(j1, 2), Rational(j2, 2), Rational(J, 2),
                                Rational(m1, 2), Rational(m2, 2), Rational(M, 2)))


def annihilation_matrix(occupation, coefficients, columns):
    """Rows: reduced determinants; column (p, q, ...) holds ... a_q a_p |Psi>.

    With this ordering, A_k^T A_k is the k-body density in the
    <a^dagger_p a^dagger_q ... a_s a_r> convention of densities.py.
    """
    below = np.cumsum(occupation, axis=1, dtype=np.int32) - occupation
    packed = np.packbits(occupation, axis=1)
    keys, cols, values = [], [], []
    for column, orbit_tuple in enumerate(columns):
        mask = np.all(occupation[:, list(orbit_tuple)] == 1, axis=1)
        rows = np.flatnonzero(mask)
        if rows.size == 0:
            continue
        reduced = packed[rows].copy()
        exponent = np.zeros(rows.size, dtype=np.int64)
        for k, p in enumerate(orbit_tuple):
            reduced[:, p // 8] &= np.uint8(0xFF ^ (0x80 >> (p % 8)))
            exponent += below[rows, p] - k  # earlier-removed orbits lie below p
        keys.append(reduced)
        cols.append(np.full(rows.size, column))
        values.append(coefficients[rows] * np.where(exponent % 2, -1.0, 1.0))
    if not keys:
        return sparse.csr_matrix((1, len(columns)))
    keys = np.concatenate(keys)
    view = np.ascontiguousarray(keys).view(np.dtype((np.void, keys.shape[1]))).ravel()
    _, row_index = np.unique(view, return_inverse=True)
    return sparse.csr_matrix((np.concatenate(values), (row_index.ravel(), np.concatenate(cols))),
                             shape=(row_index.max() + 1, len(columns)))


def _permutation_sign(order) -> float:
    order = list(order)
    sign = 1.0
    for i in range(len(order)):
        for j in range(i + 1, len(order)):
            if order[i] > order[j]:
                sign = -sign
    return sign


def mscheme_rotation(mstates, transformation: np.ndarray) -> np.ndarray:
    """U[old m-state, new m-state] = C[old orbit, new orbit] for equal m."""
    size = len(mstates)
    U = np.zeros((size, size))
    for p, (a, ma) in enumerate(mstates):
        for q, (b, mb) in enumerate(mstates):
            if ma == mb:
                U[p, q] = transformation[a, b]
    return U


def natural_amplitudes(occupation, coefficients, U, natural_columns):
    """Annihilation amplitudes for sorted natural-orbit tuples.

    a_{t'} = sum_t det(U[t, t']) a_t, so each natural column is a short
    combination of oscillator columns (U only mixes radial orbits).
    """
    support = [np.flatnonzero(np.abs(U[:, x]) > 0.0) for x in range(U.shape[1])]
    ho_index: dict[tuple, int] = {}
    rows, cols, values = [], [], []
    for column, target in enumerate(natural_columns):
        for combo in product(*(support[x] for x in target)):
            if len(set(combo)) < len(combo):
                continue
            value = 1.0
            for p, x in zip(combo, target):
                value *= U[p, x]
            order = np.argsort(combo)
            key = tuple(int(combo[i]) for i in order)
            rows.append(ho_index.setdefault(key, len(ho_index)))
            cols.append(column)
            values.append(_permutation_sign(order) * value)
    R = sparse.csr_matrix((values, (rows, cols)), shape=(len(ho_index), len(natural_columns)))
    ho_columns = sorted(ho_index, key=ho_index.get)
    return annihilation_matrix(occupation, coefficients, ho_columns) @ R


def coupled_lambda3(orbits, mstates, occupation, coefficients, transformation, occupations,
                    *, active_tolerance: float = 1e-12):
    """J-coupled rank-0 lambda3 in the natural-orbit basis.

    orbits: (n, l, 2j, 2tz) per orbit; mstates: (orbit, 2m) per m-state in
    the column order of `occupation` (determinants x m-states, oscillator
    basis); transformation: C[old orbit, new orbit]; occupations: natural
    occupation per orbit.  Returns (blocks, report): blocks are
    (2J, parity, 2Tz, states[(a, b, c, J12)], matrix) with a <= b.
    """
    occupations = np.asarray(occupations, dtype=float)
    active = [a for a in range(len(orbits))
              if active_tolerance < occupations[a] < 1.0 - active_tolerance]
    active_set = set(active)
    active_m = [k for k, (a, _) in enumerate(mstates) if a in active_set]
    by_orbit = defaultdict(dict)
    for k, (a, m2) in enumerate(mstates):
        by_orbit[a][m2] = k
    n_m = np.array([occupations[a] for a, _ in mstates])
    U = mscheme_rotation(mstates, transformation)

    # lambda2 on active natural m-scheme pairs.
    pairs = list(combinations(active_m, 2))
    pair_index = {pq: i for i, pq in enumerate(pairs)}
    A2 = natural_amplitudes(occupation, coefficients, U, pairs)
    lambda2 = (A2.T @ A2).toarray()
    for i, (x, y) in enumerate(pairs):
        lambda2[i, i] -= n_m[x] * n_m[y]

    triples = list(combinations(active_m, 3))
    triple_index = {t: i for i, t in enumerate(triples)}
    A3 = natural_amplitudes(occupation, coefficients, U, triples)

    # A(gamma1 lambda2)[t, t'] = sum over shared z of (-1)^(i+j) n_z lambda2(t\z, t'\z)
    # with i, j the positions of z in t and t'.
    nonzero = np.argwhere(np.abs(lambda2) > 0.0)
    rows, cols, values = [], [], []
    for z in active_m:
        for i, j in nonzero:
            P, Q = pairs[i], pairs[j]
            if z in P or z in Q:
                continue
            t = tuple(sorted(P + (z,)))
            u = tuple(sorted(Q + (z,)))
            sign = -1.0 if (t.index(z) + u.index(z)) % 2 else 1.0
            rows.append(triple_index[t])
            cols.append(triple_index[u])
            values.append(sign * n_m[z] * lambda2[i, j])
    g_lambda2 = sparse.csr_matrix((values, (rows, cols)), shape=(len(triples), len(triples)))
    disconnected = sparse.diags([n_m[a] * n_m[b] * n_m[c] for a, b, c in triples])

    channels = defaultdict(list)
    for a in active:
        _, la, ja, ta = orbits[a]
        for b in active:
            if b < a:
                continue
            _, lb, jb, tb = orbits[b]
            for J12 in range(abs(ja - jb) // 2, (ja + jb) // 2 + 1):
                if a == b and J12 % 2:
                    continue
                for c in active:
                    _, lc, jc, tc = orbits[c]
                    for twoJ in range(abs(2 * J12 - jc), 2 * J12 + jc + 1, 2):
                        channels[(twoJ, (la + lb + lc) % 2, ta + tb + tc)].append((a, b, c, J12))

    blocks = []
    max_asymmetry = 0.0
    for (twoJ, parity, twoTz), states in sorted(channels.items()):
        matrix = np.zeros((len(states), len(states)))
        for twoM in range(-twoJ, twoJ + 1, 2):
            P_rows, P_cols, P_values = [], [], []
            for column, (a, b, c, J12) in enumerate(states):
                ja, jb, jc = orbits[a][2], orbits[b][2], orbits[c][2]
                for ma in range(-ja, ja + 1, 2):
                    for mb in range(-jb, jb + 1, 2):
                        mc = twoM - ma - mb
                        if abs(mc) > jc:
                            continue
                        weight = (cg2(ja, ma, jb, mb, 2 * J12, ma + mb)
                                  * cg2(2 * J12, ma + mb, jc, mc, twoJ, twoM))
                        if weight == 0.0:
                            continue
                        combo = (by_orbit[a][ma], by_orbit[b][mb], by_orbit[c][mc])
                        if len(set(combo)) < 3:
                            continue
                        order = np.argsort(combo)
                        key = tuple(combo[i] for i in order)
                        P_rows.append(triple_index[key])
                        P_cols.append(column)
                        P_values.append(_permutation_sign(order) * weight)
            P = sparse.csr_matrix((P_values, (P_rows, P_cols)), shape=(len(triples), len(states)))
            B = A3 @ P
            matrix += (B.T @ B).toarray()
            matrix -= (P.T @ (disconnected @ P)).toarray()
            matrix -= (P.T @ (g_lambda2 @ P)).toarray()
        matrix /= twoJ + 1
        max_asymmetry = max(max_asymmetry, float(np.max(np.abs(matrix - matrix.T))) if matrix.size else 0.0)
        blocks.append((twoJ, parity, twoTz, states, 0.5 * (matrix + matrix.T)))
    report = {
        "active orbits": len(active),
        "lambda3 channels": len(blocks),
        "lambda3 states": sum(len(s) for _, _, _, s, _ in blocks),
        "lambda3 max asymmetry": max_asymmetry,
        "lambda3 max |element|": max((float(np.max(np.abs(m))) for *_, m in blocks if m.size), default=0.0),
    }
    return blocks, report
