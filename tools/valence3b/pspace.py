#!/usr/bin/env python3
"""Dense m-scheme valence-space solver with an optional three-body interaction (test oracle only).

Reads a KSHELL .snt (zero-body term from the "Zero body term" comment, one-body, normalized J-coupled
two-body) and, optionally, a valence three-body file written by ReadWrite::WriteValence3body
(lines "a b c d e f Jab Jde 2J V"; V is the un-normalized J-coupled element of the antisymmetric
m-scheme matrix element, (ab) coupled to Jab first, then c to J).  J-coupled elements are turned into
m-scheme elements by writing the coupled states on canonical Slater determinants and inverting that
expansion with a pseudo-inverse, which also covers repeated orbits.  Many-body matrices are built on
bit-string determinants of the valence m-states and diagonalized densely: meant for spaces of a few
hundred determinants (p shell), not for production.
"""

from itertools import combinations
from math import sqrt

import numpy as np
from sympy.physics.quantum.cg import CG
from sympy import S

_CG_CACHE = {}


def cg(j1x2, m1x2, j2x2, m2x2, Jx2, Mx2):
    key = (j1x2, m1x2, j2x2, m2x2, Jx2, Mx2)
    if key not in _CG_CACHE:
        if m1x2 + m2x2 != Mx2 or abs(Mx2) > Jx2 or Jx2 > j1x2 + j2x2 or Jx2 < abs(j1x2 - j2x2):
            _CG_CACHE[key] = 0.0
        else:
            _CG_CACHE[key] = float(CG(S(j1x2) / 2, S(m1x2) / 2, S(j2x2) / 2, S(m2x2) / 2, S(Jx2) / 2,
                                      S(Mx2) / 2).doit())
    return _CG_CACHE[key]


def read_snt(path):
    zero, rows = 0.0, []
    for line in open(path):
        if "zero body" in line.lower():
            zero = float(line.split(":")[-1])
        body = line.split("!")[0].split("#")[0].split()
        if body:
            rows.append(body)
    nprot, nneut = int(rows[0][0]), int(rows[0][1])
    orbits = {int(r[0]): (int(r[1]), int(r[2]), int(r[3]), int(r[4])) for r in rows[1:1 + nprot + nneut]}
    i = 1 + nprot + nneut
    n1 = int(rows[i][0])
    one = {(int(r[0]), int(r[1])): float(r[2]) for r in rows[i + 1:i + 1 + n1]}
    i += 1 + n1
    n2 = int(rows[i][0])
    two = [(int(r[0]), int(r[1]), int(r[2]), int(r[3]), int(r[4]), float(r[5])) for r in rows[i + 1:i + 1 + n2]]
    return {"zero": zero, "orbits": orbits, "one": one, "two": two}


def read_3b(path, snt_orbits=None):
    """Entries (a, b, c, d, e, f, Jab, Jde, 2J, V); with snt_orbits the file's orbit indices are mapped to the
    .snt indices through (n, l, 2j, 2tz) from the file header ("!  idx  n l j/2 tz/2")."""
    entries, table = [], {}
    for line in open(path):
        if line.lstrip().startswith("!"):
            f = line.lstrip("! ").split()
            if len(f) == 5 and f[0].isdigit() and "/" in f[3]:
                table[int(f[0])] = (int(f[1]), int(f[2]), int(f[3].split("/")[0]), int(f[4].split("/")[0]))
            continue
        if not line.strip():
            continue
        f = line.split()
        entries.append((*(int(x) for x in f[:9]), float(f[9])))
    if snt_orbits is not None:
        to_snt = {lab: k for k, lab in snt_orbits.items()}
        remap = {k: to_snt[lab] for k, lab in table.items()}
        entries = [(*(remap[x] for x in e[:6]), *e[6:]) for e in entries]
    return entries


class Space:
    """Valence m-states: index -> (orbit, 2m); orbits from the .snt (index -> (n, l, 2j, 2tz))."""

    def __init__(self, orbits):
        self.orbits = orbits
        self.mstates = [(o, m) for o in sorted(orbits) for m in range(-orbits[o][2], orbits[o][2] + 1, 2)]
        self.index = {s: k for k, s in enumerate(self.mstates)}
        self.n = len(self.mstates)

    def j2(self, o):
        return self.orbits[o][2]

    def coupled_rows(self, orbs, labels, M):
        """Rows: unnormalized coupled states (labels: (Jab..., 2J) chains) on canonical determinants of orbs."""
        js = [self.j2(o) for o in orbs]
        dets = {}
        rows = []
        for lab in labels:
            row = {}
            for ms in _m_tuples(js, M):
                c = _chain_cg(js, ms, lab)
                if c == 0.0:
                    continue
                g = [self.index[(o, m)] for o, m in zip(orbs, ms)]
                if len(set(g)) < len(g):
                    continue
                order = sorted(range(len(g)), key=lambda k: g[k])
                sign = _perm_sign(order)
                det = tuple(sorted(g))
                row[det] = row.get(det, 0.0) + sign * c
            rows.append(row)
            for d in row:
                dets.setdefault(d, len(dets))
        U = np.zeros((len(labels), len(dets)))
        for r, row in enumerate(rows):
            for d, c in row.items():
                U[r, dets[d]] = c
        return U, list(dets)


def _m_tuples(js, M):
    if len(js) == 1:
        return [(M,)] if abs(M) <= js[0] and (M - js[0]) % 2 == 0 else []
    out = []
    for m in range(-js[0], js[0] + 1, 2):
        out += [(m, *rest) for rest in _m_tuples(js[1:], M - m)]
    return out


def _chain_cg(js, ms, lab):
    """Two-body lab = (2J,): <j1 m1 j2 m2|J M>.  Three-body lab = (Jab, 2J): <j1 m1 j2 m2|Jab Mab><Jab Mab j3 m3|J M>."""
    if len(js) == 2:
        return cg(js[0], ms[0], js[1], ms[1], lab[0], ms[0] + ms[1])
    Mab = ms[0] + ms[1]
    return cg(js[0], ms[0], js[1], ms[1], 2 * lab[0], Mab) * cg(2 * lab[0], Mab, js[2], ms[2], lab[1], Mab + ms[2])


def _perm_sign(order):
    sign, seen = 1, [False] * len(order)
    for i in range(len(order)):
        if seen[i]:
            continue
        j, length = i, 0
        while not seen[j]:
            seen[j] = True
            j = order[j]
            length += 1
        if length % 2 == 0:
            sign = -sign
    return sign


def mscheme_from_coupled(space, blocks, rank):
    """blocks: {(bra orbits, ket orbits): {(bra label, ket label): V_unnormalized}} with labels ending in 2J.
    Returns {(bra det, ket det): element} of the antisymmetric m-scheme operator."""
    out = {}
    for (ob, ok), elems in blocks.items():
        twoJs = sorted({lb[-1] for lb, _ in elems})
        for twoJ in twoJs:
            lb_list = sorted({lb for lb, _ in elems if lb[-1] == twoJ})
            lk_list = sorted({lk for _, lk in elems if lk[-1] == twoJ})
            V = np.array([[elems.get((lb, lk), 0.0) for lk in lk_list] for lb in lb_list])
            for M in range(-twoJ, twoJ + 1, 2):
                Ub, db = space.coupled_rows(ob, lb_list, M)
                Uk, dk = space.coupled_rows(ok, lk_list, M)
                if not db or not dk:
                    continue
                # V = Ub W Uk^T on this J block; W restricted to states of total J is recovered by the
                # pseudo-inverses (the coupled states with this J span exactly that part of W).
                W = np.linalg.pinv(Ub) @ V @ np.linalg.pinv(Uk).T
                for i, d in enumerate(db):
                    for k, e in enumerate(dk):
                        if abs(W[i, k]) > 1e-12:
                            out[(d, e)] = out.get((d, e), 0.0) + W[i, k]
    return out


def two_body_blocks(snt, space):
    blocks = {}
    for a, b, c, d, J, v in snt["two"]:
        nab = sqrt((1 + (a == b)) * (1 + (c == d)))  # normalized -> unnormalized coupled element
        for (p, q, r, s, val) in ((a, b, c, d, v), (c, d, a, b, v)):
            blocks.setdefault(((p, q), (r, s)), {})[((2 * J,), (2 * J,))] = nab * val
    return blocks


def three_body_blocks(entries):
    blocks = {}
    for a, b, c, d, e, f, Jab, Jde, twoJ, v in entries:
        blocks.setdefault(((a, b, c), (d, e, f)), {})[((Jab, twoJ), (Jde, twoJ))] = v
    return blocks


class Basis:
    def __init__(self, space, Zv, Nv, M2=None):
        prot = [k for k, (o, m) in enumerate(space.mstates) if space.orbits[o][3] < 0]
        neut = [k for k, (o, m) in enumerate(space.mstates) if space.orbits[o][3] > 0]
        dets = []
        for P in combinations(prot, Zv):
            for Nn in combinations(neut, Nv):
                occ = P + Nn
                if M2 is not None and sum(space.mstates[k][1] for k in occ) != M2:
                    continue
                dets.append(sum(1 << k for k in occ))
        self.dets = dets
        self.lookup = {d: i for i, d in enumerate(dets)}
        self.space = space


def _apply(det, ann, cre):
    """a+_{cre} ... a_{ann} (ann applied right to left in the order given) on bit string det; returns (sign, det)."""
    sign = 1
    for k in ann:
        if not det >> k & 1:
            return 0, 0
        sign *= -1 if bin(det & ((1 << k) - 1)).count("1") % 2 else 1
        det ^= 1 << k
    for k in reversed(cre):
        if det >> k & 1:
            return 0, 0
        sign *= -1 if bin(det & ((1 << k) - 1)).count("1") % 2 else 1
        det |= 1 << k
    return sign, det


def many_body_matrix(basis, one=None, two=None, three=None, zero=0.0):
    """one: {(i, k): f}; two: {((i,j), (k,l)): V} with i<j, k<l; three: {((i,j,k), (l,m,n)): W}, i<j<k, l<m<n."""
    dim = len(basis.dets)
    H = np.eye(dim) * zero
    terms = []
    for ops in (one, two, three):
        if ops:
            terms.append(ops)
    for col, det in enumerate(basis.dets):
        for ops in terms:
            for (bra, ket), val in ops.items():
                bra = bra if isinstance(bra, tuple) else (bra,)
                ket = ket if isinstance(ket, tuple) else (ket,)
                # a+_bra[0] a+_bra[1] .. a_ket[-1] .. a_ket[0]: a_ket[0] acts first
                sign, new = _apply(det, list(ket), list(bra))
                if sign == 0:
                    continue
                row = basis.lookup.get(new)
                if row is not None:
                    H[row, col] += sign * val
    return H


def one_body_mscheme(snt, space):
    out = {}
    for (a, b), v in snt["one"].items():
        for (o, m), k in space.index.items():
            if o == a and (b, m) in space.index:
                out[(k, space.index[(b, m)])] = v
                if a != b:
                    out[(space.index[(b, m)], k)] = v
    return out


def hamiltonian(snt, space, three_entries=None):
    one = one_body_mscheme(snt, space)
    two = mscheme_from_coupled(space, two_body_blocks(snt, space), 2)
    three = mscheme_from_coupled(space, three_body_blocks(three_entries), 3) if three_entries else None
    return one, two, three


def jsquared(space, Zv, Nv, basis0):
    """<J^2> matrix on an M=0 basis: J^2 = J- J+ there."""
    basis1 = Basis(space, Zv, Nv, 2)
    jp = {}
    for (o, m), k in space.index.items():
        if (o, m + 2) in space.index:
            j = space.orbits[o][2] / 2
            mm = m / 2
            jp[(space.index[(o, m + 2)], k)] = sqrt(j * (j + 1) - mm * (mm + 1))
    Jp = np.zeros((len(basis1.dets), len(basis0.dets)))
    for col, det in enumerate(basis0.dets):
        for (i, k), c in jp.items():
            sign, new = _apply(det, [k], [i])
            if sign and new in basis1.lookup:
                Jp[basis1.lookup[new], col] += sign * c
    return Jp.T @ Jp


def spectrum(snt_path, Zv, Nv, three_path=None, nstates=4, operator_paths=()):
    snt = read_snt(snt_path)
    space = Space(snt["orbits"])
    basis = Basis(space, Zv, Nv, 0)
    one, two, three = hamiltonian(snt, space, read_3b(three_path, snt["orbits"]) if three_path else None)
    H = many_body_matrix(basis, one, two, three, snt["zero"])
    if np.max(np.abs(H - H.T)) > 1e-6:
        raise ValueError(f"Hamiltonian is not symmetric: {np.max(np.abs(H - H.T))}")
    E, X = np.linalg.eigh(0.5 * (H + H.T))
    J2 = jsquared(space, Zv, Nv, basis)
    states = []
    for k in range(min(nstates, len(E))):
        jj = float(X[:, k] @ J2 @ X[:, k])
        states.append({"E": float(E[k]), "J": 0.5 * (-1 + sqrt(1 + 4 * jj))})
    for name, path in operator_paths:
        op = read_snt(path)
        o1, o2, _ = hamiltonian(op, space)
        O = many_body_matrix(basis, o1, o2, None, op["zero"])
        for k, s in enumerate(states):
            s[name] = float(X[:, k] @ O @ X[:, k])
    return states, (space, basis, X, H)


def _canon(idx):
    order = sorted(range(len(idx)), key=lambda k: idx[k])
    return tuple(sorted(idx)), _perm_sign(order)


def lower_body_from_three(three, occ):
    """Vacuum (core) form of :W:_ref for a valence three-body W and valence occupations occ[m-state]:
    two-body -sum_c n_c W_{abc,dec}, one-body +1/2 sum n n W, zero-body -1/6 sum n n n W (antisymmetric elements)."""
    W = three
    occupied = [k for k, n in enumerate(occ) if n > 0]
    two, one, zero = {}, {}, 0.0
    for (bra, ket), w in W.items():
        for c in occupied:
            if c in bra and c in ket:
                i = [x for x in bra if x != c]
                k = [x for x in ket if x != c]
                # W_{i1 i2 c, k1 k2 c} with the canonical triples: sign to move c to the last slot
                sb = (-1) ** (2 - bra.index(c))
                sk = (-1) ** (2 - ket.index(c))
                key = (tuple(i), tuple(k))
                two[key] = two.get(key, 0.0) - occ[c] * sb * sk * w
        for c1, c2 in combinations(occupied, 2):
            if c1 in bra and c2 in bra and c1 in ket and c2 in ket:
                i = [x for x in bra if x not in (c1, c2)][0]
                k = [x for x in ket if x not in (c1, c2)][0]
                sb = _canon((i, c1, c2))[1]
                sk = _canon((k, c1, c2))[1]
                one[(i, k)] = one.get((i, k), 0.0) + occ[c1] * occ[c2] * sb * sk * w
        if bra == ket and all(occ[x] > 0 for x in bra):
            zero -= occ[bra[0]] * occ[bra[1]] * occ[bra[2]] * w
    return zero, one, two
