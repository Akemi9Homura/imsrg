#!/usr/bin/env python3
"""BIGSTICK NCSM reference -> mrimsrg_jref (MR-IMSRG reference input).

BIGSTICK's own two-body densities (menu '2d') fail in Nmax-truncated spaces
("Problem cannot find sector jumps PP"; they work without truncation), so the
densities are computed here from the saved wavefunction (.wfn) and its ASCII
basis file (.bas, menu 'ba' with the same space and truncation).

For a J = 0 state only scalar densities are needed:
  gamma1_ab   = <a^dagger_{a m} a_{b m}>                (any m, averaged)
  gamma2^J    = <A^dagger_{ab;JM} A_{cd;JM}>             (normalized pairs, M-averaged)
built from one- and two-annihilation amplitudes (gamma2 = A2^T A2) restricted
to the rank-0 (M-averaged) projection.  Natural orbits diagonalize gamma1 within (l, j, tz) blocks
(most occupied orbit -> lowest n slot, as prepare_natural_basis), gamma2 is
rotated in the pair basis, lambda2 = gamma2 - gamma1^gamma1 (= n_a n_b on the
pair diagonal in natural orbits), and the file is written with the same
writer as export_jref.py.  No m-scheme two-body array is formed.

With --basis-transformation the BIGSTICK run is taken to be in another
single-particle basis (e.g. the imsrg++ HF basis, whose no2bpack was given to
BIGSTICK's `no2` reader): C.mat is the imsrg++ HO -> basis matrix (written as
<intfile>C.mat) and --basis-orbits a no2bpack whose header lists the imsrg++
orbit order.  The jref then stores C_total = C_basis C_nat (HO -> reference
natural orbits), so the MR driver keeps reading HO-basis interactions.

By default the three-body cumulant lambda3 (lambda3.py, from the
three-annihilation amplitudes on the fractionally occupied orbits) and the
BIGSTICK eigenvalue are added (jref version 2); --no-lambda3 writes the
version-1 file.

Usage:
  bigstick_jref.py --bas X.bas --wfn X.wfn --interaction H.MFDn \
      --nrefmax 2 --hw 20 [--e2max E] [--state 1] [--no-lambda3] \
      [--basis-transformation C.mat --basis-orbits H_basis.no2bpack] --output X.jref
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
from math import sqrt
from pathlib import Path
import struct

import numpy as np
from scipy import sparse

try:
    from .export_jref import write_jref
    from .jcoupling import _cg
    from .lambda3 import annihilation_matrix, coupled_lambda3
except ImportError:
    from export_jref import write_jref
    from jcoupling import _cg
    from lambda3 import annihilation_matrix, coupled_lambda3


# ----------------------------------------------------------------------------
# BIGSTICK files (format version 1, as shell-model-obs bigstick_reader.cpp)

class BigstickBasis:
    """ASCII .bas: orbits, hsps (m-states) and proton/neutron sectors of SDs."""

    def __init__(self, path: Path):
        tokens = path.read_text().split()
        pos = 0

        def take(count=1):
            nonlocal pos
            values = tokens[pos:pos + count]
            if len(values) != count:
                raise ValueError(f"truncated .bas file {path}")
            pos += count
            return values if count > 1 else values[0]

        if int(take()) != 31415926 or int(take()) != 1:
            raise ValueError(f"{path}: not a version-1 BIGSTICK .bas file")
        self.np = (int(take()), int(take()))
        take()  # isospin flag
        numorb = (int(take()), int(take()))
        self.orbits = []  # (n, l, 2j, 2tz): protons first, BIGSTICK order
        for species, tz2 in ((0, -1), (1, +1)):
            for _ in range(numorb[species]):
                nr, j2, l, _par, _w = (int(x) for x in take(5))
                self.orbits.append((nr, l, j2, tz2))
        nsps = (int(take()), int(take()))
        take(8 * (nsps[0] + nsps[1]))
        self.jz2 = int(take())
        self.parity = take()
        self.max_w = int(take())
        take(3)
        self.dimension = int(take())
        nhsps = [int(x) for x in take(5)]
        proton_max = nhsps[1] + nhsps[3]
        self.hsps = {}  # nstate -> (n, l, 2j, 2m, 2tz)
        for _ in range(sum(nhsps)):
            nstate, nr, l, j2, m2, _w, _par = (int(x) for x in take(7))
            self.hsps[nstate] = (nr, l, j2, m2, -1 if nstate <= proton_max else +1)
        self.sectors = []
        for species in (0, 1):
            nsectors, _nsd = int(take()), int(take())
            sectors = {}
            for _ in range(nsectors):
                sector_id = int(take())
                take(3)  # 2Jz, parity, W
                take(2)  # first/last SD
                nsd = int(take())
                ncs = int(take())
                conjugates = [int(x) for x in take(ncs)] if ncs > 1 else ([int(take())] if ncs == 1 else [])
                take(2)  # basis start/end
                sds = []
                for _ in range(nsd):
                    _ip, offset = int(take()), int(take())
                    occupied = [int(x) for x in take(self.np[species])] if self.np[species] > 1 else (
                        [int(take())] if self.np[species] == 1 else [])
                    sds.append((offset, occupied))
                sectors[sector_id] = (conjugates, sds)
            self.sectors.append(sectors)
        if pos != len(tokens):
            raise ValueError(f"{path}: trailing data after the sectors")


def read_bigstick_wfn(path: Path, basis: BigstickBasis, state: int):
    """Return (energy, J, amplitudes in BIGSTICK basis order) of one state."""
    data = path.read_bytes()
    pos = 0

    def ints(fmt):
        nonlocal pos
        values = struct.unpack_from("<" + fmt, data, pos)
        pos += struct.calcsize("<" + fmt)
        return values

    magic, version, _zero, vector_offset = ints("4i")
    if magic != 31415926 or version != 1:
        raise ValueError(f"{path}: not a version-1 BIGSTICK .wfn file")
    ints("4i")
    if ints("2i") != basis.np:
        raise ValueError(f"{path}: particle numbers differ from the .bas file")
    ints("i")
    numorb = ints("2i")
    ints(f"{5 * (numorb[0] + numorb[1])}i")
    nsps = ints("2i")
    ints(f"{8 * (nsps[0] + nsps[1])}i")
    jz2, _parity, _maxw = ints("3i")
    ints("3i")
    (dimension,) = ints("q")
    (nkeep,) = ints("i")
    if pos != vector_offset or dimension != basis.dimension or jz2 != basis.jz2:
        raise ValueError(f"{path}: header does not match the .bas file")
    remaining = len(data) - pos
    precision = None
    for bytes_per_amplitude in (4, 8):
        if remaining == nkeep * (20 + dimension * bytes_per_amplitude):
            precision = bytes_per_amplitude
    if precision is None:
        raise ValueError(f"{path}: size matches neither real4 nor real8 Lanczos vectors")
    if not 1 <= state <= nkeep:
        raise ValueError(f"{path}: state {state} not among the {nkeep} stored states")
    record = 20 + dimension * precision
    pos += (state - 1) * record
    vmagic, index = ints("2i")
    energy, xj, _xt2 = struct.unpack_from("<3f", data, pos)
    pos += 12
    if vmagic != 27182818 or index != state:
        raise ValueError(f"{path}: corrupt eigenvector record {state}")
    amplitudes = np.frombuffer(data, dtype="<f4" if precision == 4 else "<f8",
                               count=dimension, offset=pos).astype(np.float64)
    return float(energy), float(xj), amplitudes, precision


def configurations(basis: BigstickBasis):
    """Occupations (dimension x m-states) and Fermi signs, BIGSTICK order.

    m-states are ordered by orbit (BIGSTICK order) and then 2m; determinants
    are created in increasing m-state index.  BIGSTICK creates protons, then
    neutrons, each in increasing nstate; the permutation sign is absorbed.
    """
    mstates = []
    for orbit, (n, l, j2, tz2) in enumerate(basis.orbits):
        for m2 in range(-j2, j2 + 1, 2):
            mstates.append((orbit, m2))
    index = {}
    for k, (orbit, m2) in enumerate(mstates):
        n, l, j2, tz2 = basis.orbits[orbit]
        index[(n, l, j2, m2, tz2)] = k
    to_m = {nstate: index[label] for nstate, label in basis.hsps.items()}
    occupation = np.zeros((basis.dimension, len(mstates)), dtype=np.uint8)
    sign = np.zeros(basis.dimension, dtype=np.int8)

    def permutation_sign(sequence):
        inversions = sum(1 for i in range(len(sequence)) for j in range(i + 1, len(sequence))
                         if sequence[i] > sequence[j])
        return -1 if inversions % 2 else 1

    protons, neutrons = basis.sectors
    for _sector_id, (conjugates, proton_sds) in protons.items():
        for neutron_sector in conjugates:
            for neutron_offset, neutron_occupied in neutrons[neutron_sector][1]:
                nbits = [to_m[x] for x in neutron_occupied]
                for proton_offset, proton_occupied in proton_sds:
                    pbits = [to_m[x] for x in proton_occupied]
                    row = proton_offset + neutron_offset - 1
                    if sign[row] != 0:
                        raise ValueError("duplicate configuration index in .bas")
                    occupation[row, pbits + nbits] = 1
                    sign[row] = permutation_sign(pbits + nbits)
    if np.any(sign == 0):
        raise ValueError(".bas sectors do not generate every configuration")
    return mstates, occupation, sign


# ----------------------------------------------------------------------------
# Scalar densities

def scalar_densities(orbits, mstates, occupation, coefficients):
    """gamma1[a, b] (orbit level) and gamma2^J blocks for normalized pairs."""
    norbits = len(orbits)
    by_orbit = defaultdict(dict)
    for k, (orbit, m2) in enumerate(mstates):
        by_orbit[orbit][m2] = k

    A1 = annihilation_matrix(occupation, coefficients, [(k,) for k in range(len(mstates))])
    g1m = (A1.T @ A1).toarray()
    gamma1 = np.zeros((norbits, norbits))
    for a in range(norbits):
        for b in range(norbits):
            na, la, ja, ta = orbits[a]
            nb, lb, jb, tb = orbits[b]
            if (la, ja, ta) == (lb, jb, tb):
                gamma1[a, b] = np.mean([g1m[by_orbit[a][m], by_orbit[b][m]] for m in range(-ja, ja + 1, 2)])

    # gamma2^J = (2J+1)^-1 sum_M (A2 P_M)^T (A2 P_M): the rank-0 projection
    # (all M), so the contraction identities with gamma1 hold exactly.
    pairs = [(p, q) for p in range(len(mstates)) for q in range(p + 1, len(mstates))]
    pair_index = {pq: k for k, pq in enumerate(pairs)}
    A2 = annihilation_matrix(occupation, coefficients, pairs)

    channels = {}
    for a in range(norbits):
        for b in range(a, norbits):
            _, la, ja, ta = orbits[a]
            _, lb, jb, tb = orbits[b]
            for J in range(abs(ja - jb) // 2, (ja + jb) // 2 + 1):
                if a == b and J % 2:
                    continue
                channels.setdefault((J, (la + lb) % 2, (ta + tb) // 2), []).append((a, b))
    blocks = {}
    for (J, parity, Tz), pair_list in sorted(channels.items()):
        matrix = np.zeros((len(pair_list), len(pair_list)))
        for M2 in range(-2 * J, 2 * J + 1, 2):
            rows, cols, values = [], [], []
            for column, (a, b) in enumerate(pair_list):
                ja, jb = orbits[a][2], orbits[b][2]
                norm = sqrt(2.0) if a == b else 1.0
                for ma in range(-ja, ja + 1, 2):
                    mb = M2 - ma
                    if abs(mb) > jb:
                        continue
                    c = _cg(ja, jb, J, ma, mb)
                    if c == 0.0:
                        continue
                    p, q = by_orbit[a][ma], by_orbit[b][mb]
                    if p == q:
                        continue
                    sign = 1.0
                    if p > q:
                        p, q, sign = q, p, -1.0
                    rows.append(pair_index[(p, q)])
                    cols.append(column)
                    values.append(sign * c / norm)
            P = sparse.csr_matrix((values, (rows, cols)), shape=(len(pairs), len(pair_list)))
            B = (A2 @ P).toarray() if A2.shape[0] * len(pair_list) < 5e7 else (A2 @ P)
            matrix += np.asarray(B.T @ B if isinstance(B, np.ndarray) else (B.T @ B).toarray())
        blocks[(J, parity, Tz)] = (pair_list, matrix / (2 * J + 1))
    return gamma1, blocks


# ----------------------------------------------------------------------------
# Natural orbits, lambda2

def natural_orbits(orbits, gamma1, tolerance=1e-10):
    """C[old, new]: per (l, j, tz) block, connected components of gamma1,
    occupations descending onto slots of ascending n, positive pivot."""
    norbits = len(orbits)
    C = np.eye(norbits)
    groups = defaultdict(list)
    for a, (n, l, j2, tz2) in enumerate(orbits):
        groups[(l, j2, tz2)].append(a)
    for members in groups.values():
        members = sorted(members, key=lambda a: orbits[a][0])
        remaining = set(members)
        while remaining:
            start = min(remaining, key=lambda a: orbits[a][0])
            component, stack = [], [start]
            remaining.discard(start)
            while stack:
                a = stack.pop()
                component.append(a)
                for b in list(remaining):
                    if abs(gamma1[a, b]) > tolerance:
                        remaining.discard(b)
                        stack.append(b)
            if len(component) == 1:
                continue
            component.sort(key=lambda a: orbits[a][0])
            values, vectors = np.linalg.eigh(gamma1[np.ix_(component, component)])
            vectors = vectors[:, np.argsort(values)[::-1]]
            for column in range(vectors.shape[1]):
                pivot = int(np.argmax(np.abs(vectors[:, column])))
                if vectors[pivot, column] < 0.0:
                    vectors[:, column] *= -1.0
            C[np.ix_(component, component)] = vectors
    return C


def rotate_pairs(orbits, pair_list, J, C):
    """D[old pair, new pair] in the normalized antisymmetric pair basis."""
    D = np.zeros((len(pair_list), len(pair_list)))
    for i, (p, q) in enumerate(pair_list):
        phase = -1.0 if ((orbits[p][2] + orbits[q][2]) // 2 + J + 1) % 2 else 1.0
        for k, (r, s) in enumerate(pair_list):
            value = C[p, r] * C[q, s]
            if p != q:
                value += C[q, r] * C[p, s] * phase
            if p == q:
                value *= sqrt(2.0)
            if r == s:
                value /= sqrt(2.0)
            D[i, k] = value
    return D


def read_arma_binary(path: Path) -> np.ndarray:
    """Armadillo arma_binary double matrix (Mat::save default)."""
    data = Path(path).read_bytes()
    first = data.index(b"\n")
    second = data.index(b"\n", first + 1)
    if data[:first] != b"ARMA_MAT_BIN_FN008":
        raise ValueError(f"{path}: not an Armadillo binary double matrix")
    rows, cols = (int(x) for x in data[first + 1:second].split())
    values = np.frombuffer(data, dtype="<f8", count=rows * cols, offset=second + 1)
    return values.reshape((cols, rows)).T.copy()


def read_no2bpack_orbits(path: Path):
    """Orbit labels (n, l, 2j, 2tz) in file order from a no2bpack header."""
    with Path(path).open("rb") as stream:
        _hw, _emax, norbits, _nobme, _ntbme = struct.unpack("<d4i", stream.read(24))
        return [tuple(struct.unpack("<4i", stream.read(16))) for _ in range(norbits)]


def basis_transformation_for(orbits, c_path: Path, labels_path: Path) -> np.ndarray:
    """imsrg++ HO -> basis matrix reordered to the BIGSTICK orbit list."""
    matrix = read_arma_binary(c_path)
    labels = read_no2bpack_orbits(labels_path)
    if matrix.shape != (len(labels), len(labels)):
        raise ValueError("basis transformation and orbit table have different sizes")
    if sorted(labels) != sorted(orbits):
        raise ValueError("the BIGSTICK orbit space must equal the full imsrg++ orbit space "
                         "(use auto with the emax of the basis Hamiltonian)")
    index = {label: k for k, label in enumerate(labels)}
    order = [index[label] for label in orbits]
    reordered = matrix[np.ix_(order, order)]
    if np.max(np.abs(reordered.T @ reordered - np.eye(len(orbits)))) > 1e-10:
        raise ValueError("basis transformation is not orthogonal")
    return reordered


def sha256_files(*paths):
    digest = hashlib.sha256()
    for path in paths:
        digest.update(Path(path).read_bytes())
    return digest.hexdigest()


def convert(bas_path, wfn_path, interaction_path, output, *, nrefmax, hw, e2max=None,
            state=1, tolerance=1e-10, lambda3=True, basis_transformation=None, basis_orbits=None):
    basis = BigstickBasis(Path(bas_path))
    energy, xj, amplitudes, precision = read_bigstick_wfn(Path(wfn_path), basis, state)
    mstates, occupation, sign = configurations(basis)
    coefficients = sign * amplitudes
    coefficients /= np.linalg.norm(coefficients)
    orbits = basis.orbits
    Z, N = basis.np
    A = Z + N
    J2 = int(round(2 * xj))
    if J2 != 0 or basis.parity != "+":
        raise ValueError(f"state {state} has 2J={2 * xj:.3f}, parity {basis.parity}; need J=0+")

    gamma1, blocks = scalar_densities(orbits, mstates, occupation, coefficients)
    j2 = np.array([o[2] for o in orbits])
    proton = np.array([o[3] == -1 for o in orbits])
    traces = ((j2 + 1) * np.diag(gamma1))
    pair_trace = sum((2 * J + 1) * np.trace(matrix) for (J, _, _), (_, matrix) in blocks.items())
    checks = {
        "Z": float(traces[proton].sum()), "N": float(traces[~proton].sum()),
        "pairs": float(pair_trace), "gamma1 asymmetry": float(np.max(np.abs(gamma1 - gamma1.T))),
    }
    if abs(checks["Z"] - Z) > 1e-6 or abs(checks["N"] - N) > 1e-6 or abs(pair_trace - A * (A - 1) / 2) > 1e-5:
        raise ValueError(f"density traces inconsistent: {checks}")
    gamma1 = 0.5 * (gamma1 + gamma1.T)

    C = natural_orbits(orbits, gamma1, tolerance)
    occupations = np.diag(C.T @ gamma1 @ C).copy()
    occupations[np.abs(occupations) < 1e-14] = 0.0
    lambda_blocks = []
    for (J, parity, Tz), (pair_list, matrix) in blocks.items():
        D = rotate_pairs(orbits, pair_list, J, C)
        natural = D.T @ (0.5 * (matrix + matrix.T)) @ D
        disconnected = np.diag([occupations[a] * occupations[b] for a, b in pair_list])
        lambda_blocks.append((J, parity, Tz, pair_list, natural - disconnected))

    lambda3_blocks = None
    reference_energy = None
    lambda3_report = {}
    if lambda3:
        lambda3_blocks, lambda3_report = coupled_lambda3(
            orbits, mstates, occupation, coefficients, C, occupations)
        reference_energy = (energy, f"BIGSTICK eigenvalue (float32), state {state}")

    transformation = C
    if basis_transformation is not None:
        transformation = basis_transformation_for(orbits, basis_transformation, basis_orbits) @ C

    emax = max(2 * n + l for n, l, _, _ in orbits)
    rdm_digest = hashlib.sha256(np.asarray(gamma1, dtype="<f8").tobytes())
    for _, _, _, _, matrix in lambda_blocks:
        rdm_digest.update(np.asarray(matrix, dtype="<f8").tobytes())
    write_jref(
        Path(output),
        header=(A, Z, nrefmax, J2, 1, emax, e2max if e2max is not None else 2 * emax, hw),
        orbits=[(k, n, l, j2_, tz2) for k, (n, l, j2_, tz2) in enumerate(orbits)],
        occupations=occupations,
        transformation=transformation,
        lambda_blocks=lambda_blocks,
        digests=(sha256_files(interaction_path), rdm_digest.hexdigest(),
                 sha256_files(bas_path, wfn_path)),
        lambda3_blocks=lambda3_blocks,
        reference_energy=reference_energy,
    )
    return {**lambda3_report, "energy": energy, "2J": 2 * xj, "precision": f"real{precision}", "dimension": basis.dimension,
            "A": A, "Z": Z, "Nrefmax": nrefmax, "emax": emax, **checks,
            "min occupation > 0": float(min(o for o in occupations if o > 0))}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bas", required=True, type=Path)
    parser.add_argument("--wfn", required=True, type=Path)
    parser.add_argument("--interaction", required=True, type=Path,
                        help="interaction file BIGSTICK used (its SHA-256 goes into the jref)")
    parser.add_argument("--nrefmax", required=True, type=int)
    parser.add_argument("--hw", required=True, type=float)
    parser.add_argument("--e2max", type=int)
    parser.add_argument("--state", type=int, default=1)
    parser.add_argument("--no-lambda3", action="store_true",
                        help="skip lambda3 and the reference energy (write jref version 1)")
    parser.add_argument("--basis-transformation", type=Path,
                        help="imsrg++ <intfile>C.mat: HO -> basis of the BIGSTICK Hamiltonian")
    parser.add_argument("--basis-orbits", type=Path,
                        help="no2bpack whose header gives the imsrg++ orbit order of C.mat")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if (args.basis_transformation is None) != (args.basis_orbits is None):
        parser.error("--basis-transformation and --basis-orbits go together")
    summary = convert(args.bas, args.wfn, args.interaction, args.output, nrefmax=args.nrefmax,
                      hw=args.hw, e2max=args.e2max, state=args.state, lambda3=not args.no_lambda3,
                      basis_transformation=args.basis_transformation, basis_orbits=args.basis_orbits)
    for key, value in summary.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
