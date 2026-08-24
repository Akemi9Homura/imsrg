"""Exact one- and two-body densities from a small NCSM wavefunction.

Conventions are explicit throughout:

    gamma1[p,q]       = <Psi| a^dagger_p a_q |Psi>
    gamma2[p,q,r,s]   = <Psi| a^dagger_p a^dagger_q a_s a_r |Psi>
    lambda2           = gamma2 - A(gamma1 gamma1)

The annihilated-state Gram construction avoids an O(dim*norb**4) loop and is
still simple enough to verify directly in tiny Fock spaces.
"""

from __future__ import annotations

from dataclasses import dataclass
import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy import sparse

try:
    from .reference_io import ReferenceData, load_reference
except ImportError:  # Direct execution: python prototype/mrimsrg/densities.py ...
    from reference_io import ReferenceData, load_reference


@dataclass(frozen=True)
class Densities:
    gamma1: np.ndarray
    gamma2: np.ndarray
    lambda2: np.ndarray


def _occupied_orbits(det: int, norb: int) -> list[int]:
    return [p for p in range(norb) if (det >> p) & 1]


def _annihilate(det: int, orbit: int) -> tuple[int, int]:
    if not ((det >> orbit) & 1):
        raise ValueError("cannot annihilate an unoccupied orbital")
    lower_mask = (1 << orbit) - 1
    phase = -1 if (det & lower_mask).bit_count() % 2 else 1
    return det ^ (1 << orbit), phase


def _determinant_ints(occupations: np.ndarray) -> list[int]:
    norb = int(occupations.shape[1])
    if norb > 64:
        raise ValueError("the rapid prototype supports at most 64 m-scheme orbitals")
    weights = np.left_shift(np.uint64(1), np.arange(norb, dtype=np.uint64))
    return [int(value) for value in occupations.astype(np.uint64) @ weights]


def _annihilation_amplitudes(
    determinants: list[int], coefficients: np.ndarray, norb: int, rank: int
) -> tuple[sparse.csr_matrix, list[tuple[int, ...]]]:
    columns = list(combinations(range(norb), rank))
    column_index = {orbits: index for index, orbits in enumerate(columns)}
    reduced_index: dict[int, int] = {}
    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []

    for det, coefficient in zip(determinants, coefficients, strict=True):
        occupied = _occupied_orbits(det, norb)
        for orbits in combinations(occupied, rank):
            reduced = det
            phase = 1
            # a_{p_rank} ... a_{p_2} a_{p_1}: the smallest orbital acts first.
            for orbit in orbits:
                reduced, local_phase = _annihilate(reduced, orbit)
                phase *= local_phase
            row = reduced_index.setdefault(reduced, len(reduced_index))
            rows.append(row)
            cols.append(column_index[orbits])
            values.append(float(coefficient) * phase)

    matrix = sparse.coo_matrix(
        (values, (rows, cols)), shape=(len(reduced_index), len(columns)), dtype=np.float64
    ).tocsr()
    return matrix, columns


def compute_densities(determinants: np.ndarray, coefficients: np.ndarray) -> Densities:
    norb = int(determinants.shape[1])
    det_ints = _determinant_ints(determinants)

    one_amplitudes, one_columns = _annihilation_amplitudes(det_ints, coefficients, norb, 1)
    if one_columns != [(p,) for p in range(norb)]:
        raise AssertionError("unexpected one-body column ordering")
    gamma1 = (one_amplitudes.T @ one_amplitudes).toarray()

    pair_amplitudes, pairs = _annihilation_amplitudes(det_ints, coefficients, norb, 2)
    pair_density = (pair_amplitudes.T @ pair_amplitudes).toarray()
    gamma2 = np.zeros((norb, norb, norb, norb), dtype=np.float64)
    pair_array = np.asarray(pairs, dtype=np.intp)
    p = pair_array[:, 0]
    q = pair_array[:, 1]
    gamma2[p[:, None], q[:, None], p[None, :], q[None, :]] = pair_density
    gamma2[q[:, None], p[:, None], p[None, :], q[None, :]] = -pair_density
    gamma2[p[:, None], q[:, None], q[None, :], p[None, :]] = -pair_density
    gamma2[q[:, None], p[:, None], q[None, :], p[None, :]] = pair_density

    disconnected = np.einsum("pr,qs->pqrs", gamma1, gamma1) - np.einsum(
        "ps,qr->pqrs", gamma1, gamma1
    )
    lambda2 = gamma2 - disconnected
    return Densities(gamma1=gamma1, gamma2=gamma2, lambda2=lambda2)


def validate_densities(densities: Densities, particle_number: int, tolerance: float = 2e-10) -> None:
    gamma1 = densities.gamma1
    gamma2 = densities.gamma2
    lambda2 = densities.lambda2
    if abs(float(np.trace(gamma1)) - particle_number) > tolerance:
        raise ValueError("Tr(gamma1) does not equal A")
    if np.max(np.abs(gamma1 - gamma1.T)) > tolerance:
        raise ValueError("gamma1 is not Hermitian")
    if np.max(np.abs(gamma2 + gamma2.swapaxes(0, 1))) > tolerance:
        raise ValueError("gamma2 is not antisymmetric in creation indices")
    if np.max(np.abs(gamma2 + gamma2.swapaxes(2, 3))) > tolerance:
        raise ValueError("gamma2 is not antisymmetric in annihilation indices")
    if np.max(np.abs(gamma2 - gamma2.transpose(2, 3, 0, 1))) > tolerance:
        raise ValueError("gamma2 is not Hermitian")
    contraction = np.einsum("pqrq->pr", gamma2)
    if np.max(np.abs(contraction - (particle_number - 1) * gamma1)) > tolerance:
        raise ValueError("gamma2 contraction does not reproduce (A-1) gamma1")
    if np.max(np.abs(lambda2 + lambda2.swapaxes(0, 1))) > tolerance:
        raise ValueError("lambda2 is not antisymmetric in creation indices")
    if np.max(np.abs(lambda2 + lambda2.swapaxes(2, 3))) > tolerance:
        raise ValueError("lambda2 is not antisymmetric in annihilation indices")


def reference_energy(reference: ReferenceData, densities: Densities) -> float:
    zero_body = float(reference.metadata["zero_body"])
    one_body = np.einsum("pq,pq->", reference.one_body, densities.gamma1, optimize=True)
    two_body = 0.25 * np.einsum(
        "pqrs,pqrs->", reference.two_body, densities.gamma2, optimize=True
    )
    return float(zero_body + one_body + two_body)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    args = parser.parse_args()
    reference = load_reference(args.reference)
    densities = compute_densities(reference.determinants, reference.coefficients)
    validate_densities(densities, int(reference.metadata["A"]))
    contracted_energy = reference_energy(reference, densities)
    expected_energy = float(reference.metadata["reference_energy"])
    if abs(contracted_energy - expected_energy) > 1e-8:
        raise ValueError(
            f"density contraction energy {contracted_energy:.12f} differs from NCSM {expected_energy:.12f}"
        )
    print(f"Tr(gamma1) = {np.trace(densities.gamma1):.12f}")
    print(f"max|lambda2| = {np.max(np.abs(densities.lambda2)):.12e}")
    print(f"Eref = {contracted_energy:.12f} MeV")


if __name__ == "__main__":
    main()
