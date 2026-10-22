from __future__ import annotations

import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from densities import compute_densities
from normal_order import MRHamiltonian, VacuumHamiltonian, normal_order, to_vacuum, validate_hermitian


def random_hamiltonian(norb: int, seed: int = 1927) -> VacuumHamiltonian:
    rng = np.random.default_rng(seed)
    one_body = rng.normal(size=(norb, norb))
    one_body = 0.5 * (one_body + one_body.T)
    raw = rng.normal(size=(norb, norb, norb, norb))
    two_body = raw - raw.swapaxes(0, 1) - raw.swapaxes(2, 3) + raw.transpose(1, 0, 3, 2)
    two_body = 0.5 * (two_body + two_body.transpose(2, 3, 0, 1))
    return VacuumHamiltonian(0.731, one_body, two_body)


def occupations(*occupied_sets: tuple[int, ...], norb: int) -> np.ndarray:
    result = np.zeros((len(occupied_sets), norb), dtype=np.uint8)
    for row, occupied in enumerate(occupied_sets):
        result[row, list(occupied)] = 1
    return result


def apply_annihilation(det: int, orbit: int) -> tuple[int, int] | None:
    if not ((det >> orbit) & 1):
        return None
    phase = -1 if bin(det & ((1 << orbit) - 1)).count("1") % 2 else 1
    return det ^ (1 << orbit), phase


def apply_creation(det: int, orbit: int) -> tuple[int, int] | None:
    if (det >> orbit) & 1:
        return None
    phase = -1 if bin(det & ((1 << orbit) - 1)).count("1") % 2 else 1
    return det | (1 << orbit), phase


def explicit_matrix_element(hamiltonian: VacuumHamiltonian, bra: int, ket: int, norb: int) -> float:
    value = hamiltonian.zero_body if bra == ket else 0.0
    for p in range(norb):
        for q in range(norb):
            state = apply_annihilation(ket, q)
            if state is None:
                continue
            intermediate, phase = state
            state = apply_creation(intermediate, p)
            if state is not None and state[0] == bra:
                value += hamiltonian.one_body[p, q] * phase * state[1]
    for p in range(norb):
        for q in range(norb):
            for r in range(norb):
                for s in range(norb):
                    state = apply_annihilation(ket, r)
                    if state is None:
                        continue
                    intermediate, phase = state
                    state = apply_annihilation(intermediate, s)
                    if state is None:
                        continue
                    intermediate, local_phase = state
                    phase *= local_phase
                    state = apply_creation(intermediate, q)
                    if state is None:
                        continue
                    intermediate, local_phase = state
                    phase *= local_phase
                    state = apply_creation(intermediate, p)
                    if state is not None and state[0] == bra:
                        value += 0.25 * hamiltonian.two_body[p, q, r, s] * phase * state[1]
    return float(value)


class NormalOrderTests(unittest.TestCase):
    def test_correlated_reference_round_trip(self) -> None:
        norb = 5
        determinants = occupations((0, 1), (0, 2), (3, 4), norb=norb)
        coefficients = np.array([0.5, -0.5, np.sqrt(0.5)])
        densities = compute_densities(determinants, coefficients)
        vacuum = random_hamiltonian(norb)
        mr = normal_order(vacuum, densities)
        recovered = to_vacuum(mr, densities)
        validate_hermitian(mr)
        np.testing.assert_allclose(recovered.zero_body, vacuum.zero_body, atol=2e-13)
        np.testing.assert_allclose(recovered.one_body, vacuum.one_body, atol=2e-13)
        np.testing.assert_allclose(recovered.two_body, vacuum.two_body, atol=2e-13)

    def test_zero_body_is_explicit_reference_expectation(self) -> None:
        norb = 4
        determinant = (1 << 0) | (1 << 2)
        densities = compute_densities(occupations((0, 2), norb=norb), np.array([1.0]))
        vacuum = random_hamiltonian(norb, seed=9)
        mr = normal_order(vacuum, densities)
        direct = explicit_matrix_element(vacuum, determinant, determinant, norb)
        self.assertAlmostEqual(mr.zero_body, direct, places=12)

    def test_one_body_index_order_matches_single_excitation(self) -> None:
        norb = 4
        reference_det = (1 << 0) | (1 << 1)
        densities = compute_densities(occupations((0, 1), norb=norb), np.array([1.0]))
        vacuum = random_hamiltonian(norb, seed=28)
        mr = normal_order(vacuum, densities)

        annihilated, phase_ann = apply_annihilation(reference_det, 0)
        excited_det, phase_cre = apply_creation(annihilated, 2)
        excitation_phase = phase_ann * phase_cre
        direct = explicit_matrix_element(vacuum, excited_det, reference_det, norb)
        self.assertAlmostEqual(direct, excitation_phase * mr.one_body[2, 0], places=12)

    def test_empty_reference_leaves_hamiltonian_unchanged(self) -> None:
        norb = 3
        densities = compute_densities(occupations((), norb=norb), np.array([1.0]))
        vacuum = random_hamiltonian(norb, seed=71)
        mr = normal_order(vacuum, densities)
        np.testing.assert_allclose(mr.one_body, vacuum.one_body)
        np.testing.assert_allclose(mr.two_body, vacuum.two_body)
        self.assertEqual(mr.zero_body, vacuum.zero_body)


if __name__ == "__main__":
    unittest.main()
