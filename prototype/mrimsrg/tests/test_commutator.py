import sys
from itertools import combinations
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from commutator import commutator
from densities import compute_densities
from normal_order import MRHamiltonian, to_vacuum


def occupations(*occupied_sets: tuple[int, ...], norb: int) -> np.ndarray:
    result = np.zeros((len(occupied_sets), norb), dtype=np.uint8)
    for row, occupied in enumerate(occupied_sets):
        result[row, list(occupied)] = 1
    return result


def antisymmetrize_two_body(raw: np.ndarray) -> np.ndarray:
    return 0.25 * (
        raw
        - raw.swapaxes(0, 1)
        - raw.swapaxes(2, 3)
        + raw.transpose(1, 0, 3, 2)
    )


def random_operator(norb: int, seed: int, adjoint_sign: int | None) -> MRHamiltonian:
    rng = np.random.default_rng(seed)
    one_body = rng.normal(size=(norb, norb))
    two_body = antisymmetrize_two_body(rng.normal(size=(norb,) * 4))
    if adjoint_sign is not None:
        one_body = 0.5 * (one_body + adjoint_sign * one_body.T)
        two_body = 0.5 * (
            two_body + adjoint_sign * two_body.transpose(2, 3, 0, 1)
        )
    return MRHamiltonian(0.0, one_body, two_body)


def annihilate(det: int, orbit: int) -> tuple[int, int] | None:
    if not ((det >> orbit) & 1):
        return None
    phase = -1 if (det & ((1 << orbit) - 1)).bit_count() % 2 else 1
    return det ^ (1 << orbit), phase


def create(det: int, orbit: int) -> tuple[int, int] | None:
    if (det >> orbit) & 1:
        return None
    phase = -1 if (det & ((1 << orbit) - 1)).bit_count() % 2 else 1
    return det | (1 << orbit), phase


def vacuum_matrix(operator: MRHamiltonian, densities, particle_number: int) -> tuple[np.ndarray, list[int]]:
    vacuum = to_vacuum(operator, densities)
    norb = vacuum.one_body.shape[0]
    basis = [sum(1 << p for p in occupied) for occupied in combinations(range(norb), particle_number)]
    index = {det: i for i, det in enumerate(basis)}
    matrix = np.zeros((len(basis), len(basis)))
    matrix += np.eye(len(basis)) * vacuum.zero_body

    for ket_index, ket in enumerate(basis):
        for q in range(norb):
            state = annihilate(ket, q)
            if state is None:
                continue
            reduced, phase_q = state
            for p in range(norb):
                final = create(reduced, p)
                if final is not None:
                    matrix[index[final[0]], ket_index] += (
                        vacuum.one_body[p, q] * phase_q * final[1]
                    )
        for r in range(norb):
            state_r = annihilate(ket, r)
            if state_r is None:
                continue
            for s in range(norb):
                state_s = annihilate(state_r[0], s)
                if state_s is None:
                    continue
                for q in range(norb):
                    state_q = create(state_s[0], q)
                    if state_q is None:
                        continue
                    for p in range(norb):
                        final = create(state_q[0], p)
                        if final is not None:
                            phase = state_r[1] * state_s[1] * state_q[1] * final[1]
                            matrix[index[final[0]], ket_index] += (
                                0.25 * vacuum.two_body[p, q, r, s] * phase
                            )
    return matrix, basis


class CommutatorTests(unittest.TestCase):
    def test_antisymmetry_for_correlated_reference(self) -> None:
        norb = 4
        densities = compute_densities(
            occupations((0, 1), (2, 3), norb=norb),
            np.array([1.0, -1.0]) / np.sqrt(2.0),
        )
        left = random_operator(norb, 19, None)
        right = random_operator(norb, 83, None)
        xy = commutator(left, right, densities)
        yx = commutator(right, left, densities)
        self.assertAlmostEqual(xy.zero_body, -yx.zero_body, places=12)
        np.testing.assert_allclose(xy.one_body, -yx.one_body, atol=2e-12)
        np.testing.assert_allclose(xy.two_body, -yx.two_body, atol=2e-12)

    def test_antihermitian_hermitian_commutator_is_hermitian(self) -> None:
        norb = 4
        densities = compute_densities(
            occupations((0, 1), (2, 3), norb=norb),
            np.array([1.0, 1.0]) / np.sqrt(2.0),
        )
        eta = random_operator(norb, 117, -1)
        hamiltonian = random_operator(norb, 215, 1)
        result = commutator(eta, hamiltonian, densities)
        np.testing.assert_allclose(result.one_body, result.one_body.T, atol=2e-12)
        np.testing.assert_allclose(
            result.two_body, result.two_body.transpose(2, 3, 0, 1), atol=2e-12
        )
        np.testing.assert_allclose(
            result.two_body, -result.two_body.swapaxes(0, 1), atol=2e-12
        )
        np.testing.assert_allclose(
            result.two_body, -result.two_body.swapaxes(2, 3), atol=2e-12
        )

    def test_correlated_zero_body_matches_explicit_fock_commutator(self) -> None:
        norb = 4
        coefficients = np.array([np.sqrt(0.3), -np.sqrt(0.7)])
        determinant_array = occupations((0, 1), (2, 3), norb=norb)
        densities = compute_densities(determinant_array, coefficients)
        left = random_operator(norb, 907, None)
        # Restrict the left operator to 1B.  The exact [1B,1B+2B]
        # expectation depends only on gamma1/lambda2 and is therefore a valid
        # oracle for the lambda3=0 equations.  A general [2B,2B] expectation
        # also contains lambda3, which need not vanish even when A=2.
        left = MRHamiltonian(left.zero_body, left.one_body, np.zeros_like(left.two_body))
        right = random_operator(norb, 1117, None)
        result = commutator(left, right, densities)
        left_matrix, basis = vacuum_matrix(left, densities, 2)
        right_matrix, _ = vacuum_matrix(right, densities, 2)
        reference = np.zeros(len(basis))
        for occupied, coefficient in zip(((0, 1), (2, 3)), coefficients, strict=True):
            det = sum(1 << p for p in occupied)
            reference[basis.index(det)] = coefficient
        expected = reference @ (left_matrix @ right_matrix - right_matrix @ left_matrix) @ reference
        self.assertAlmostEqual(result.zero_body, expected, places=11)

    def test_single_reference_limit_matches_explicit_single_double_excitations(self) -> None:
        norb = 5
        occupied = (0, 2)
        determinant_array = occupations(occupied, norb=norb)
        densities = compute_densities(determinant_array, np.array([1.0]))
        np.testing.assert_allclose(densities.lambda2, 0.0, atol=1e-14)
        left = random_operator(norb, 2718, None)
        right = random_operator(norb, 3141, None)
        result = commutator(left, right, densities)
        left_matrix, basis = vacuum_matrix(left, densities, len(occupied))
        right_matrix, _ = vacuum_matrix(right, densities, len(occupied))
        exact = left_matrix @ right_matrix - right_matrix @ left_matrix
        reference_det = sum(1 << p for p in occupied)
        reference_index = basis.index(reference_det)
        self.assertAlmostEqual(result.zero_body, exact[reference_index, reference_index], places=11)

        particles = tuple(p for p in range(norb) if p not in occupied)
        for hole in occupied:
            state_h = annihilate(reference_det, hole)
            for particle in particles:
                state_p = create(state_h[0], particle)
                excitation_index = basis.index(state_p[0])
                phase = state_h[1] * state_p[1]
                self.assertAlmostEqual(
                    phase * result.one_body[particle, hole],
                    exact[excitation_index, reference_index],
                    places=11,
                )

        for holes in combinations(occupied, 2):
            state = (reference_det, 1)
            for hole in holes:
                local = annihilate(state[0], hole)
                state = (local[0], state[1] * local[1])
            for particles_pair in combinations(particles, 2):
                created = state
                for particle in reversed(particles_pair):
                    local = create(created[0], particle)
                    created = (local[0], created[1] * local[1])
                excitation_index = basis.index(created[0])
                p, q = particles_pair
                i, j = holes
                self.assertAlmostEqual(
                    created[1] * result.two_body[p, q, i, j],
                    exact[excitation_index, reference_index],
                    places=11,
                )

    def test_rejects_nonnatural_basis(self) -> None:
        norb = 3
        densities = compute_densities(
            occupations((0,), (1,), norb=norb), np.array([3.0, 4.0]) / 5.0
        )
        operator = random_operator(norb, 44, 1)
        with self.assertRaisesRegex(ValueError, "natural-orbital basis"):
            commutator(operator, operator, densities)


if __name__ == "__main__":
    unittest.main()
