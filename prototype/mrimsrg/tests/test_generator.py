import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from commutator import commutator
from densities import compute_densities
from generator import (
    _safe_denominator,
    decoupling_masks,
    decoupling_matrix_elements,
    epstein_nesbet_denominators,
    masked_decoupling_residual,
    oscillator_quanta_from_orbits,
    white_generator,
    white_ncsm_matrix_elements,
    white_ncsm_numerator_residual,
)
from normal_order import MRHamiltonian


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


def random_hermitian_hamiltonian(norb: int, seed: int) -> MRHamiltonian:
    rng = np.random.default_rng(seed)
    one_body = rng.normal(size=(norb, norb))
    one_body = 0.5 * (one_body + one_body.T)
    two_body = antisymmetrize_two_body(rng.normal(size=(norb,) * 4))
    two_body = 0.5 * (two_body + two_body.transpose(2, 3, 0, 1))
    return MRHamiltonian(-0.4, one_body, two_body)


class GeneratorTests(unittest.TestCase):
    def test_masks_use_ho_quanta_and_preserve_internal_blocks(self) -> None:
        orbits = np.array(
            [
                [0, 0, 0, 1, -1, -1],
                [1, 0, 1, 1, -1, -1],
                [2, 1, 0, 1, -1, -1],
            ]
        )
        quanta = oscillator_quanta_from_orbits(orbits)
        np.testing.assert_array_equal(quanta, [0, 1, 2])
        mask1, mask2 = decoupling_masks(quanta)
        self.assertFalse(mask1[1, 1])
        self.assertTrue(mask1[0, 2])
        self.assertFalse(mask2[0, 2, 1, 1])  # 0+2 == 1+1
        self.assertTrue(mask2[0, 1, 1, 2])

    def test_generator_is_antihermitian_and_masked(self) -> None:
        norb = 4
        densities = compute_densities(
            occupations((0, 1), (2, 3), norb=norb),
            np.array([np.sqrt(0.4), -np.sqrt(0.6)]),
        )
        hamiltonian = random_hermitian_hamiltonian(norb, 531)
        quanta = np.array([0, 0, 1, 2])
        eta = white_generator(hamiltonian, densities, quanta)
        mask1, mask2 = decoupling_masks(quanta)
        np.testing.assert_allclose(eta.one_body, -eta.one_body.T, atol=2e-12)
        np.testing.assert_allclose(
            eta.two_body, -eta.two_body.transpose(2, 3, 0, 1), atol=2e-12
        )
        np.testing.assert_allclose(eta.two_body, -eta.two_body.swapaxes(0, 1), atol=2e-12)
        np.testing.assert_allclose(eta.two_body, -eta.two_body.swapaxes(2, 3), atol=2e-12)
        np.testing.assert_allclose(eta.one_body[~mask1], 0.0, atol=1e-14)
        np.testing.assert_allclose(eta.two_body[~mask2], 0.0, atol=1e-14)

    def test_selected_elements_equal_explicit_decoupling_residuals(self) -> None:
        norb = 4
        densities = compute_densities(
            occupations((0, 1), (2, 3), norb=norb),
            np.array([np.sqrt(0.25), np.sqrt(0.75)]),
        )
        hamiltonian = random_hermitian_hamiltonian(norb, 8128)
        quanta = np.arange(norb)
        residual = masked_decoupling_residual(hamiltonian, densities, quanta)

        for i, j in ((0, 1), (2, 0), (3, 1)):
            basis = np.zeros((norb, norb))
            basis[i, j] = 1.0
            elementary = MRHamiltonian(0.0, basis, np.zeros((norb,) * 4))
            expected = commutator(hamiltonian, elementary, densities).zero_body
            self.assertAlmostEqual(residual.one_body[i, j], expected, places=11)

        for i, j, k, l in ((0, 1, 2, 3), (0, 2, 1, 3), (1, 3, 0, 2)):
            basis = np.zeros((norb,) * 4)
            # An elementary antisymmetric coefficient tensor represents the
            # single normal string :A^{ij}_{kl}: after the 1/4 prefactor.
            basis[i, j, k, l] = 1.0
            basis[j, i, k, l] = -1.0
            basis[i, j, l, k] = -1.0
            basis[j, i, l, k] = 1.0
            elementary = MRHamiltonian(0.0, np.zeros((norb, norb)), basis)
            expected = commutator(hamiltonian, elementary, densities).zero_body
            self.assertAlmostEqual(residual.two_body[i, j, k, l], expected, places=11)

    def test_single_reference_one_body_limit_matches_white_generator(self) -> None:
        densities = compute_densities(occupations((0,), norb=2), np.array([1.0]))
        one_body = np.array([[0.0, 0.2], [0.2, 2.0]])
        hamiltonian = MRHamiltonian(0.0, one_body, np.zeros((2,) * 4))
        eta = white_generator(hamiltonian, densities, np.array([0, 1]))
        self.assertAlmostEqual(eta.one_body[1, 0], 0.1)
        self.assertAlmostEqual(eta.one_body[0, 1], -0.1)
        derivative = commutator(eta, hamiltonian, densities)
        self.assertLess(derivative.one_body[1, 0], 0.0)

    def test_directional_two_body_elements_have_required_antisymmetry(self) -> None:
        norb = 4
        densities = compute_densities(
            occupations((0, 1), (2, 3), norb=norb),
            np.array([np.sqrt(0.4), -np.sqrt(0.6)]),
        )
        hamiltonian = random_hermitian_hamiltonian(norb, 614)
        _, d2 = decoupling_matrix_elements(hamiltonian, densities)
        np.testing.assert_allclose(d2, -d2.swapaxes(0, 1), atol=2e-12)
        np.testing.assert_allclose(d2, -d2.swapaxes(2, 3), atol=2e-12)

    def test_white_ncsm_numerator_omits_irreducible_density_terms(self) -> None:
        norb = 4
        densities = compute_densities(
            occupations((0, 1), (2, 3), norb=norb),
            np.array([np.sqrt(0.4), -np.sqrt(0.6)]),
        )
        hamiltonian = random_hermitian_hamiltonian(norb, 901)
        d1, d2 = white_ncsm_matrix_elements(hamiltonian, densities)
        n = np.diag(densities.gamma1)
        nbar = 1.0 - n
        np.testing.assert_allclose(
            d1,
            n[None, :] * nbar[:, None] * hamiltonian.one_body.T,
            atol=2e-12,
        )
        np.testing.assert_allclose(
            d2,
            nbar[:, None, None, None]
            * nbar[None, :, None, None]
            * n[None, None, :, None]
            * n[None, None, None, :]
            * hamiltonian.two_body.transpose(2, 3, 0, 1),
            atol=2e-12,
        )

        # The published White-NCSM truncation uses gamma1 occupations but no
        # lambda2 terms in the generator numerator.
        altered = type(densities)(
            gamma1=densities.gamma1,
            gamma2=densities.gamma2,
            lambda2=3.0 * densities.lambda2,
        )
        changed1, changed2 = white_ncsm_matrix_elements(hamiltonian, altered)
        np.testing.assert_array_equal(changed1, d1)
        np.testing.assert_array_equal(changed2, d2)
        numerator = white_ncsm_numerator_residual(hamiltonian, densities)
        np.testing.assert_allclose(numerator.one_body, d1 - d1.T, atol=2e-12)
        np.testing.assert_allclose(
            numerator.two_body,
            d2 - d2.transpose(2, 3, 0, 1),
            atol=2e-12,
        )

    def test_epstein_nesbet_denominator_reduces_to_imsrg_sr_formula(self) -> None:
        norb = 6
        occupied = (0, 1)
        particles = (2, 3, 4, 5)
        densities = compute_densities(
            occupations(occupied, norb=norb), np.array([1.0])
        )
        hamiltonian = random_hermitian_hamiltonian(norb, 1701)
        delta1, delta2 = epstein_nesbet_denominators(hamiltonian, densities)
        f = hamiltonian.one_body
        gamma = hamiltonian.two_body
        for particle in particles:
            for hole in occupied:
                expected = (
                    f[particle, particle]
                    - f[hole, hole]
                    - gamma[particle, hole, particle, hole]
                )
                self.assertAlmostEqual(delta1[particle, hole], expected, places=12)
        for particle1, particle2 in ((2, 3), (4, 5)):
            hole1, hole2 = occupied
            expected = (
                f[particle1, particle1]
                + f[particle2, particle2]
                - f[hole1, hole1]
                - f[hole2, hole2]
                + gamma[particle1, particle2, particle1, particle2]
                + gamma[hole1, hole2, hole1, hole2]
                - gamma[particle1, hole1, particle1, hole1]
                - gamma[particle1, hole2, particle1, hole2]
                - gamma[particle2, hole1, particle2, hole1]
                - gamma[particle2, hole2, particle2, hole2]
            )
            self.assertAlmostEqual(
                delta2[particle1, particle2, hole1, hole2], expected, places=12
            )

    def test_correlated_two_body_denominator_is_pair_symmetric(self) -> None:
        norb = 4
        densities = compute_densities(
            occupations((0, 1), (2, 3), norb=norb),
            np.array([np.sqrt(0.3), np.sqrt(0.7)]),
        )
        hamiltonian = random_hermitian_hamiltonian(norb, 701)
        _, delta2 = epstein_nesbet_denominators(hamiltonian, densities)
        np.testing.assert_allclose(delta2, delta2.swapaxes(0, 1), atol=2e-12)
        np.testing.assert_allclose(delta2, delta2.swapaxes(2, 3), atol=2e-12)

    def test_denominator_cutoff_preserves_nonzero_sign(self) -> None:
        values = np.array([-1e-12, 0.0, 1e-12, -2.0, 3.0])
        np.testing.assert_array_equal(
            _safe_denominator(values, 1e-6),
            np.array([-1e-6, 1e-6, 1e-6, -2.0, 3.0]),
        )


if __name__ == "__main__":
    unittest.main()
