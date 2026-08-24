import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from commutator import commutator
from densities import compute_densities
from generator import brillouin_generator, decoupling_masks, oscillator_quanta_from_orbits
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
        eta = brillouin_generator(hamiltonian, densities, quanta)
        mask1, mask2 = decoupling_masks(quanta)
        np.testing.assert_allclose(eta.one_body, -eta.one_body.T, atol=2e-12)
        np.testing.assert_allclose(
            eta.two_body, -eta.two_body.transpose(2, 3, 0, 1), atol=2e-12
        )
        np.testing.assert_allclose(eta.two_body, -eta.two_body.swapaxes(0, 1), atol=2e-12)
        np.testing.assert_allclose(eta.two_body, -eta.two_body.swapaxes(2, 3), atol=2e-12)
        np.testing.assert_allclose(eta.one_body[~mask1], 0.0, atol=1e-14)
        np.testing.assert_allclose(eta.two_body[~mask2], 0.0, atol=1e-14)

    def test_selected_elements_equal_explicit_brillouin_residuals(self) -> None:
        norb = 4
        densities = compute_densities(
            occupations((0, 1), (2, 3), norb=norb),
            np.array([np.sqrt(0.25), np.sqrt(0.75)]),
        )
        hamiltonian = random_hermitian_hamiltonian(norb, 8128)
        quanta = np.arange(norb)
        eta = brillouin_generator(hamiltonian, densities, quanta)

        for i, j in ((0, 1), (2, 0), (3, 1)):
            basis = np.zeros((norb, norb))
            basis[i, j] = 1.0
            elementary = MRHamiltonian(0.0, basis, np.zeros((norb,) * 4))
            expected = commutator(hamiltonian, elementary, densities).zero_body
            self.assertAlmostEqual(eta.one_body[i, j], expected, places=11)

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
            self.assertAlmostEqual(eta.two_body[i, j, k, l], expected, places=11)

    def test_single_reference_one_body_limit_matches_imaginary_time_sign(self) -> None:
        densities = compute_densities(occupations((0,), norb=2), np.array([1.0]))
        one_body = np.array([[0.0, 0.2], [0.2, 2.0]])
        hamiltonian = MRHamiltonian(0.0, one_body, np.zeros((2,) * 4))
        eta = brillouin_generator(hamiltonian, densities, np.array([0, 1]))
        self.assertAlmostEqual(eta.one_body[1, 0], 0.2)
        self.assertAlmostEqual(eta.one_body[0, 1], -0.2)
        derivative = commutator(eta, hamiltonian, densities)
        self.assertLess(derivative.one_body[1, 0], 0.0)


if __name__ == "__main__":
    unittest.main()
