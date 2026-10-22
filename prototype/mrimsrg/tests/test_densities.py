import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from densities import compute_densities, validate_densities


def occupations(*occupied_sets: tuple[int, ...], norb: int = 4) -> np.ndarray:
    result = np.zeros((len(occupied_sets), norb), dtype=np.uint8)
    for row, occupied in enumerate(occupied_sets):
        result[row, list(occupied)] = 1
    return result


class DensityTests(unittest.TestCase):
    def test_single_slater_has_zero_two_body_cumulant(self) -> None:
        determinants = occupations((0, 2))
        densities = compute_densities(determinants, np.array([1.0]))
        validate_densities(densities, 2)
        np.testing.assert_allclose(densities.lambda2, 0.0, atol=1e-14)
        np.testing.assert_allclose(np.diag(densities.gamma1), [1.0, 0.0, 1.0, 0.0])

    def test_correlated_pair_superposition_has_nonzero_cumulant(self) -> None:
        determinants = occupations((0, 1), (2, 3))
        coefficients = np.array([1.0, 1.0]) / np.sqrt(2.0)
        densities = compute_densities(determinants, coefficients)
        validate_densities(densities, 2)
        np.testing.assert_allclose(np.diag(densities.gamma1), 0.5, atol=1e-14)
        self.assertGreater(np.max(np.abs(densities.lambda2)), 0.49)

    def test_coherent_one_particle_superposition_keeps_offdiagonal_density(self) -> None:
        determinants = occupations((0,), (1,))
        coefficients = np.array([3.0, 4.0]) / 5.0
        densities = compute_densities(determinants, coefficients)
        validate_densities(densities, 1)
        expected = np.zeros((4, 4))
        expected[:2, :2] = [[0.36, 0.48], [0.48, 0.64]]
        np.testing.assert_allclose(densities.gamma1, expected, atol=1e-14)
        np.testing.assert_allclose(densities.gamma2, 0.0, atol=1e-14)


if __name__ == "__main__":
    unittest.main()
