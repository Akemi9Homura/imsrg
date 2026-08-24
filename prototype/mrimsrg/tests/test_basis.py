import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from basis import prepare_natural_basis, transform_array, transform_hamiltonian
from densities import compute_densities
from normal_order import MRHamiltonian, VacuumHamiltonian, normal_order


def antisymmetrize_two_body(raw: np.ndarray) -> np.ndarray:
    return 0.25 * (
        raw
        - raw.swapaxes(0, 1)
        - raw.swapaxes(2, 3)
        + raw.transpose(1, 0, 3, 2)
    )


class BasisTests(unittest.TestCase):
    def test_sparse_tensor_transform_matches_dense_contraction_and_round_trip(self) -> None:
        angle = 0.37
        vectors = np.eye(4)
        vectors[np.ix_([0, 2], [0, 2])] = np.array(
            [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
        )
        rng = np.random.default_rng(77)
        raw = rng.normal(size=(4,) * 4)
        tensor = antisymmetrize_two_body(raw)
        transformed = transform_array(tensor, vectors, to_natural=True)
        expected = np.einsum(
            "pa,qb,rc,sd,pqrs->abcd",
            vectors,
            vectors,
            vectors,
            vectors,
            tensor,
            optimize=True,
        )
        np.testing.assert_allclose(transformed, expected, atol=2e-13)
        np.testing.assert_allclose(
            transform_array(transformed, vectors, to_natural=False),
            tensor,
            atol=2e-13,
        )

    def test_naturalization_is_covariant_with_mr_normal_ordering(self) -> None:
        determinants = np.array([[1, 0], [0, 1]], dtype=np.uint8)
        coefficients = np.array([3.0, 4.0]) / 5.0
        densities = compute_densities(determinants, coefficients)
        basis = prepare_natural_basis(densities)
        self.assertFalse(basis.is_identity)
        np.testing.assert_allclose(
            basis.densities.gamma1,
            np.diag(np.diag(basis.densities.gamma1)),
            atol=2e-14,
        )

        rng = np.random.default_rng(910)
        one_body = rng.normal(size=(2, 2))
        one_body = 0.5 * (one_body + one_body.T)
        two_body = antisymmetrize_two_body(rng.normal(size=(2,) * 4))
        two_body = 0.5 * (two_body + two_body.transpose(2, 3, 0, 1))
        vacuum = VacuumHamiltonian(0.2, one_body, two_body)
        mr_original = normal_order(vacuum, densities)
        vacuum_natural = VacuumHamiltonian(
            vacuum.zero_body,
            transform_array(vacuum.one_body, basis.vectors, to_natural=True),
            transform_array(vacuum.two_body, basis.vectors, to_natural=True),
        )
        mr_natural = normal_order(vacuum_natural, basis.densities)
        expected = transform_hamiltonian(
            mr_original, basis.vectors, to_natural=True
        )
        self.assertAlmostEqual(mr_natural.zero_body, expected.zero_body, places=13)
        np.testing.assert_allclose(mr_natural.one_body, expected.one_body, atol=2e-13)
        np.testing.assert_allclose(mr_natural.two_body, expected.two_body, atol=2e-13)

    def test_identity_density_does_not_transform(self) -> None:
        densities = compute_densities(
            np.array([[1, 0]], dtype=np.uint8), np.array([1.0])
        )
        basis = prepare_natural_basis(densities)
        self.assertTrue(basis.is_identity)
        np.testing.assert_array_equal(basis.vectors, np.eye(2))


if __name__ == "__main__":
    unittest.main()
