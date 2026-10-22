from __future__ import annotations

import sys
from pathlib import Path
import unittest

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from densities import compute_densities
from magnus import add, bch_product, bch_transform, max_abs, scale
from normal_order import MRHamiltonian
from test_commutator import occupations, random_operator, vacuum_matrix


def assert_adjoint(test: unittest.TestCase, operator: MRHamiltonian, sign: int) -> None:
    test.assertLess(abs(operator.zero_body - sign * operator.zero_body), 2e-11)
    np.testing.assert_allclose(operator.one_body, sign * operator.one_body.T, atol=2e-11)
    np.testing.assert_allclose(
        operator.two_body,
        sign * operator.two_body.transpose(2, 3, 0, 1),
        atol=2e-11,
    )
    np.testing.assert_allclose(
        operator.two_body, -operator.two_body.swapaxes(0, 1), atol=2e-11
    )
    np.testing.assert_allclose(
        operator.two_body, -operator.two_body.swapaxes(2, 3), atol=2e-11
    )


class MagnusTests(unittest.TestCase):
    def setUp(self) -> None:
        self.norb = 4
        self.densities = compute_densities(
            occupations((0, 1), (2, 3), norb=self.norb),
            np.array([np.sqrt(0.4), -np.sqrt(0.6)]),
        )

    def test_bch_terms_are_recursive_factorial_scaled_mr_commutators(self) -> None:
        omega = scale(random_operator(self.norb, 771, -1), 0.025)
        hamiltonian = random_operator(self.norb, 991, 1)
        result = bch_transform(hamiltonian, omega, self.densities, max_order=6)
        self.assertEqual(len(result.terms), 7)
        reconstructed = result.terms[0]
        for order, term in enumerate(result.terms):
            assert_adjoint(self, term, 1)
            if order:
                reconstructed = add(reconstructed, term)
        self.assertLess(max_abs(add(result.operator, scale(reconstructed, -1.0))), 3e-15)
        assert_adjoint(self, result.operator, 1)

    def test_bch_product_orientation_and_first_step(self) -> None:
        omega = scale(random_operator(self.norb, 1234, -1), 0.03)
        d_omega = scale(random_operator(self.norb, 5678, -1), 0.002)
        product = bch_product(d_omega, omega, self.densities, threshold=0.0)
        self.assertEqual(product.contributions[0][0], "domega")
        self.assertEqual(product.contributions[1][0], "omega")
        self.assertIn("bernoulli_1", dict(product.contributions))
        assert_adjoint(self, product.omega, -1)

        zero = scale(omega, 0.0)
        first = bch_product(d_omega, zero, self.densities, threshold=0.0)
        self.assertLess(max_abs(add(first.omega, scale(d_omega, -1.0))), 1e-15)

    def test_one_body_bch_matches_explicit_fock_space_similarity_transform(self) -> None:
        # The commutator of number-conserving one-body operators closes exactly,
        # so this checks the exponent/sign convention independently of the BCH
        # recursion and of MR normal-ordering constants.
        omega = random_operator(self.norb, 808, -1)
        omega = MRHamiltonian(0.0, 0.08 * omega.one_body, np.zeros_like(omega.two_body))
        hamiltonian = random_operator(self.norb, 909, 1)
        hamiltonian = MRHamiltonian(
            0.31, hamiltonian.one_body, np.zeros_like(hamiltonian.two_body)
        )
        transformed = bch_transform(
            hamiltonian, omega, self.densities, max_order=30
        ).operator

        omega_matrix, basis = vacuum_matrix(omega, self.densities, 2)
        initial_matrix, _ = vacuum_matrix(hamiltonian, self.densities, 2)
        actual_matrix, _ = vacuum_matrix(transformed, self.densities, 2)
        expected = expm(omega_matrix) @ initial_matrix @ expm(-omega_matrix)
        self.assertEqual(len(basis), 6)
        np.testing.assert_allclose(actual_matrix, expected, atol=3e-13)

    def test_one_body_bch_product_has_production_exponential_order(self) -> None:
        omega = random_operator(self.norb, 812, -1)
        d_omega = random_operator(self.norb, 813, -1)
        zero2 = np.zeros_like(omega.two_body)
        omega = MRHamiltonian(0.0, 0.02 * omega.one_body, zero2)
        d_omega = MRHamiltonian(0.0, 0.001 * d_omega.one_body, zero2)
        product = bch_product(d_omega, omega, self.densities, threshold=0.0).omega
        product_matrix, _ = vacuum_matrix(product, self.densities, 2)
        omega_matrix, _ = vacuum_matrix(omega, self.densities, 2)
        d_omega_matrix, _ = vacuum_matrix(d_omega, self.densities, 2)
        np.testing.assert_allclose(
            expm(product_matrix),
            expm(d_omega_matrix) @ expm(omega_matrix),
            atol=2e-13,
        )


if __name__ == "__main__":
    unittest.main()
