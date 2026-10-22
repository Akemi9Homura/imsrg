import sys
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from densities import compute_densities
from basis import prepare_natural_basis, transform_hamiltonian
from flow import FlowSettings, _pack, _unpack, integrate_flow
from generator import (
    decoupling_masks,
    masked_residual_norm,
    white_ncsm_numerator_residual,
)
from normal_order import MRHamiltonian, VacuumHamiltonian, normal_order


def occupations(*occupied_sets: tuple[int, ...], norb: int) -> np.ndarray:
    result = np.zeros((len(occupied_sets), norb), dtype=np.uint8)
    for row, occupied in enumerate(occupied_sets):
        result[row, list(occupied)] = 1
    return result


class FlowTests(unittest.TestCase):
    def test_pack_round_trip(self) -> None:
        rng = np.random.default_rng(12)
        one_body = rng.normal(size=(3, 3))
        one_body = 0.5 * (one_body + one_body.T)
        raw = rng.normal(size=(3, 3, 3, 3))
        two_body = 0.25 * (
            raw
            - raw.swapaxes(0, 1)
            - raw.swapaxes(2, 3)
            + raw.transpose(1, 0, 3, 2)
        )
        two_body = 0.5 * (two_body + two_body.transpose(2, 3, 0, 1))
        hamiltonian = MRHamiltonian(
            0.7, one_body, two_body
        )
        recovered = _unpack(_pack(hamiltonian), 3)
        self.assertEqual(recovered.zero_body, hamiltonian.zero_body)
        np.testing.assert_array_equal(recovered.one_body, hamiltonian.one_body)
        np.testing.assert_allclose(recovered.two_body, hamiltonian.two_body, atol=2e-16)

    def test_two_level_direct_flow_decouples_and_preserves_spectrum(self) -> None:
        densities = compute_densities(occupations((0,), norb=2), np.array([1.0]))
        one_body = np.array([[0.0, 0.2], [0.2, 2.0]])
        initial = MRHamiltonian(0.0, one_body, np.zeros((2,) * 4))
        settings = FlowSettings(
            smax=25.0,
            relative_tolerance=1e-10,
            absolute_tolerance=1e-12,
            max_step=0.5,
            residual_ratio=1e-8,
            checkpoint_s=0.5,
        )
        result = integrate_flow(initial, densities, np.array([0, 1]), settings)
        self.assertTrue(result.converged, result.message)
        self.assertLessEqual(result.trajectory[-1].residual_ratio, 1e-8)
        np.testing.assert_allclose(
            np.linalg.eigvalsh(result.hamiltonian.one_body),
            np.linalg.eigvalsh(initial.one_body),
            atol=2e-10,
        )
        self.assertLess(abs(result.hamiltonian.one_body[0, 1]), 2e-9)
        self.assertEqual(len(result.checkpoints), 1)
        self.assertEqual(result.checkpoints[0].point.s, 0.5)
        self.assertLess(
            abs(result.checkpoints[0].hamiltonian.one_body[0, 1]),
            abs(initial.one_body[0, 1]),
        )

    def test_continuation_preserves_original_residual_normalization(self) -> None:
        densities = compute_densities(occupations((0,), norb=2), np.array([1.0]))
        initial = MRHamiltonian(
            0.0, np.array([[0.0, 0.2], [0.2, 2.0]]), np.zeros((2,) * 4)
        )
        first = integrate_flow(
            initial,
            densities,
            np.array([0, 1]),
            FlowSettings(
                smax=0.4,
                relative_tolerance=1e-10,
                absolute_tolerance=1e-12,
                max_step=0.1,
                residual_ratio=1e-12,
            ),
        )
        initial_point = first.trajectory[0]
        normalization = (
            initial_point.residual,
            initial_point.generator_residual,
            initial_point.generator_numerator_residual,
        )
        continued = integrate_flow(
            first.hamiltonian,
            densities,
            np.array([0, 1]),
            FlowSettings(
                smax=0.8,
                relative_tolerance=1e-10,
                absolute_tolerance=1e-12,
                max_step=0.1,
                residual_ratio=1e-12,
            ),
            start_s=first.trajectory[-1].s,
            residual_normalization=normalization,
        )
        direct = integrate_flow(
            initial,
            densities,
            np.array([0, 1]),
            FlowSettings(
                smax=0.8,
                relative_tolerance=1e-10,
                absolute_tolerance=1e-12,
                max_step=0.1,
                residual_ratio=1e-12,
            ),
        )
        self.assertEqual(continued.trajectory[0].s, 0.4)
        self.assertAlmostEqual(
            continued.trajectory[0].generator_residual_ratio,
            first.trajectory[-1].generator_residual_ratio,
            places=13,
        )
        np.testing.assert_allclose(
            continued.hamiltonian.one_body,
            direct.hamiltonian.one_body,
            atol=2e-11,
        )

    def test_continuation_checkpoint_must_follow_start(self) -> None:
        densities = compute_densities(occupations((0,), norb=2), np.array([1.0]))
        initial = MRHamiltonian(
            0.0, np.diag([0.0, 2.0]), np.zeros((2,) * 4)
        )
        with self.assertRaisesRegex(ValueError, "strictly between start_s"):
            integrate_flow(
                initial,
                densities,
                np.array([0, 1]),
                FlowSettings(smax=2.0, checkpoint_s=0.5),
                start_s=1.0,
            )

    def test_equal_ho_quanta_are_not_decoupled(self) -> None:
        densities = compute_densities(occupations((0,), norb=2), np.array([1.0]))
        initial = MRHamiltonian(
            0.0, np.array([[0.0, 0.2], [0.2, 2.0]]), np.zeros((2,) * 4)
        )
        result = integrate_flow(initial, densities, np.array([0, 0]))
        self.assertTrue(result.converged)
        self.assertEqual(result.function_evaluations, 0)
        np.testing.assert_array_equal(result.hamiltonian.one_body, initial.one_body)

    def test_published_generator_fixed_point_keeps_strict_diagnostic_visible(self) -> None:
        densities = compute_densities(
            occupations((0, 1), (2, 3), norb=4),
            np.array([np.sqrt(0.4), np.sqrt(0.6)]),
        )
        # A diagonal 1B Hamiltonian has a zero lambda-free White-NCSM
        # numerator, while its f-lambda2 terms give a nonzero strict D2.
        initial = MRHamiltonian(
            0.0,
            np.diag([0.0, 1.0, 3.0, 7.0]),
            np.zeros((4,) * 4),
        )
        result = integrate_flow(
            initial,
            densities,
            np.arange(4),
            FlowSettings(smax=1.0, max_step=1.0),
        )
        self.assertTrue(result.converged)
        self.assertAlmostEqual(result.trajectory[-1].residual_ratio, 1.0)
        self.assertEqual(
            result.trajectory[-1].generator_numerator_residual_ratio, 0.0
        )
        np.testing.assert_array_equal(result.hamiltonian.one_body, initial.one_body)
        np.testing.assert_array_equal(result.hamiltonian.two_body, initial.two_body)

    def test_acceptance_uses_actual_generator_not_unweighted_numerator(self) -> None:
        densities = compute_densities(occupations((0,), norb=2), np.array([1.0]))
        initial = MRHamiltonian(
            0.0, np.array([[0.0, 0.2], [0.2, 2.0]]), np.zeros((2,) * 4)
        )
        zero_generator = MRHamiltonian(
            0.0, np.zeros((2, 2)), np.zeros((2,) * 4)
        )
        with patch("flow.white_generator_unmasked", return_value=zero_generator):
            result = integrate_flow(initial, densities, np.array([0, 1]))
        self.assertTrue(result.converged)
        self.assertEqual(result.function_evaluations, 0)
        self.assertEqual(result.trajectory[0].generator_residual_ratio, 0.0)
        self.assertEqual(
            result.trajectory[0].generator_numerator_residual_ratio, 1.0
        )

    def test_nondiagonal_gamma_is_evaluated_and_masked_in_natural_basis(self) -> None:
        determinants = occupations((0,), (1,), norb=2)
        densities = compute_densities(
            determinants, np.array([3.0, 4.0]) / 5.0
        )
        self.assertGreater(abs(densities.gamma1[0, 1]), 0.1)
        vacuum = VacuumHamiltonian(
            0.0, np.diag([0.0, 2.0]), np.zeros((2,) * 4)
        )
        initial = normal_order(vacuum, densities)

        # Equal HO quanta must suppress the flow even though the temporary
        # natural orbitals are mixtures of the two original orbitals.
        blocked = integrate_flow(initial, densities, np.array([0, 0]))
        self.assertTrue(blocked.converged)
        self.assertEqual(blocked.function_evaluations, 0)
        np.testing.assert_allclose(blocked.hamiltonian.one_body, initial.one_body)

        active = integrate_flow(
            initial,
            densities,
            np.array([0, 1]),
            FlowSettings(
                smax=0.2,
                max_step=0.1,
                residual_ratio=1e-12,
            ),
        )
        self.assertGreater(active.function_evaluations, 0)
        np.testing.assert_allclose(
            active.hamiltonian.one_body,
            active.hamiltonian.one_body.T,
            atol=2e-12,
        )

    def test_initial_residual_uses_published_natural_basis_mask(self) -> None:
        determinants = occupations((0,), (1,), norb=3)
        densities = compute_densities(
            determinants, np.array([3.0, 4.0]) / 5.0
        )
        vacuum = VacuumHamiltonian(
            0.0,
            np.array(
                [
                    [0.0, 0.3, 0.7],
                    [0.3, 2.0, -0.4],
                    [0.7, -0.4, 5.0],
                ]
            ),
            np.zeros((3,) * 4),
        )
        initial = normal_order(vacuum, densities)
        quanta = np.array([0, 1, 1])
        natural = prepare_natural_basis(densities)
        working = transform_hamiltonian(
            initial, natural.vectors, to_natural=True
        )
        raw = white_ncsm_numerator_residual(working, natural.densities)
        mask1, mask2 = decoupling_masks(quanta)
        expected = masked_residual_norm(
            MRHamiltonian(0.0, raw.one_body * mask1, raw.two_body * mask2)
        )

        # The obsolete HO-basis projection is intentionally distinguishable
        # for this three-orbit example.
        raw_ho = transform_hamiltonian(raw, natural.vectors, to_natural=False)
        projected_ho = masked_residual_norm(
            MRHamiltonian(
                0.0, raw_ho.one_body * mask1, raw_ho.two_body * mask2
            )
        )
        self.assertGreater(abs(expected - projected_ho), 1e-4)

        result = integrate_flow(
            initial,
            densities,
            quanta,
            FlowSettings(smax=1e-4, max_step=1e-4, residual_ratio=1e-12),
        )
        self.assertAlmostEqual(
            result.trajectory[0].generator_numerator_residual,
            expected,
            places=13,
        )


if __name__ == "__main__":
    unittest.main()
