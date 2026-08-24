import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from densities import compute_densities
from flow import FlowSettings, _pack, _unpack, integrate_flow
from normal_order import MRHamiltonian


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

    def test_equal_ho_quanta_are_not_decoupled(self) -> None:
        densities = compute_densities(occupations((0,), norb=2), np.array([1.0]))
        initial = MRHamiltonian(
            0.0, np.array([[0.0, 0.2], [0.2, 2.0]]), np.zeros((2,) * 4)
        )
        result = integrate_flow(initial, densities, np.array([0, 0]))
        self.assertTrue(result.converged)
        self.assertEqual(result.function_evaluations, 0)
        np.testing.assert_array_equal(result.hamiltonian.one_body, initial.one_body)


if __name__ == "__main__":
    unittest.main()
