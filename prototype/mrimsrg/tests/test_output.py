import json
import sys
from pathlib import Path
import tempfile
import unittest
import struct

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from densities import compute_densities
from flow import FlowPoint, FlowResult, FlowSettings
from normal_order import VacuumHamiltonian, normal_order
from output import save_flow_output
from reference_io import ReferenceData


class OutputTests(unittest.TestCase):
    def test_saves_explicit_vacuum_representation_and_refuses_overwrite(self) -> None:
        orbits = np.array([[0, 0, 0, 1, -1, -1], [1, 0, 1, 1, 1, 1]])
        determinants = np.array([[1, 0]], dtype=np.uint8)
        coefficients = np.array([1.0])
        densities = compute_densities(determinants, coefficients)
        vacuum = VacuumHamiltonian(0.3, np.diag([1.0, 2.0]), np.zeros((2,) * 4))
        initial = normal_order(vacuum, densities)
        point = FlowPoint(0, 0.0, initial.zero_body, 0.0, 0.0, 0.0, 0.0, 0.0)
        result = FlowResult(initial, (point,), 0, True, "test")
        reference = ReferenceData(
            metadata={"A": 1, "interaction_sha256": "test"},
            orbits=orbits,
            one_body=vacuum.one_body,
            two_body=vacuum.two_body,
            determinants=determinants,
            coefficients=coefficients,
        )
        with tempfile.TemporaryDirectory() as temporary:
            output_path = Path(temporary) / "result"
            recovered = save_flow_output(
                output_path, "reference", reference, densities, initial, result, FlowSettings()
            )
            np.testing.assert_allclose(np.load(output_path / "final_vacuum_one_body.npy"), vacuum.one_body)
            self.assertAlmostEqual(recovered.zero_body, vacuum.zero_body)
            metadata = json.loads((output_path / "metadata.json").read_text())
            self.assertEqual(metadata["density_approximation"], "lambda3=0")
            self.assertAlmostEqual(metadata["final_vacuum_zero_body"], vacuum.zero_body)
            with (output_path / "vacuum_mscheme.bin").open("rb") as stream:
                self.assertEqual(stream.read(16), b"mrimsrg_m_v1\0\0\0\0")
                norb, zero_body = struct.unpack("<Qd", stream.read(16))
                self.assertEqual(norb, 2)
                self.assertAlmostEqual(zero_body, vacuum.zero_body)
            with self.assertRaises(FileExistsError):
                save_flow_output(
                    output_path, "reference", reference, densities, initial, result, FlowSettings()
                )


if __name__ == "__main__":
    unittest.main()
