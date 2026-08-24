import hashlib
import json
import sys
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from reference_io import load_reference


class ReferenceIOTests(unittest.TestCase):
    def test_relocated_interaction_is_still_verified_by_digest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            interaction = root / "relocated.minipack"
            interaction.write_bytes(b"fixed interaction")
            digest = hashlib.sha256(interaction.read_bytes()).hexdigest()
            reference = root / "reference"
            reference.mkdir()
            metadata = {
                "schema": "mrimsrg_reference_v1",
                "A": 1,
                "Z": 1,
                "N": 0,
                "hw": 20,
                "emax": 2,
                "e2max": 4,
                "J2": 0,
                "parity": 1,
                "interaction": "/unavailable/original/path",
                "interaction_sha256": digest,
                "one_body_convention": "t[p,q] a^dagger_p a_q",
                "two_body_convention": "(1/4) V[p,q,r,s] a^dagger_p a^dagger_q a_s a_r",
            }
            (reference / "metadata.json").write_text(
                json.dumps(metadata), encoding="utf-8"
            )
            np.save(reference / "orbits.npy", np.array([[0, 0, 0, 1, 1, -1]]))
            np.save(reference / "one_body.npy", np.zeros((1, 1)))
            np.save(reference / "two_body.npy", np.zeros((1, 1, 1, 1)))
            np.save(reference / "determinants.npy", np.ones((1, 1), dtype=np.uint8))
            np.save(reference / "coefficients.npy", np.ones(1))

            with patch("reference_io.FIXED_INTERACTION_SHA256", digest):
                data = load_reference(reference, interaction_path=interaction)
            self.assertEqual(
                data.metadata["validated_interaction_path"], str(interaction.resolve())
            )
            interaction.write_bytes(b"different interaction")
            with patch("reference_io.FIXED_INTERACTION_SHA256", digest):
                with self.assertRaisesRegex(ValueError, "unexpected interaction SHA-256"):
                    load_reference(reference, interaction_path=interaction)


if __name__ == "__main__":
    unittest.main()
