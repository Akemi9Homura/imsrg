import sys
from pathlib import Path
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from gen_flow_job import JobSettings, generate_job


class FlowJobTests(unittest.TestCase):
    def test_generates_one_checked_point7_job(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = (
                root
                / "prototype"
                / "mrimsrg"
                / "data"
                / "He4_Nrefmax0_final"
            )
            reference.mkdir(parents=True)
            (reference / "metadata.json").write_text("{}\n", encoding="utf-8")
            interaction = root / "interaction.minipack"
            interaction.write_bytes(b"test")
            script = generate_job(
                root,
                root / "result",
                JobSettings(
                    nucleus="He4", interaction=interaction, nodelist="node2"
                ),
            )
            contents = script.read_text(encoding="utf-8")
            self.assertIn("#SBATCH --partition=c128m512", contents)
            self.assertIn("#SBATCH --qos=low", contents)
            self.assertIn("#SBATCH --cpus-per-task=64", contents)
            self.assertIn("#SBATCH --nodelist=node2", contents)
            self.assertIn("source /opt/modules/init/bash", contents)
            self.assertIn("source ./sourceme.sh", contents)
            self.assertIn("--checkpoint-s 1", contents)
            self.assertIn(f"--interaction {interaction}", contents)
            self.assertIn("--residual-ratio 9.9999999999999995e-07", contents)
            self.assertNotIn("%N", contents)

    def test_rejects_checkpoint_outside_flow_interval(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            interaction = root / "interaction.minipack"
            interaction.write_bytes(b"test")
            with self.assertRaises(ValueError):
                generate_job(
                    root,
                    root / "result",
                    JobSettings(
                        nucleus="He4",
                        interaction=interaction,
                        smax=1.0,
                        checkpoint_s=1.0,
                    ),
                )


if __name__ == "__main__":
    unittest.main()
