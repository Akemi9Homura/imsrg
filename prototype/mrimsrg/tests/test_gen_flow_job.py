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
                    nucleus="He4",
                    interaction=interaction,
                    label="smax5",
                    nodelist="node2",
                ),
            )
            self.assertIn("smax5", script.name)
            contents = script.read_text(encoding="utf-8")
            self.assertIn("#SBATCH --partition=c128m512", contents)
            self.assertIn("#SBATCH --qos=low", contents)
            self.assertIn("#SBATCH --cpus-per-task=64", contents)
            self.assertIn("#SBATCH --nodelist=node2", contents)
            self.assertIn("source /opt/modules/init/bash", contents)
            self.assertIn("source ./sourceme.sh", contents)
            self.assertIn("--smax 25000", contents)
            self.assertIn("--checkpoint-s 40", contents)
            self.assertIn("--max-step 10", contents)
            self.assertIn(f"--interaction {interaction}", contents)
            self.assertIn("/opt/library/miniconda-3.12.9/bin/python3 -u", contents)
            self.assertIn("import sys, numpy, scipy", contents)
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

    def test_generates_single_continuation_job(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = (
                root
                / "prototype"
                / "mrimsrg"
                / "data"
                / "He4_Nrefmax2"
            )
            reference.mkdir(parents=True)
            (reference / "metadata.json").write_text("{}\n", encoding="utf-8")
            interaction = root / "interaction.minipack"
            interaction.write_bytes(b"test")
            resume = root / "prior" / "flow"
            resume.mkdir(parents=True)
            (resume / "metadata.json").write_text(
                '{"trajectory": [{"s": 25000.0}]}\n', encoding="utf-8"
            )
            script = generate_job(
                root,
                root / "result",
                JobSettings(
                    nucleus="He4",
                    nrefmax=2,
                    interaction=interaction,
                    label="continue50k",
                    checkpoint_s=40000.0,
                    smax=50000.0,
                    resume_from=resume,
                ),
            )
            contents = script.read_text(encoding="utf-8")
            self.assertIn("--smax 50000", contents)
            self.assertIn("--checkpoint-s 40000", contents)
            self.assertIn(f"--resume-from {resume}", contents)

    def test_rejects_continuation_checkpoint_before_start(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            interaction = root / "interaction.minipack"
            interaction.write_bytes(b"test")
            resume = root / "prior" / "flow"
            resume.mkdir(parents=True)
            (resume / "metadata.json").write_text(
                '{"trajectory": [{"s": 25000.0}]}\n', encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "between start_s"):
                generate_job(
                    root,
                    root / "result",
                    JobSettings(
                        nucleus="He4",
                        interaction=interaction,
                        smax=50000.0,
                        checkpoint_s=40.0,
                        resume_from=resume,
                    ),
                )

    def test_rejects_unplanned_open_shell_nrefmax2_reference(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            interaction = root / "interaction.minipack"
            interaction.write_bytes(b"test")
            with self.assertRaisesRegex(ValueError, "unsupported reference"):
                generate_job(
                    root,
                    root / "result",
                    JobSettings(
                        nucleus="Be8",
                        nrefmax=2,
                        interaction=interaction,
                    ),
                )


if __name__ == "__main__":
    unittest.main()
