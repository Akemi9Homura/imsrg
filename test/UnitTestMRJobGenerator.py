#!/usr/bin/env python3
"""Regression checks for the single-point C++ J-scheme MR job generator."""

from dataclasses import replace
import json
from pathlib import Path
import sys
import tempfile


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))

from gen_job import MRJschemeSettings, generate_mr_jscheme_slurm  # noqa: E402


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: UnitTestMRJobGenerator.py <imsrg++ executable>")
    executable = Path(sys.argv[1]).resolve()

    with tempfile.TemporaryDirectory(prefix="mr-job-generator-") as temporary:
        root = Path(temporary)
        interaction = root / "bare.jcoupled64"
        reference = root / "reference.jref"
        interaction.write_bytes(b"test interaction\n")
        reference.write_bytes(b"test reference\n")
        settings = MRJschemeSettings(
            nucleus="He4",
            nrefmax=2,
            emax=4,
            interaction=interaction,
            reference_file=reference,
            start_s=40.0,
            target_s=100.0,
            method="flow",
            ode_tolerance=1e-9,
            partition="compute_C",
            cpus=64,
            label="restart",
            result_root=root / "result",
            executable=executable,
        )
        script = generate_mr_jscheme_slurm(settings)
        contents = script.read_text(encoding="utf-8")
        metadata = json.loads((script.parent / "metadata.json").read_text(encoding="utf-8"))

        require(metadata["schema"] == "mrimsrg_cpp_jscheme_slurm_v1",
                "unexpected manifest schema")
        require(metadata["cumulative_start_s"] == 40.0, "lost cumulative start s")
        require(metadata["cumulative_target_s"] == 100.0, "lost cumulative target s")
        require(metadata["segment_smax"] == 60.0, "restart segment length is wrong")
        required_tokens = (
            "#SBATCH --partition=compute_C",
            "#SBATCH --cpus-per-task=64",
            "fmt2=jcoupled64",
            "3bme=none",
            "reference=He4",
            "valence_space=He4",
            "A=4",
            "hw=20",
            "emax=4",
            "basis=oscillator",
            "method=flow",
            "nsteps=1",
            "core_generator=white-ncsm",
            "denominator_partitioning=Epstein_Nesbet",
            "denominator_delta=0",
            "BetaCM=0",
            "nucleon_mass_correction=false",
            "mr_validation_tolerance=1e-10",
            "smax=60",
            "ode_tolerance=1.0000000000000001e-09",
            "write_H_jcoupled64=",
            "write_H_no2bpack=",
            "sha256sum -c -",
            "git diff --quiet HEAD --",
            "ldd ",
        )
        for token in required_tokens:
            require(token in contents, "generated script is missing: " + token)
        require(len(metadata["executable_sha256"]) == 64,
                "executable SHA-256 was not recorded")

        try:
            generate_mr_jscheme_slurm(settings)
        except FileExistsError:
            pass
        else:
            raise AssertionError("generator silently overwrote an existing result")

        invalid = replace(settings, nucleus="Be8", nrefmax=2,
                          result_root=root / "invalid")
        try:
            generate_mr_jscheme_slurm(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError("unsupported Be8 Nrefmax=2 was accepted")

    print("MR J-scheme job-generator regression passed")


if __name__ == "__main__":
    main()
