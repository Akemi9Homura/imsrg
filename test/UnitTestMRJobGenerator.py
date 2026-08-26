#!/usr/bin/env python3
"""Regression checks for the single-point C++ J-scheme MR job generator."""

from dataclasses import replace
import json
from pathlib import Path
import subprocess
import sys
import tempfile


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))

from gen_job import (  # noqa: E402
    MRJschemeInputSettings,
    MRJschemeSettings,
    MRNCSMReadbackSettings,
    generate_mr_jscheme_input_slurm,
    generate_mr_jscheme_slurm,
    generate_mr_ncsm_readback_slurm,
)


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
            partition="compute_C",
            cpus=64,
            label="restart",
            result_root=root / "result",
            executable=executable,
        )
        script = generate_mr_jscheme_slurm(settings)
        subprocess.run(["bash", "-n", str(script)], check=True)
        contents = script.read_text(encoding="utf-8")
        metadata = json.loads((script.parent / "metadata.json").read_text(encoding="utf-8"))

        require(metadata["schema"] == "mrimsrg_cpp_jscheme_slurm_v3",
                "unexpected manifest schema")
        require(metadata["cumulative_start_s"] == 40.0, "lost cumulative start s")
        require(metadata["cumulative_target_s"] == 100.0, "lost cumulative target s")
        require(metadata["segment_smax"] == 60.0, "restart segment length is wrong")
        require(metadata["solver_method"] == "magnus",
                "MR production solver is not fixed to the existing Magnus path")
        require(metadata["ode_parameter_source"] == "imsrg++ Magnus runtime defaults",
                "MR Magnus does not record the source of its ODE parameters")
        require(metadata["ode_parameter_overrides"] == [],
                "MR flow still records wrapper-level ODE overrides")
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
            "method=magnus",
            "nsteps=1",
            "core_generator=white-ncsm",
            "denominator_partitioning=Epstein_Nesbet",
            "denominator_delta=0",
            "BetaCM=0",
            "nucleon_mass_correction=false",
            "mr_validation_tolerance=1e-10",
            "smax=60",
            "write_H_jcoupled64=",
            "write_H_no2bpack=",
            "write_omega=true",
            "scratch=",
            "compgen -G",
            "echo repository_commit=",
            "/proc/$flow_pid/status",
            "echo wall_seconds=",
            "echo maximum_rss_kib=",
            "resource_usage.txt",
            "sha256sum -c -",
            "ldd ",
        )
        for token in required_tokens:
            require(token in contents, "generated script is missing: " + token)
        for token in ("ds_0=", "dsmax=", "ode_tolerance="):
            require(token not in contents,
                    "generated MR script overrides an imsrg++ ODE default: " + token)
        require(len(metadata["executable_sha256"]) == 64,
                "executable SHA-256 was not recorded")
        require(len(metadata["shared_library_sha256"]) == 64,
                "libIMSRG SHA-256 was not recorded")
        require(len(metadata["environment_script_sha256"]) == 64,
                "sourceme.sh SHA-256 was not recorded")
        require(metadata["repository_state_files"],
                "repository state files were not recorded")
        require("git rev-parse" not in contents,
                "flow job still requires Git on the compute node")
        for state_file in metadata["repository_state_files"]:
            require(len(state_file["sha256"]) == 64,
                    "repository-state SHA-256 was not recorded")
            require(state_file["path"] in contents,
                    "flow job does not verify a repository-state file")
        require(metadata["shared_library"] in contents,
                "generated script does not verify libIMSRG")
        require(metadata["environment_script"] in contents,
                "generated script does not verify sourceme.sh")
        require(metadata["resource_usage"] in contents,
                "generated script does not preserve sampled flow resources")
        require(metadata["omega_represents"] == "segment",
                "Magnus metadata does not identify a segment Omega")
        require(metadata["omega_segment_start_s"] == 40.0,
                "Magnus metadata lost its segment start")
        require(metadata["omega_segment_target_s"] == 100.0,
                "Magnus metadata lost its segment target")
        require(metadata["omega_prefix"] in contents,
                "generated script does not verify written Omega files")

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

        pyimsrg_dir = root / "pyimsrg"
        pyimsrg_dir.mkdir()
        (pyimsrg_dir / "pyIMSRG.so").write_bytes(b"test pyimsrg module\n")
        input_settings = MRJschemeInputSettings(
            nucleus="He4",
            nrefmax=2,
            emax=12,
            minipack=interaction,
            source_reference=reference,
            partition="c128m1024",
            cpus=128,
            label="capacity",
            result_root=root / "input-result",
            converter=executable,
            pyimsrg_dir=pyimsrg_dir,
        )
        input_script = generate_mr_jscheme_input_slurm(input_settings)
        subprocess.run(["bash", "-n", str(input_script)], check=True)
        input_contents = input_script.read_text(encoding="utf-8")
        input_metadata = json.loads(
            (input_script.parent / "metadata.json").read_text(encoding="utf-8")
        )
        require(input_metadata["schema"] == "mrimsrg_jscheme_input_slurm_v1",
                "unexpected MR input manifest schema")
        for token in (
            "#SBATCH --partition=c128m1024",
            "#SBATCH --cpus-per-task=128",
            "echo repository_commit=",
            "run_with_rss()",
            "/proc/$sampled_pid/status",
            "echo wall_seconds=",
            "echo maximum_rss_kib=",
            "run_with_rss converter",
            "run_with_rss embed_reference",
            "--interaction",
            "--output",
            "--A 4",
            "embed_jref.py",
            "jref_format.py",
            "pyimsrg_utils.py",
            "--source",
            "--pyimsrg-dir",
            "converter.time",
            "embed.time",
            "sha256sum -c -",
            "ldd ",
        ):
            require(token in input_contents,
                    "generated MR input script is missing: " + token)
        require("/usr/bin/time" not in input_contents,
                "point7 MR input script still requires unavailable /usr/bin/time")
        require("git rev-parse" not in input_contents,
                "input job still requires Git on the compute node")
        require(input_metadata["repository_state_files"],
                "input job lost repository-state hashes")
        for key in (
            "minipack_sha256", "source_reference_sha256", "converter_sha256",
            "pyimsrg_sha256", "embedder_sha256", "format_module_sha256",
            "pyimsrg_loader_sha256", "environment_script_sha256",
        ):
            require(len(input_metadata[key]) == 64,
                    "MR input metadata lost SHA-256: " + key)
        for key in ("format_module", "pyimsrg_loader"):
            require(input_metadata[key] in input_contents,
                    "MR input script does not verify dependency: " + key)
        blocked_sympy_import = r'''\
import builtins
original_import = builtins.__import__
def blocked_import(name, *args, **kwargs):
    if name == "sympy" or name.startswith("sympy."):
        raise ModuleNotFoundError("blocked by dependency-isolation regression")
    return original_import(name, *args, **kwargs)
builtins.__import__ = blocked_import
import prototype.mrimsrg.embed_jref
'''
        subprocess.run(
            [sys.executable, "-c", blocked_sympy_import],
            cwd=REPOSITORY,
            check=True,
        )
        try:
            generate_mr_jscheme_input_slurm(input_settings)
        except FileExistsError:
            pass
        else:
            raise AssertionError("MR input generator silently overwrote a result")

        packed = root / "flow.no2bpack"
        packed.write_bytes(b"packed Hamiltonian\n")
        ncsm_settings = MRNCSMReadbackSettings(
            no2bpack=packed,
            proton_number=2,
            neutron_number=2,
            nmax=10,
            states=3,
            max_iter=400,
            partition="compute_C",
            cpus=64,
            label="flow_s0p02",
            result_root=root / "ncsm-result",
            executable=executable,
        )
        ncsm_script = generate_mr_ncsm_readback_slurm(ncsm_settings)
        ncsm_contents = ncsm_script.read_text(encoding="utf-8")
        ncsm_metadata = json.loads(
            (ncsm_script.parent / "metadata.json").read_text(encoding="utf-8")
        )
        require(ncsm_metadata["schema"] == "mrimsrg_ncsm_readback_slurm_v1",
                "unexpected NCSM manifest schema")
        for token in (
            "#SBATCH --partition=compute_C",
            "#SBATCH --cpus-per-task=64",
            "echo repository_commit=",
            "start_seconds=$SECONDS",
            "maximum_rss_kib=0",
            "/proc/$ncsm_pid/status",
            "echo wall_seconds=",
            "echo maximum_rss_kib=",
            "--no2bpack",
            "--Z 2",
            "--N 2",
            "--nmax 10",
            "--states 3",
            "--max-iter 400",
            "sha256sum -c -",
            "ldd ",
        ):
            require(token in ncsm_contents,
                    "generated NCSM script is missing: " + token)
        require(len(ncsm_metadata["no2bpack_sha256"]) == 64,
                "NCSM input SHA-256 was not recorded")
        require(len(ncsm_metadata["executable_sha256"]) == 64,
                "NCSM executable SHA-256 was not recorded")
        require("git rev-parse" not in ncsm_contents,
                "NCSM job still requires Git on the compute node")
        require(ncsm_metadata["repository_state_files"],
                "NCSM job lost repository-state hashes")

        foreground_settings = replace(
            ncsm_settings,
            label="flow_s0p02_foreground",
            sample_rss=False,
            result_root=root / "ncsm-foreground-result",
        )
        foreground_script = generate_mr_ncsm_readback_slurm(foreground_settings)
        foreground_contents = foreground_script.read_text(encoding="utf-8")
        foreground_metadata = json.loads(
            (foreground_script.parent / "metadata.json").read_text(encoding="utf-8")
        )
        require("/proc/$ncsm_pid/status" not in foreground_contents,
                "foreground NCSM script unexpectedly samples /proc")
        require("maximum_rss_kib=not_sampled" in foreground_contents,
                "foreground NCSM script does not label RSS as unsampled")
        require("--nmax 10" in foreground_contents,
                "foreground NCSM script lost the requested calculation")
        require(foreground_metadata["settings"]["sample_rss"] is False,
                "foreground NCSM mode is not recorded in metadata")

        try:
            generate_mr_ncsm_readback_slurm(ncsm_settings)
        except FileExistsError:
            pass
        else:
            raise AssertionError("NCSM generator silently overwrote an existing result")

    print("MR J-scheme job-generator regression passed")


if __name__ == "__main__":
    main()
