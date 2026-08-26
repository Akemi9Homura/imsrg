#!/usr/bin/env python3

"""Real He4/O16 exact single-reference degeneration through both drivers."""

import os
from pathlib import Path
import json
import subprocess
import sys
import tempfile

import pyIMSRG


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))
sys.path.insert(0, str(REPOSITORY / "test"))

from UnitTestMRDriver import payload_error, vacuum_operator  # noqa: E402
from prototype.mrimsrg.export_jref import export_reference  # noqa: E402
from prototype.mrimsrg.reference_io import load_reference  # noqa: E402


def run_driver(
    executable,
    root,
    nucleus,
    A,
    input_file,
    reference_file,
    smax,
    mr,
    step=1e-4,
    method="flow_RK4",
):
    label = (
        f"{nucleus}_{'mr' if mr else 'sr'}_{method}_"
        f"{smax:.8g}_ds{step:.8g}"
    )
    output = root / f"{label}.jcoupled64"
    int_prefix = root / ("int_" + label)
    command = [
        str(executable),
        f"2bme={input_file}",
        "fmt2=jcoupled64",
        "3bme=none",
        f"reference={nucleus}",
        f"valence_space={nucleus}",
        f"A={A}",
        "hw=20",
        "emax=2",
        "basis=oscillator",
        f"method={method}",
        f"core_generator={'white-ncsm' if mr else 'white'}",
        f"smax={smax:.8g}",
        f"ds_0={step:.8g}",
        f"dsmax={step:.8g}",
        "eta_criterion=0",
        "BetaCM=0",
        f"flowfile={root / ('flow_' + label + '.txt')}",
        f"intfile={int_prefix}",
        f"write_H_jcoupled64={output}",
    ]
    if method.startswith("magnus"):
        scratch = root / ("scratch_" + label)
        scratch.mkdir()
        command.extend((f"scratch={scratch}", "write_omega=true"))
    if mr:
        command.append(f"mr_reference_file={reference_file}")
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    completed = subprocess.run(
        command,
        text=True,
        capture_output=True,
        env=environment,
        check=False,
    )
    if completed.returncode != 0:
        print(completed.stdout)
        print(completed.stderr, file=sys.stderr)
        raise RuntimeError(f"{label} failed")
    if method.startswith("magnus"):
        omega_files = tuple(root.glob(int_prefix.name + "_Omega_*"))
        if not omega_files or any(path.stat().st_size == 0 for path in omega_files):
            raise RuntimeError(f"{label} did not materialize Omega")
        omega_validator = executable.parent / "imsrg_operator_validate"
        completed = subprocess.run(
            [
                str(omega_validator),
                "--emax", "2",
                "--hw", "20",
                "--reference", nucleus,
                *(str(path) for path in sorted(omega_files)),
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        if completed.returncode != 0:
            print(completed.stdout)
            print(completed.stderr, file=sys.stderr)
            raise RuntimeError(f"{label} wrote an invalid Omega payload")
        report = json.loads(completed.stdout.splitlines()[-1])
        if not report["passed"] or report["files"] != len(omega_files):
            raise RuntimeError(f"{label} Omega validation report is inconsistent")
        if mr and nucleus == "He4" and method == "magnus_adaptive" and smax == 0.0:
            invalid_path = root / "invalid_antihermitian_omega.bin"
            invalid_modelspace = pyIMSRG.ModelSpace(2, nucleus, nucleus)
            invalid_omega = pyIMSRG.Operator(invalid_modelspace)
            invalid_omega.ReadBinary(str(sorted(omega_files)[0]))
            invalid_omega.SetOneBodyME(0, 0, 1e-3)
            invalid_omega.WriteBinary(str(invalid_path))
            rejected = subprocess.run(
                [
                    str(omega_validator),
                    "--emax", "2",
                    "--hw", "20",
                    "--reference", nucleus,
                    str(invalid_path),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            rejected_report = json.loads(rejected.stdout.splitlines()[-1])
            if rejected.returncode == 0 or rejected_report["passed"]:
                raise RuntimeError("Omega validator accepted a numerical symmetry violation")
    return output


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: UnitTestMRSRDriver.py /path/to/imsrg++")
    executable = Path(sys.argv[1]).resolve()
    fixtures = (
        ("He4", "He4_Nrefmax0_final"),
        ("O16", "O16_Nrefmax0_final"),
    )
    with tempfile.TemporaryDirectory() as temporary_directory:
        root = Path(temporary_directory)
        for nucleus, fixture in fixtures:
            data = load_reference(REPOSITORY / "prototype/mrimsrg/data" / fixture)
            A = int(data.metadata["A"])
            reference_file = root / f"{nucleus}.jref"
            summary = export_reference(
                REPOSITORY / "prototype/mrimsrg/data" / fixture,
                reference_file,
            )
            assert summary["Nrefmax"] == 0
            assert summary["natural_basis_is_identity"]

            modelspace = pyIMSRG.ModelSpace(2, nucleus, nucleus)
            modelspace.SetHbarOmega(20.0)
            input_file = root / f"{nucleus}_bare.jcoupled64"
            pyIMSRG.ReadWrite().Write_jcoupled64(
                str(input_file), vacuum_operator(modelspace, data)
            )

            density_modelspace = pyIMSRG.ModelSpace(2, nucleus, nucleus)
            density_modelspace.SetHbarOmega(20.0)
            density_modelspace.SetReferenceOcc(
                pyIMSRG.MRReference.ReadOccupationMap(
                    density_modelspace, str(reference_file)
                )
            )
            reference = pyIMSRG.MRReference.ReadBinary(
                density_modelspace, str(reference_file)
            )
            assert reference.Lambda2.Norm() == 0.0

            for method in ("flow_RK4", "magnus", "magnus_adaptive"):
                for smax in (0.0, 1e-4):
                    mr_output = run_driver(
                        executable,
                        root,
                        nucleus,
                        A,
                        input_file,
                        reference_file,
                        smax,
                        True,
                        method=method,
                    )
                    sr_output = run_driver(
                        executable,
                        root,
                        nucleus,
                        A,
                        input_file,
                        reference_file,
                        smax,
                        False,
                        method=method,
                    )
                    errors = payload_error(sr_output, mr_output)
                    print(
                        f"{nucleus} SR degeneration {method} s={smax:.1e}: "
                        f"zero={errors[0]:.3e} one={errors[1]:.3e} "
                        f"two={errors[2]:.3e}"
                    )
                    assert max(errors) < 1e-10


if __name__ == "__main__":
    main()
