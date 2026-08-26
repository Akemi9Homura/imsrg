#!/usr/bin/env python3

"""Finite-flow MR-Magnus step/BCH and downstream spectral stability gate."""

from pathlib import Path
import re
import subprocess
import sys
import tempfile

import numpy as np
import pyIMSRG


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))
sys.path.insert(0, str(REPOSITORY / "test"))

from UnitTestMRDriver import payload_error, vacuum_operator  # noqa: E402
from UnitTestMRSRDriver import run_driver  # noqa: E402
from prototype.mrimsrg.export_jref import export_reference  # noqa: E402
from prototype.mrimsrg.reference_io import load_reference  # noqa: E402


INTERACTION = Path(
    "/home/mengziyan/Forces/N2LO_opt/"
    "TwBME_N2LO_opt_hw20_emax2_e2max4.minipack"
)
ENERGY_PATTERN = re.compile(r"^state=(\d+) E=([-+0-9.eE]+)", re.MULTILINE)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def run_solver(initial, reference, *, dsmax, ode_tolerance,
               transform_threshold, product_threshold, smax):
    solver = pyIMSRG.IMSRGSolver(initial)
    solver.SetMRReference(reference)
    solver.SetGenerator("white-ncsm")
    solver.SetMethod("magnus_adaptive")
    solver.SetEtaCriterion(0.0)
    solver.SetDs(0.5)
    solver.SetDsmax(dsmax)
    solver.SetODETolerance(ode_tolerance)
    solver.SetOmegaNormMax(0.25)
    solver.SetSmax(smax)
    pyIMSRG.BCH.Set_BCH_Transform_Threshold(transform_threshold)
    pyIMSRG.BCH.Set_BCH_Product_Threshold(product_threshold)
    try:
        solver.Solve()
    finally:
        pyIMSRG.BCH.Set_BCH_Transform_Threshold(1e-9)
        pyIMSRG.BCH.Set_BCH_Product_Threshold(1e-4)
    require(abs(solver.GetS() - smax) < 1e-12, "Magnus flow did not reach smax")
    require(solver.GetH_s().IsHermitian(), "MR-Magnus H lost Hermitian metadata")
    omega = tuple(solver.GetOmega())
    require(omega, "MR-Magnus did not retain an Omega segment")
    require(all(operator.IsAntiHermitian() for operator in omega),
            "MR-Magnus Omega lost anti-Hermitian metadata")
    require(all(abs(float(operator.ZeroBody)) < 1e-12 for operator in omega),
            "anti-Hermitian scalar Omega acquired a real zero-body term")
    return solver


def materialize(root, label, solver, reference):
    vacuum = reference.UndoNormalOrder(solver.GetH_s()).TransformOneAndTwoBody(
        reference.NaturalOrbitTransformation.t()
    )
    require(vacuum.IsHermitian(), "vacuum Hamiltonian lost Hermitian metadata")
    j64 = root / f"{label}.jcoupled64"
    packed = root / f"{label}.no2bpack"
    pyIMSRG.ReadWrite().Write_jcoupled64(str(j64), vacuum)
    exporter = REPOSITORY / "prototype/mrimsrg/build/mrimsrg_export_no2bpack"
    completed = subprocess.run(
        [
            str(exporter),
            "--interaction", str(INTERACTION),
            "--jcoupled64", str(j64),
            "--output", str(packed),
            "--Z", "2",
            "--N", "2",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    require(completed.returncode == 0, completed.stdout + completed.stderr)
    require(j64.stat().st_size > 0 and packed.stat().st_size > 0,
            "Hamiltonian materialization produced an empty file")
    return j64, packed


def diagonalize(path, *, packed):
    validator = REPOSITORY / "prototype/mrimsrg/build/mrimsrg_validate"
    selector = "--no2bpack" if packed else "--jcoupled64"
    command = [str(validator), selector, str(path)]
    if not packed:
        command.extend(("--interaction", str(INTERACTION)))
    command.extend(("--Z", "2", "--N", "2", "--nmax", "8", "--states", "3"))
    completed = subprocess.run(
        command, text=True, capture_output=True, check=False
    )
    require(completed.returncode == 0, completed.stdout + completed.stderr)
    energies = tuple(
        float(energy) for _, energy in ENERGY_PATTERN.findall(completed.stdout)
    )
    require(len(energies) == 3, "NCSM readback did not return three states")
    return energies


def main():
    if len(sys.argv) > 3:
        raise SystemExit(
            "usage: UnitTestMRMagnusStability.py [smax [imsrg-executable]]"
        )
    smax = float(sys.argv[1]) if len(sys.argv) == 2 else 1.0
    if len(sys.argv) == 3:
        smax = float(sys.argv[1])
        executable = Path(sys.argv[2]).resolve()
    else:
        executable = None
    require(INTERACTION.is_file(), f"missing interaction: {INTERACTION}")
    reference_directory = REPOSITORY / "prototype/mrimsrg/data/He4_Nrefmax2"
    data = load_reference(reference_directory)

    with tempfile.TemporaryDirectory() as temporary_directory:
        root = Path(temporary_directory)
        reference_file = root / "He4_Nrefmax2.jref"
        export_reference(reference_directory, reference_file)

        modelspace = pyIMSRG.ModelSpace(2, "He4", "He4")
        modelspace.SetHbarOmega(20.0)
        modelspace.SetReferenceOcc(
            pyIMSRG.MRReference.ReadOccupationMap(
                modelspace, str(reference_file)
            )
        )
        reference = pyIMSRG.MRReference.ReadBinary(
            modelspace, str(reference_file)
        )
        vacuum_input = vacuum_operator(modelspace, data)
        input_file = root / "He4_bare.jcoupled64"
        pyIMSRG.ReadWrite().Write_jcoupled64(str(input_file), vacuum_input)
        initial = reference.NormalOrder(
            vacuum_input.TransformOneAndTwoBody(
                reference.NaturalOrbitTransformation
            )
        )

        default = run_solver(
            initial, reference, dsmax=0.5, ode_tolerance=1e-6,
            transform_threshold=1e-9, product_threshold=1e-4, smax=smax,
        )
        tight = run_solver(
            initial, reference, dsmax=0.05, ode_tolerance=1e-7,
            transform_threshold=1e-10, product_threshold=1e-5, smax=smax,
        )
        default_j64, default_packed = materialize(
            root, "default", default, reference
        )
        tight_j64, tight_packed = materialize(root, "tight", tight, reference)
        if executable is not None:
            driver_j64 = run_driver(
                executable, root, "He4", 4, input_file, reference_file,
                smax, True, step=0.5, method="magnus_adaptive",
            )
            driver_errors = payload_error(default_j64, driver_j64)
            print(
                "adaptive Magnus driver-vs-library H: "
                f"zero={driver_errors[0]:.3e} one={driver_errors[1]:.3e} "
                f"two={driver_errors[2]:.3e} MeV"
            )
            require(max(driver_errors) < 2e-10,
                    "adaptive Magnus driver differs from the library path")
        errors = payload_error(default_j64, tight_j64)
        print(
            f"s={smax:g} default-vs-tight H: zero={errors[0]:.6e} "
            f"one={errors[1]:.6e} two={errors[2]:.6e} MeV"
        )

        default_j64_energies = diagonalize(default_j64, packed=False)
        tight_j64_energies = diagonalize(tight_j64, packed=False)
        default_packed_energies = diagonalize(default_packed, packed=True)
        tight_packed_energies = diagonalize(tight_packed, packed=True)
        numerical_spectral_error = max(
            abs(left - right)
            for left, right in zip(default_j64_energies, tight_j64_energies)
        )
        default_packing_error = max(
            abs(left - right)
            for left, right in zip(default_j64_energies, default_packed_energies)
        )
        tight_packing_error = max(
            abs(left - right)
            for left, right in zip(tight_j64_energies, tight_packed_energies)
        )
        print(f"default J64 energies={default_j64_energies}")
        print(f"tight   J64 energies={tight_j64_energies}")
        print(
            f"spectral numerical delta={numerical_spectral_error:.6e} MeV; "
            f"packing deltas={default_packing_error:.6e}/"
            f"{tight_packing_error:.6e} MeV"
        )
        require(max(errors) < 1e-3,
                "tenfold Magnus/BCH refinement changed a matrix element by >=1 keV")
        require(numerical_spectral_error < 1e-3,
                "tenfold Magnus/BCH refinement changed an NCSM level by >=1 keV")
        require(max(default_packing_error, tight_packing_error) < 2e-6,
                "no2bpack float32 readback exceeded the established 2 eV window")


if __name__ == "__main__":
    main()
