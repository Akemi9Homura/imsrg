#!/usr/bin/env python3

"""Real correlated-reference E/f/Gamma and short-flow driver gates."""

from pathlib import Path
import sys
import tempfile

import numpy as np
import pyIMSRG


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))
sys.path.insert(0, str(REPOSITORY / "test"))

from UnitTestMRDriver import payload_error, vacuum_operator  # noqa: E402
from UnitTestMRSRDriver import run_driver  # noqa: E402
from prototype.mrimsrg.basis import (  # noqa: E402
    prepare_natural_basis,
    transform_array,
)
from prototype.mrimsrg.densities import compute_densities  # noqa: E402
from prototype.mrimsrg.commutator import commutator  # noqa: E402
from prototype.mrimsrg.export_jref import export_reference  # noqa: E402
from prototype.mrimsrg.generator import (  # noqa: E402
    oscillator_quanta_from_orbits,
    spherical_orbit_groups_from_orbits,
    white_generator,
)
from prototype.mrimsrg.normal_order import (  # noqa: E402
    VacuumHamiltonian,
    normal_order,
)
from prototype.mrimsrg.reference_io import load_reference  # noqa: E402
from prototype.mrimsrg.sr_imsrgpp_check import operator_to_mscheme  # noqa: E402


def main():
    if len(sys.argv) != 2:
        raise SystemExit(
            "usage: UnitTestMRCorrelatedDriver.py /path/to/imsrg++"
        )
    executable = Path(sys.argv[1]).resolve()
    fixtures = (
        ("He4", "He4_Nrefmax2"),
        ("Be8", "Be8_Nrefmax0_final"),
        ("C12", "C12_Nrefmax0_final"),
        ("O16", "O16_Nrefmax2"),
    )
    with tempfile.TemporaryDirectory() as temporary_directory:
        root = Path(temporary_directory)
        for nucleus, fixture in fixtures:
            fixture_path = REPOSITORY / "prototype/mrimsrg/data" / fixture
            data = load_reference(fixture_path)
            A = int(data.metadata["A"])
            reference_file = root / f"{nucleus}.jref"
            export_reference(fixture_path, reference_file)

            vacuum_modelspace = pyIMSRG.ModelSpace(2, nucleus, nucleus)
            vacuum_modelspace.SetHbarOmega(20.0)
            input_file = root / f"{nucleus}_bare.jcoupled64"
            pyIMSRG.ReadWrite().Write_jcoupled64(
                str(input_file), vacuum_operator(vacuum_modelspace, data)
            )

            modelspace = pyIMSRG.ModelSpace(2, nucleus, nucleus)
            modelspace.SetHbarOmega(20.0)
            modelspace.SetReferenceOcc(
                pyIMSRG.MRReference.ReadOccupationMap(
                    modelspace, str(reference_file)
                )
            )
            reference = pyIMSRG.MRReference.ReadBinary(
                modelspace, str(reference_file)
            )
            assert reference.Lambda2.Norm() > 1e-6
            vacuum_fractional = pyIMSRG.Operator(modelspace)
            pyIMSRG.ReadWrite().Read_jcoupled64(
                str(input_file), vacuum_fractional
            )
            initial = reference.NormalOrder(
                vacuum_fractional.TransformOneAndTwoBody(
                    reference.NaturalOrbitTransformation
                )
            )

            natural = prepare_natural_basis(
                compute_densities(data.determinants, data.coefficients)
            )
            expected_initial = normal_order(
                VacuumHamiltonian(
                    0.0,
                    transform_array(
                        data.one_body, natural.vectors, to_natural=True
                    ),
                    transform_array(
                        data.two_body, natural.vectors, to_natural=True
                    ),
                ),
                natural.densities,
            )
            actual_initial = operator_to_mscheme(initial, data.orbits)
            initial_errors = (
                abs(expected_initial.zero_body - actual_initial.zero_body),
                float(
                    np.max(
                        np.abs(
                            expected_initial.one_body - actual_initial.one_body
                        )
                    )
                ),
                float(
                    np.max(
                        np.abs(
                            expected_initial.two_body - actual_initial.two_body
                        )
                    )
                ),
            )
            print(
                f"{nucleus} E/f/Gamma: zero={initial_errors[0]:.3e} "
                f"one={initial_errors[1]:.3e} two={initial_errors[2]:.3e}"
            )
            assert max(initial_errors) < 2e-10

            expected_eta = white_generator(
                expected_initial,
                natural.densities,
                oscillator_quanta_from_orbits(data.orbits),
                spherical_orbit_groups=spherical_orbit_groups_from_orbits(
                    data.orbits
                ),
            )
            actual_eta = pyIMSRG.Operator(modelspace)
            actual_eta.SetAntiHermitian()
            generator = pyIMSRG.Generator()
            generator.SetType("white-ncsm")
            generator.Update(initial, actual_eta)
            actual_eta_m = operator_to_mscheme(actual_eta, data.orbits)
            eta_errors = (
                abs(expected_eta.zero_body - actual_eta_m.zero_body),
                float(
                    np.max(
                        np.abs(expected_eta.one_body - actual_eta_m.one_body)
                    )
                ),
                float(
                    np.max(
                        np.abs(expected_eta.two_body - actual_eta_m.two_body)
                    )
                ),
            )
            print(
                f"{nucleus} eta(s=0): zero={eta_errors[0]:.3e} "
                f"one={eta_errors[1]:.3e} two={eta_errors[2]:.3e}"
            )
            assert max(eta_errors) < 2e-10

            expected_rhs = commutator(
                expected_eta, expected_initial, natural.densities
            )
            actual_rhs = pyIMSRG.MR_Commutator(
                actual_eta, initial, reference
            )
            actual_rhs_m = operator_to_mscheme(actual_rhs, data.orbits)
            rhs_errors = (
                abs(expected_rhs.zero_body - actual_rhs_m.zero_body),
                float(
                    np.max(
                        np.abs(expected_rhs.one_body - actual_rhs_m.one_body)
                    )
                ),
                float(
                    np.max(
                        np.abs(expected_rhs.two_body - actual_rhs_m.two_body)
                    )
                ),
            )
            print(
                f"{nucleus} RHS(s=0): zero={rhs_errors[0]:.3e} "
                f"one={rhs_errors[1]:.3e} two={rhs_errors[2]:.3e}"
            )
            assert max(rhs_errors) < 1e-10

            for label, smax in (("s0", 0.0), ("rk4", 1e-4)):
                if smax == 0.0:
                    final_mr = pyIMSRG.Operator(initial)
                else:
                    solver = pyIMSRG.IMSRGSolver(initial)
                    solver.SetMRReference(reference)
                    solver.SetGenerator("white-ncsm")
                    solver.SetMethod("flow_RK4")
                    solver.SetEtaCriterion(0.0)
                    solver.SetDsmax(1e-4)
                    solver.SetSmax(smax)
                    solver.Solve()
                    final_mr = solver.GetH_s()
                expected_vacuum = reference.UndoNormalOrder(
                    final_mr
                ).TransformOneAndTwoBody(
                    reference.NaturalOrbitTransformation.t()
                )
                expected_path = root / f"{nucleus}_expected_{label}.jcoupled64"
                pyIMSRG.ReadWrite().Write_jcoupled64(
                    str(expected_path), expected_vacuum
                )
                actual_path = run_driver(
                    executable,
                    root,
                    nucleus,
                    A,
                    input_file,
                    reference_file,
                    smax,
                    True,
                )
                errors = payload_error(expected_path, actual_path)
                print(
                    f"{nucleus} driver {label}: zero={errors[0]:.3e} "
                    f"one={errors[1]:.3e} two={errors[2]:.3e}"
                )
                assert max(errors) < 2e-10


if __name__ == "__main__":
    main()
