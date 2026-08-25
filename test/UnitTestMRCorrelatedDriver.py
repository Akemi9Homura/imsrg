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
from prototype.mrimsrg.commutator import (  # noqa: E402
    commutator,
    commutator_contributions,
)
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
from prototype.mrimsrg.jcoupling import couple_scalar_two_body  # noqa: E402
from prototype.mrimsrg.sr_imsrgpp_check import (  # noqa: E402
    _mr_rk4_step,
    operator_to_mscheme,
)


def mr_operator_errors(operator, expected, orbits):
    actual = operator_to_mscheme(operator, orbits)
    return (
        abs(expected.zero_body - actual.zero_body),
        float(np.max(np.abs(expected.one_body - actual.one_body))),
        float(np.max(np.abs(expected.two_body - actual.two_body))),
    )


def coupled_term_errors(operator, expected, data, modelspace):
    zero_error = abs(float(operator.ZeroBody) - expected.zero_body)
    one_error = 0.0
    for p, row_p in enumerate(data.orbits):
        a, _, lp, jp, mp, tzp = (int(value) for value in row_p)
        for q, row_q in enumerate(data.orbits):
            b, _, lq, jq, mq, tzq = (int(value) for value in row_q)
            actual = 0.0
            if (lp, jp, mp, tzp) == (lq, jq, mq, tzq):
                actual = operator.OneBody(a, b)
            one_error = max(one_error, abs(actual - expected.one_body[p, q]))

    two_error = 0.0
    expected_blocks = couple_scalar_two_body(expected.two_body, data.orbits)
    for block in expected_blocks:
        channel = modelspace.GetTwoBodyChannelIndex(
            block.J, block.parity, block.Tz
        )
        two_body_channel = modelspace.GetTwoBodyChannel(channel)
        for ibra, (a, b) in enumerate(block.pairs):
            local_bra = two_body_channel.GetLocalIndex(a, b)
            for iket, (c, d) in enumerate(block.pairs):
                local_ket = two_body_channel.GetLocalIndex(c, d)
                actual = operator.TwoBody.GetTBME_norm_chij(
                    channel, channel, local_bra, local_ket
                )
                two_error = max(
                    two_error, abs(actual - block.matrix[ibra, iket])
                )
    return zero_error, one_error, two_error


def one_body_matrix_error(matrix, expected, data):
    maximum = 0.0
    for p, row_p in enumerate(data.orbits):
        a, _, lp, jp, mp, tzp = (int(value) for value in row_p)
        for q, row_q in enumerate(data.orbits):
            b, _, lq, jq, mq, tzq = (int(value) for value in row_q)
            actual = 0.0
            if (lp, jp, mp, tzp) == (lq, jq, mq, tzq):
                actual = matrix(a, b)
            maximum = max(maximum, abs(actual - expected[p, q]))
    return maximum


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

            expected_terms = commutator_contributions(
                expected_eta, expected_initial, natural.densities
            )
            named_maximum = 0.0
            for name in (
                "comm110ss",
                "comm220ss",
                "comm111ss",
                "comm121ss",
                "comm221ss",
                "comm122ss",
                "comm222_pp_hhss",
                "comm222_phss",
            ):
                term = pyIMSRG.Operator(modelspace)
                term.SetHermitian()
                getattr(pyIMSRG.Commutator, name)(actual_eta, initial, term)
                errors = coupled_term_errors(
                    term, expected_terms[name], data, modelspace
                )
                named_maximum = max(named_maximum, *errors)

            mr_one_body = pyIMSRG.MR_comm221_lambda2(
                actual_eta, initial, reference
            ).Total()
            named_maximum = max(
                named_maximum,
                one_body_matrix_error(
                    mr_one_body,
                    expected_terms["mr_lambda2_one_body"].one_body,
                    data,
                ),
            )
            sr_rhs = pyIMSRG.Commutator.Commutator(actual_eta, initial)
            named_maximum = max(
                named_maximum,
                abs(
                    reference.ContractLambda2(sr_rhs.TwoBody)
                    - expected_terms["mr_lambda2_zero_body"].zero_body
                ),
            )
            print(f"{nucleus} named RHS max={named_maximum:.3e}")
            assert named_maximum < 1e-10

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

            groups = spherical_orbit_groups_from_orbits(data.orbits)
            quanta = oscillator_quanta_from_orbits(data.orbits)

            def expected_rhs_at(hamiltonian):
                eta = white_generator(
                    hamiltonian,
                    natural.densities,
                    quanta,
                    spherical_orbit_groups=groups,
                )
                return commutator(eta, hamiltonian, natural.densities)

            checkpoint_step = 1e-3
            expected_checkpoint = expected_initial
            for checkpoint_index in range(1, 4):
                expected_checkpoint = _mr_rk4_step(
                    expected_checkpoint, checkpoint_step, expected_rhs_at
                )
                checkpoint_s = checkpoint_index * checkpoint_step
                actual_path = run_driver(
                    executable,
                    root,
                    nucleus,
                    A,
                    input_file,
                    reference_file,
                    checkpoint_s,
                    True,
                    step=checkpoint_step,
                )
                actual_vacuum = pyIMSRG.Operator(modelspace)
                pyIMSRG.ReadWrite().Read_jcoupled64(
                    str(actual_path), actual_vacuum
                )
                actual_checkpoint = reference.NormalOrder(
                    actual_vacuum.TransformOneAndTwoBody(
                        reference.NaturalOrbitTransformation
                    )
                )
                h_errors = mr_operator_errors(
                    actual_checkpoint, expected_checkpoint, data.orbits
                )

                expected_checkpoint_eta = white_generator(
                    expected_checkpoint,
                    natural.densities,
                    quanta,
                    spherical_orbit_groups=groups,
                )
                actual_checkpoint_eta = pyIMSRG.Operator(modelspace)
                actual_checkpoint_eta.SetAntiHermitian()
                generator.Update(actual_checkpoint, actual_checkpoint_eta)
                eta_checkpoint_errors = mr_operator_errors(
                    actual_checkpoint_eta,
                    expected_checkpoint_eta,
                    data.orbits,
                )

                expected_checkpoint_rhs = commutator(
                    expected_checkpoint_eta,
                    expected_checkpoint,
                    natural.densities,
                )
                actual_checkpoint_rhs = pyIMSRG.MR_Commutator(
                    actual_checkpoint_eta, actual_checkpoint, reference
                )
                rhs_checkpoint_errors = mr_operator_errors(
                    actual_checkpoint_rhs,
                    expected_checkpoint_rhs,
                    data.orbits,
                )
                print(
                    f"{nucleus} checkpoint s={checkpoint_s:.3f}: "
                    f"H={max(h_errors):.3e} "
                    f"eta={max(eta_checkpoint_errors):.3e} "
                    f"RHS={max(rhs_checkpoint_errors):.3e}"
                )
                assert max(h_errors) < 2e-10
                assert max(eta_checkpoint_errors) < 2e-10
                assert max(rhs_checkpoint_errors) < 2e-10


if __name__ == "__main__":
    main()
