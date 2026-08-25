#!/usr/bin/env python3

"""End-to-end explicit MR driver regression in a real He4 model space."""

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import pyIMSRG


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))

from prototype.mrimsrg.export_jref import export_reference  # noqa: E402
from prototype.mrimsrg.basis import (  # noqa: E402
    prepare_natural_basis,
    transform_array,
)
from prototype.mrimsrg.densities import compute_densities  # noqa: E402
from prototype.mrimsrg.jcoupling import (  # noqa: E402
    couple_scalar_two_body,
    extract_j_orbits,
)
from prototype.mrimsrg.reference_io import load_reference  # noqa: E402
from prototype.mrimsrg.normal_order import (  # noqa: E402
    VacuumHamiltonian,
    normal_order,
)
from prototype.mrimsrg.sr_imsrgpp_check import (  # noqa: E402
    operator_to_mscheme,
    read_jcoupled64,
)


def vacuum_operator(modelspace, data):
    operator = pyIMSRG.Operator(modelspace)
    operator.SetHermitian()
    j_orbits = extract_j_orbits(data.orbits)
    file_to_model = {
        orbit.index: modelspace.GetOrbitIndex(
            orbit.n, orbit.l, orbit.j2, orbit.tz2
        )
        for orbit in j_orbits
    }
    assert all(file_index == model_index for file_index, model_index in file_to_model.items())

    for first in j_orbits:
        for second in j_orbits:
            if (first.l, first.j2, first.tz2) != (
                second.l,
                second.j2,
                second.tz2,
            ):
                continue
            operator.SetOneBody(
                file_to_model[first.index],
                file_to_model[second.index],
                float(data.one_body[first.substates[0], second.substates[0]]),
            )

    for block in couple_scalar_two_body(data.two_body, data.orbits):
        channel = modelspace.GetTwoBodyChannelIndex(
            block.J, block.parity, block.Tz
        )
        two_body_channel = modelspace.GetTwoBodyChannel(channel)
        for ibra, (a, b) in enumerate(block.pairs):
            local_bra = two_body_channel.GetLocalIndex(
                file_to_model[a], file_to_model[b]
            )
            for iket in range(ibra, len(block.pairs)):
                c, d = block.pairs[iket]
                local_ket = two_body_channel.GetLocalIndex(
                    file_to_model[c], file_to_model[d]
                )
                operator.TwoBody.SetTBME_chij(
                    channel,
                    channel,
                    local_bra,
                    local_ket,
                    float(block.matrix[ibra, iket]),
                )
    return operator


def payload_error(expected_path, actual_path):
    expected = read_jcoupled64(expected_path)
    actual = read_jcoupled64(actual_path)
    assert np.array_equal(expected.orbits, actual.orbits)
    expected_two_body = {record[:5]: record[5] for record in expected.records}
    actual_two_body = {record[:5]: record[5] for record in actual.records}
    assert expected_two_body.keys() == actual_two_body.keys()
    return (
        abs(expected.zero_body - actual.zero_body),
        float(np.max(np.abs(expected.one_body - actual.one_body))),
        max(
            abs(expected_two_body[key] - actual_two_body[key])
            for key in expected_two_body
        ),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--imsrg-executable", type=Path, required=True)
    args = parser.parse_args()
    reference_directory = REPOSITORY / "prototype/mrimsrg/data/He4_Nrefmax2"
    data = load_reference(reference_directory)

    with tempfile.TemporaryDirectory() as temporary_directory:
        root = Path(temporary_directory)
        reference_file = root / "He4_Nrefmax2.jref"
        input_file = root / "He4_bare.jcoupled64"
        export_reference(reference_directory, reference_file)

        vacuum_modelspace = pyIMSRG.ModelSpace(2, "He4", "He4")
        vacuum_modelspace.SetHbarOmega(20.0)
        vacuum = vacuum_operator(vacuum_modelspace, data)
        pyIMSRG.ReadWrite().Write_jcoupled64(str(input_file), vacuum)

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
        vacuum_fractional = pyIMSRG.Operator(modelspace)
        pyIMSRG.ReadWrite().Read_jcoupled64(
            str(input_file), vacuum_fractional
        )
        initial = reference.NormalOrder(
            vacuum_fractional.TransformOneAndTwoBody(
                reference.NaturalOrbitTransformation
            )
        )

        # Compare the real NNLOopt He4 E/f/Gamma tensors to the independent
        # dense m-scheme prototype before accepting any driver flow result.
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
                    np.abs(expected_initial.one_body - actual_initial.one_body)
                )
            ),
            float(
                np.max(
                    np.abs(expected_initial.two_body - actual_initial.two_body)
                )
            ),
        )
        print(
            f"MR real E/f/Gamma: zero={initial_errors[0]:.3e} "
            f"one={initial_errors[1]:.3e} two={initial_errors[2]:.3e}"
        )
        assert max(initial_errors) < 2e-10

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
            expected_path = root / f"expected_{label}.jcoupled64"
            actual_path = root / f"driver_{label}.jcoupled64"
            no2bpack_path = root / f"driver_{label}.no2bpack"
            pyIMSRG.ReadWrite().Write_jcoupled64(
                str(expected_path), expected_vacuum
            )

            command = [
                str(args.imsrg_executable),
                f"2bme={input_file}",
                "fmt2=jcoupled64",
                "3bme=none",
                "reference=He4",
                "valence_space=He4",
                "A=4",
                "hw=20",
                "emax=2",
                "basis=oscillator",
                "method=flow_RK4",
                "core_generator=white-ncsm",
                f"mr_reference_file={reference_file}",
                f"smax={smax:.8g}",
                "ds_0=0.0001",
                "dsmax=0.0001",
                "eta_criterion=0",
                "BetaCM=0",
                f"flowfile={root / ('flow_' + label + '.txt')}",
                f"intfile={root / ('int_' + label)}",
                f"write_H_jcoupled64={actual_path}",
                f"write_H_no2bpack={no2bpack_path}",
            ]
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
                raise RuntimeError(f"MR driver {label} failed")
            assert no2bpack_path.stat().st_size > 0
            errors = payload_error(expected_path, actual_path)
            print(
                f"MR driver {label}: zero={errors[0]:.3e} "
                f"one={errors[1]:.3e} two={errors[2]:.3e}"
            )
            assert max(errors) < 2e-11

        mismatched = [
            "reference=O16" if item == "reference=He4" else
            "valence_space=O16" if item == "valence_space=He4" else item
            for item in command
        ]
        rejected = subprocess.run(
            mismatched,
            text=True,
            capture_output=True,
            env=environment,
            check=False,
        )
        assert rejected.returncode != 0
        assert "MR reference A/Z does not match the driver target" in rejected.stderr


if __name__ == "__main__":
    main()
