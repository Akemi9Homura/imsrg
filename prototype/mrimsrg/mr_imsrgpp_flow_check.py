#!/usr/bin/env python3

"""Compare a materialized Python MR flow with the production C++ J-scheme flow."""

import argparse
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

try:
    from .basis import prepare_natural_basis, transform_hamiltonian
    from .commutator import commutator
    from .densities import compute_densities, validate_densities
    from .export_jref import export_reference
    from .generator import (
        WHITE_DENOMINATOR_CUTOFF,
        oscillator_quanta_from_orbits,
        spherical_orbit_groups_from_orbits,
        white_generator,
    )
    from .jcoupling import couple_scalar_two_body, extract_j_orbits
    from .normal_order import MRHamiltonian, VacuumHamiltonian, normal_order, to_vacuum
    from .reference_io import load_reference
    from .sr_imsrgpp_check import _operator_comparison, import_pyimsrg, operator_to_mscheme
except ImportError:
    from basis import prepare_natural_basis, transform_hamiltonian
    from commutator import commutator
    from densities import compute_densities, validate_densities
    from export_jref import export_reference
    from generator import (
        WHITE_DENOMINATOR_CUTOFF,
        oscillator_quanta_from_orbits,
        spherical_orbit_groups_from_orbits,
        white_generator,
    )
    from jcoupling import couple_scalar_two_body, extract_j_orbits
    from normal_order import MRHamiltonian, VacuumHamiltonian, normal_order, to_vacuum
    from reference_io import load_reference
    from sr_imsrgpp_check import _operator_comparison, import_pyimsrg, operator_to_mscheme


NUCLEI = {"He4": (4, 2), "Be8": (8, 4), "C12": (12, 6), "O16": (16, 8)}


def vacuum_operator(pyimsrg, modelspace, reference):
    operator = pyimsrg.Operator(modelspace)
    operator.SetHermitian()
    j_orbits = extract_j_orbits(reference.orbits)
    file_to_model = {
        orbit.index: modelspace.GetOrbitIndex(orbit.n, orbit.l, orbit.j2, orbit.tz2)
        for orbit in j_orbits
    }
    for first in j_orbits:
        for second in j_orbits:
            if (first.l, first.j2, first.tz2) != (second.l, second.j2, second.tz2):
                continue
            operator.SetOneBody(
                file_to_model[first.index],
                file_to_model[second.index],
                float(reference.one_body[first.substates[0], second.substates[0]]),
            )
    for block in couple_scalar_two_body(reference.two_body, reference.orbits):
        channel_index = modelspace.GetTwoBodyChannelIndex(
            block.J, block.parity, block.Tz
        )
        channel = modelspace.GetTwoBodyChannel(channel_index)
        for ibra, (a, b) in enumerate(block.pairs):
            local_bra = channel.GetLocalIndex(file_to_model[a], file_to_model[b])
            for iket in range(ibra, len(block.pairs)):
                c, d = block.pairs[iket]
                local_ket = channel.GetLocalIndex(file_to_model[c], file_to_model[d])
                operator.TwoBody.SetTBME_chij(
                    channel_index,
                    channel_index,
                    local_bra,
                    local_ket,
                    float(block.matrix[ibra, iket]),
                )
    return operator


def load_flow(path, reference):
    metadata = json.loads((path / "metadata.json").read_text(encoding="utf-8"))
    if metadata.get("schema") != "mrimsrg_flow_v1":
        raise ValueError("unsupported flow schema")
    for key in ("A", "Z", "Nrefmax"):
        if metadata["reference_metadata"].get(key) != reference.metadata.get(key):
            raise ValueError(f"flow and reference differ in {key}")
    trajectory = metadata.get("trajectory", [])
    if not trajectory:
        raise ValueError("flow trajectory is empty")
    target_s = float(trajectory[-1]["s"])
    if not np.isclose(target_s, metadata["flow_settings"]["smax"], atol=1e-10):
        raise ValueError("Python flow did not reach its fixed target s")
    if metadata.get("flow_converged"):
        raise ValueError("comparison flow stopped early instead of at fixed s")
    final = MRHamiltonian(
        float(metadata["final_mr_zero_body"]),
        np.load(path / "final_mr_one_body.npy", allow_pickle=False),
        np.load(path / "final_mr_two_body.npy", allow_pickle=False),
    )
    vacuum = VacuumHamiltonian(
        float(metadata["final_vacuum_zero_body"]),
        np.load(path / "final_vacuum_one_body.npy", allow_pickle=False),
        np.load(path / "final_vacuum_two_body.npy", allow_pickle=False),
    )
    return metadata, target_s, final, vacuum


def maximum(report, names):
    return max(
        report[name][rank]["max_abs"]
        for name in names
        for rank in ("zero_body", "one_body", "two_body")
    )


def run(args):
    pyimsrg = import_pyimsrg(args.pyimsrg_dir)
    reference = load_reference(args.reference, interaction_path=args.interaction)
    if (reference.metadata["A"], reference.metadata["Z"]) != NUCLEI[args.nucleus]:
        raise ValueError("reference A/Z does not match nucleus")
    densities = compute_densities(reference.determinants, reference.coefficients)
    validate_densities(densities, int(reference.metadata["A"]))
    metadata, target_s, production_final, production_vacuum = load_flow(
        args.production_flow, reference
    )
    if not np.isclose(metadata["flow_settings"]["relative_tolerance"], args.ode_tolerance):
        raise ValueError("Python rtol differs from requested C++ tolerance")
    if not np.isclose(metadata["flow_settings"]["absolute_tolerance"], args.ode_tolerance):
        raise ValueError("Python atol differs from requested C++ tolerance")

    with tempfile.TemporaryDirectory() as temporary_directory:
        jref = Path(temporary_directory) / f"{args.nucleus}.jref"
        export_reference(args.reference, jref, interaction_path=args.interaction)
        modelspace = pyimsrg.ModelSpace(2, args.nucleus, args.nucleus)
        modelspace.SetHbarOmega(20.0)
        modelspace.SetReferenceOcc(pyimsrg.MRReference.ReadOccupationMap(modelspace, str(jref)))
        mr_reference = pyimsrg.MRReference.ReadBinary(modelspace, str(jref))
        vacuum_j = vacuum_operator(pyimsrg, modelspace, reference)
        initial_j = mr_reference.NormalOrder(
            vacuum_j.TransformOneAndTwoBody(mr_reference.NaturalOrbitTransformation)
        )

        solver = pyimsrg.IMSRGSolver(initial_j)
        solver.SetMRReference(mr_reference)
        solver.SetMethod("flow")
        solver.SetGenerator("white-ncsm")
        solver.SetDenominatorPartitioning("Epstein_Nesbet")
        solver.SetDenominatorCutoff(WHITE_DENOMINATOR_CUTOFF)
        solver.SetEtaCriterion(0.0)
        solver.SetODETolerance(args.ode_tolerance)
        solver.SetDs(args.initial_step)
        solver.SetSmax(target_s)
        solver.Solve()
        final_nat_j = solver.GetH_s()

        final_ho_j = final_nat_j.TransformOneAndTwoBody(
            mr_reference.NaturalOrbitTransformation.t()
        )
        final_ho_m = operator_to_mscheme(final_ho_j, reference.orbits)
        eta_nat_j = pyimsrg.Operator(modelspace)
        eta_nat_j.SetAntiHermitian()
        generator = pyimsrg.Generator()
        generator.SetType("white-ncsm")
        generator.SetDenominatorPartitioning("Epstein_Nesbet")
        generator.SetDenominatorCutoff(WHITE_DENOMINATOR_CUTOFF)
        generator.Update(final_nat_j, eta_nat_j)
        rhs_nat_j = pyimsrg.MR_Commutator(eta_nat_j, final_nat_j, mr_reference)
        eta_ho_m = operator_to_mscheme(
            eta_nat_j.TransformOneAndTwoBody(mr_reference.NaturalOrbitTransformation.t()),
            reference.orbits,
        )
        rhs_ho_m = operator_to_mscheme(
            rhs_nat_j.TransformOneAndTwoBody(mr_reference.NaturalOrbitTransformation.t()),
            reference.orbits,
        )
        quanta = oscillator_quanta_from_orbits(reference.orbits)
        groups = spherical_orbit_groups_from_orbits(reference.orbits)
        natural = prepare_natural_basis(densities)
        production_final_nat = (
            production_final
            if natural.is_identity
            else transform_hamiltonian(
                production_final, natural.vectors, to_natural=True
            )
        )
        production_eta_nat = white_generator(
            production_final_nat,
            natural.densities,
            quanta,
            spherical_orbit_groups=groups,
        )
        production_rhs_nat = commutator(
            production_eta_nat, production_final_nat, natural.densities
        )
        production_eta = (
            production_eta_nat
            if natural.is_identity
            else transform_hamiltonian(
                production_eta_nat, natural.vectors, to_natural=False
            )
        )
        production_rhs = (
            production_rhs_nat
            if natural.is_identity
            else transform_hamiltonian(
                production_rhs_nat, natural.vectors, to_natural=False
            )
        )
        vacuum_final_j = mr_reference.UndoNormalOrder(final_nat_j).TransformOneAndTwoBody(
            mr_reference.NaturalOrbitTransformation.t()
        )
        vacuum_final_m = operator_to_mscheme(vacuum_final_j, reference.orbits)
        if args.output_jcoupled64 is not None:
            pyimsrg.ReadWrite().Write_jcoupled64(
                str(args.output_jcoupled64), vacuum_final_j
            )

    initial_vacuum = VacuumHamiltonian(
        float(reference.metadata["zero_body"]), reference.one_body, reference.two_body
    )
    initial_mr = normal_order(initial_vacuum, densities)
    report = {
        "schema": "mrimsrg_cpp_jscheme_full_flow_v1",
        "nucleus": args.nucleus,
        "nrefmax": int(reference.metadata["Nrefmax"]),
        "target_s": target_s,
        "ode_tolerance": args.ode_tolerance,
        "python_ode_method": metadata["ode_method"],
        "cpp_ode_method": "boost::odeint runge_kutta_dopri5",
        "initial_h": _operator_comparison(
            operator_to_mscheme(
                initial_j.TransformOneAndTwoBody(mr_reference.NaturalOrbitTransformation.t()),
                reference.orbits,
            ),
            initial_mr,
        ),
        "h": _operator_comparison(final_ho_m, production_final),
        "eta": _operator_comparison(eta_ho_m, production_eta),
        "rhs": _operator_comparison(rhs_ho_m, production_rhs),
        "vacuum_h": _operator_comparison(vacuum_final_m, production_vacuum),
        "python_final_zero_body": production_final.zero_body,
        "cpp_final_zero_body": final_ho_m.zero_body,
        "cpp_vacuum_jcoupled64": (
            None if args.output_jcoupled64 is None else str(args.output_jcoupled64.resolve())
        ),
    }
    report["maximum_full_flow_abs_mev"] = maximum(
        report, ("h", "eta", "rhs", "vacuum_h")
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nucleus", choices=tuple(NUCLEI), required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--interaction", type=Path, required=True)
    parser.add_argument("--production-flow", type=Path, required=True)
    parser.add_argument("--pyimsrg-dir", type=Path, required=True)
    parser.add_argument("--ode-tolerance", type=float, default=1e-9)
    parser.add_argument("--initial-step", type=float, default=1e-2)
    parser.add_argument("--flow-tolerance", type=float, default=1e-5)
    parser.add_argument("--output-jcoupled64", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    report = run(args)
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.json is not None:
        if args.json.exists():
            raise FileExistsError(f"refusing to overwrite existing report: {args.json}")
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered + "\n", encoding="utf-8")
    if report["maximum_full_flow_abs_mev"] > args.flow_tolerance:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
