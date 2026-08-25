#!/usr/bin/env python3

"""Benchmark the real correlated emax=2 MR RHS and J-scheme storage scaling."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import statistics
import tempfile
import time

import numpy as np

try:
    from .basis import prepare_natural_basis
    from .commutator import commutator
    from .densities import compute_densities
    from .export_jref import export_reference
    from .generator import WHITE_DENOMINATOR_CUTOFF
    from .mr_imsrgpp_flow_check import vacuum_operator
    from .reference_io import load_reference
    from .sr_imsrgpp_check import import_pyimsrg, operator_to_mscheme
except ImportError:
    from basis import prepare_natural_basis
    from commutator import commutator
    from densities import compute_densities
    from export_jref import export_reference
    from generator import WHITE_DENOMINATOR_CUTOFF
    from mr_imsrgpp_flow_check import vacuum_operator
    from reference_io import load_reference
    from sr_imsrgpp_check import import_pyimsrg, operator_to_mscheme


def median_runtime(function, repetitions: int) -> tuple[float, object]:
    samples = []
    result = None
    for _ in range(repetitions):
        start = time.perf_counter()
        result = function()
        samples.append(time.perf_counter() - start)
    return float(statistics.median(samples)), result


def storage_row(pyimsrg, emax: int) -> dict[str, int]:
    modelspace = pyimsrg.ModelSpace(emax, "He4", "He4")
    modelspace.SetHbarOmega(20.0)
    j_orbits = modelspace.GetNumberOrbits()
    m_orbits = sum(modelspace.GetOrbit(i).j2 + 1 for i in range(j_orbits))
    two_body_elements = sum(
        modelspace.GetTwoBodyChannel(channel).GetNumberKets() ** 2
        for channel in range(modelspace.GetNumberTwoBodyChannels())
    )
    return {
        "emax": emax,
        "j_orbits": j_orbits,
        "m_orbits": m_orbits,
        "j_operator_bytes": 8 * (1 + j_orbits**2 + two_body_elements),
        "dense_m_two_body_bytes": 8 * m_orbits**4,
    }


def run(args) -> dict[str, object]:
    pyimsrg = import_pyimsrg(args.pyimsrg_dir)
    reference = load_reference(args.reference, interaction_path=args.interaction)
    if (reference.metadata["A"], reference.metadata["Z"], reference.metadata["Nrefmax"]) != (4, 2, 2):
        raise ValueError("the timing gate requires the real He4 Nrefmax=2 fixture")
    densities = compute_densities(reference.determinants, reference.coefficients)
    natural = prepare_natural_basis(densities)

    with tempfile.TemporaryDirectory() as temporary:
        jref = Path(temporary) / "He4_Nrefmax2.jref"
        export_reference(
            args.reference,
            jref,
            interaction_path=args.interaction,
        )
        modelspace = pyimsrg.ModelSpace(2, "He4", "He4")
        modelspace.SetHbarOmega(20.0)
        modelspace.SetReferenceOcc(
            pyimsrg.MRReference.ReadOccupationMap(modelspace, str(jref))
        )
        mr_reference = pyimsrg.MRReference.ReadBinary(modelspace, str(jref))
        vacuum = vacuum_operator(pyimsrg, modelspace, reference)
        h_j = mr_reference.NormalOrder(
            vacuum.TransformOneAndTwoBody(mr_reference.NaturalOrbitTransformation)
        )
        eta_j = pyimsrg.Operator(modelspace)
        eta_j.SetAntiHermitian()
        generator = pyimsrg.Generator()
        generator.SetType("white-ncsm")
        generator.SetDenominatorPartitioning("Epstein_Nesbet")
        generator.SetDenominatorCutoff(WHITE_DENOMINATOR_CUTOFF)
        generator.Update(h_j, eta_j)

        pyimsrg.MR_Commutator(eta_j, h_j, mr_reference)
        cpp_seconds, cpp_rhs = median_runtime(
            lambda: pyimsrg.MR_Commutator(eta_j, h_j, mr_reference),
            args.repetitions,
        )
        h_m = operator_to_mscheme(h_j, reference.orbits)
        eta_m = operator_to_mscheme(eta_j, reference.orbits)
        commutator(eta_m, h_m, natural.densities)
        python_seconds, python_rhs = median_runtime(
            lambda: commutator(eta_m, h_m, natural.densities),
            args.repetitions,
        )
        cpp_rhs_m = operator_to_mscheme(cpp_rhs, reference.orbits)
        rhs_max_abs = max(
            abs(cpp_rhs_m.zero_body - python_rhs.zero_body),
            float(np.max(np.abs(cpp_rhs_m.one_body - python_rhs.one_body))),
            float(np.max(np.abs(cpp_rhs_m.two_body - python_rhs.two_body))),
        )
        cpp_operator_input_bytes = int(h_j.Size() + eta_j.Size())
        cpp_reference_bytes = int(mr_reference.DataSize())
        python_operator_input_bytes = int(
            h_m.one_body.nbytes
            + h_m.two_body.nbytes
            + eta_m.one_body.nbytes
            + eta_m.two_body.nbytes
        )
        python_reference_bytes = int(
            natural.densities.gamma1.nbytes + natural.densities.lambda2.nbytes
        )

    return {
        "schema": "mrimsrg_jscheme_performance_v1",
        "case": "He4_Nrefmax2_NNLOopt_hw20_emax2_e2max4",
        "omp_num_threads": int(os.environ.get("OMP_NUM_THREADS", "1")),
        "repetitions": args.repetitions,
        "rhs_max_abs_mev": rhs_max_abs,
        "timing_seconds_median": {
            "cpp_jscheme": cpp_seconds,
            "python_dense_mscheme": python_seconds,
            "speedup": python_seconds / cpp_seconds,
        },
        "input_storage_bytes": {
            "cpp_two_operators": cpp_operator_input_bytes,
            "cpp_reference": cpp_reference_bytes,
            "cpp_total": cpp_operator_input_bytes + cpp_reference_bytes,
            "python_two_operators": python_operator_input_bytes,
            "python_reference_gamma1_lambda2": python_reference_bytes,
            "python_total": python_operator_input_bytes + python_reference_bytes,
            "ratio": (python_operator_input_bytes + python_reference_bytes)
            / (cpp_operator_input_bytes + cpp_reference_bytes),
        },
        "storage_scaling": [storage_row(pyimsrg, emax) for emax in range(2, 15, 2)],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--interaction", required=True, type=Path)
    parser.add_argument("--pyimsrg-dir", required=True, type=Path)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    if args.repetitions < 1:
        parser.error("--repetitions must be positive")
    report = run(args)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if args.json is not None:
        if args.json.exists():
            raise FileExistsError(f"refusing to overwrite report: {args.json}")
        args.json.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
