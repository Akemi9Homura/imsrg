#!/usr/bin/env python3

"""Resolve a slow White-NCSM tail into its selected J-scheme channels.

The input Hamiltonians are ordinary vacuum-normal-ordered ``jcoupled64``
files in the original HO orbit order.  Each one is transformed to the
reference's spherical natural basis, MR normal ordered, and passed through
the exact production ``Generator::ConstructGenerator_WhiteNCSM`` path.

Besides the generator norm, the report separates

* the masked, unweighted Hamiltonian matrix elements;
* the lambda-free White-NCSM numerator ``D-D^dagger``; and
* the denominator-weighted anti-Hermitian generator.

This distinction is essential when small fractional occupations make a
perfectly ordinary Hamiltonian matrix element decay very slowly under the
published White-NCSM generator.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

try:
    from .basis import prepare_natural_basis, transform_hamiltonian
    from .densities import compute_densities, validate_densities
    from .generator import (
        masked_decoupling_residual,
        masked_residual_norm,
        oscillator_quanta_from_orbits,
    )
    from .normal_order import VacuumHamiltonian, normal_order
    from .reference_io import load_reference
except ImportError:
    from basis import prepare_natural_basis, transform_hamiltonian
    from densities import compute_densities, validate_densities
    from generator import (
        masked_decoupling_residual,
        masked_residual_norm,
        oscillator_quanta_from_orbits,
    )
    from normal_order import VacuumHamiltonian, normal_order
    from reference_io import load_reference


REPOSITORY = Path(__file__).resolve().parents[2]
WHITE_DENOMINATOR_CUTOFF_MEV = 1.0e-6


def import_pyimsrg(path: Path | None):
    if path is not None:
        sys.path.insert(0, str(path.resolve()))
    try:
        import pyIMSRG  # pylint: disable=import-outside-toplevel
    except ImportError as error:
        raise RuntimeError(
            "cannot import pyIMSRG; source ./sourceme.sh or pass --pyimsrg-dir"
        ) from error
    return pyIMSRG


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_labeled_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("expected LABEL=/path/to/H.jcoupled64")
    label, raw_path = value.split("=", 1)
    if not label:
        raise argparse.ArgumentTypeError("Hamiltonian label must not be empty")
    path = Path(raw_path).expanduser().resolve()
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"Hamiltonian does not exist: {path}")
    return label, path


def orbit_record(modelspace, index: int) -> dict[str, int | float | str]:
    orbit = modelspace.GetOrbit(index)
    letters = "spdfghijklmno"
    letter = letters[orbit.l] if orbit.l < len(letters) else f"l{orbit.l}"
    return {
        "index": index,
        "label": f"{orbit.n}{letter}{orbit.j2}/2_tz{orbit.tz2:+d}",
        "n": orbit.n,
        "l": orbit.l,
        "j2": orbit.j2,
        "tz2": orbit.tz2,
        "oscillator_quanta": 2 * orbit.n + orbit.l,
        "occupation": orbit.occ,
    }


def quantiles(values: list[float]) -> dict[str, float]:
    if not values:
        return {}
    array = np.asarray(values, dtype=np.float64)
    return {
        key: float(value)
        for key, value in zip(
            ("minimum", "p01", "p10", "median", "p90", "p99", "maximum"),
            np.quantile(array, (0.0, 0.01, 0.10, 0.50, 0.90, 0.99, 1.0)),
        )
    }


def add_entry(
    entries: list[dict[str, object]],
    *,
    kind: str,
    eta: float,
    hamiltonian: float,
    forward_weight: float,
    reverse_weight: float,
    forward_denominator: float,
    reverse_denominator: float,
    norm_squared_weight: int,
    quantum_numbers: dict[str, int],
    bra: list[dict[str, object]],
    ket: list[dict[str, object]],
) -> None:
    numerator = hamiltonian * (forward_weight - reverse_weight)
    coefficient = eta / hamiltonian if hamiltonian != 0.0 else 0.0
    entries.append(
        {
            "kind": kind,
            "eta_mev": eta,
            "hamiltonian_mev": hamiltonian,
            "numerator_mev": numerator,
            "forward_weight": forward_weight,
            "reverse_weight": reverse_weight,
            "forward_denominator_mev": forward_denominator,
            "reverse_denominator_mev": reverse_denominator,
            "eta_over_h_inverse_mev": coefficient,
            "inverse_eta_over_h_mev": (
                1.0 / abs(coefficient) if coefficient != 0.0 else math.inf
            ),
            # Both off-diagonal orientations are present in the Operator
            # Frobenius norm.  imsrg++ additionally weights a scalar two-body
            # channel by (2J+1)^2 in TwoBodyME::Norm().
            "eta_norm_squared_contribution": 2.0 * norm_squared_weight * eta**2,
            "numerator_norm_squared_contribution": (
                2.0 * norm_squared_weight * numerator**2
            ),
            "hamiltonian_norm_squared_contribution": (
                2.0 * norm_squared_weight * hamiltonian**2
            ),
            "quantum_numbers": quantum_numbers,
            "bra": bra,
            "ket": ket,
        }
    )


def strict_residual_norm(vacuum, source_reference, natural_basis) -> float:
    """Evaluate the lambda2-dependent masked D-D^dagger norm in m-scheme."""

    try:
        from .sr_imsrgpp_check import operator_to_mscheme
    except ImportError:
        from sr_imsrgpp_check import operator_to_mscheme

    vacuum_m = operator_to_mscheme(vacuum, source_reference.orbits)
    vacuum_nat_m = transform_hamiltonian(
        vacuum_m, natural_basis.vectors, to_natural=True
    )
    hamiltonian_nat_m = normal_order(
        VacuumHamiltonian(
            vacuum_nat_m.zero_body,
            vacuum_nat_m.one_body,
            vacuum_nat_m.two_body,
        ),
        natural_basis.densities,
    )
    residual = masked_decoupling_residual(
        hamiltonian_nat_m,
        natural_basis.densities,
        oscillator_quanta_from_orbits(source_reference.orbits),
    )
    return masked_residual_norm(residual)


def diagnose_stage(
    pyimsrg,
    modelspace,
    reference,
    label: str,
    path: Path,
    top: int,
    *,
    strict_context=None,
):
    vacuum = pyimsrg.Operator(modelspace)
    pyimsrg.ReadWrite().Read_jcoupled64(str(path), vacuum)
    hamiltonian = reference.NormalOrder(
        vacuum.TransformOneAndTwoBody(reference.NaturalOrbitTransformation)
    )
    eta = pyimsrg.Operator(modelspace)
    eta.SetAntiHermitian()
    generator = pyimsrg.Generator()
    generator.SetType("white-ncsm")
    generator.SetDenominatorPartitioning("Epstein_Nesbet")
    generator.SetDenominatorCutoff(WHITE_DENOMINATOR_CUTOFF_MEV)
    generator.Update(hamiltonian, eta)

    entries: list[dict[str, object]] = []
    number_orbits = modelspace.GetNumberOrbits()
    for i in range(number_orbits):
        oi = modelspace.GetOrbit(i)
        for j in range(i + 1, number_orbits):
            oj = modelspace.GetOrbit(j)
            if (oi.l, oi.j2, oi.tz2) != (oj.l, oj.j2, oj.tz2):
                continue
            if 2 * oi.n + oi.l == 2 * oj.n + oj.l:
                continue
            forward = (1.0 - oi.occ) * oj.occ
            reverse = (1.0 - oj.occ) * oi.occ
            add_entry(
                entries,
                kind="one_body",
                eta=eta.GetOneBody(i, j),
                hamiltonian=hamiltonian.GetOneBody(i, j),
                forward_weight=forward,
                reverse_weight=reverse,
                forward_denominator=generator.Get1bDenominatorWhiteNCSM(i, j),
                reverse_denominator=generator.Get1bDenominatorWhiteNCSM(j, i),
                norm_squared_weight=oi.j2 + 1,
                quantum_numbers={"l": oi.l, "j2": oi.j2, "tz2": oi.tz2},
                bra=[orbit_record(modelspace, i)],
                ket=[orbit_record(modelspace, j)],
            )

    for channel_index in range(modelspace.GetNumberTwoBodyChannels()):
        channel = modelspace.GetTwoBodyChannel(channel_index)
        number_kets = channel.GetNumberKets()
        for ibra in range(number_kets):
            bra = channel.GetKet(ibra)
            oi = modelspace.GetOrbit(bra.p)
            oj = modelspace.GetOrbit(bra.q)
            bra_quanta = 2 * oi.n + oi.l + 2 * oj.n + oj.l
            for iket in range(ibra + 1, number_kets):
                ket = channel.GetKet(iket)
                ok = modelspace.GetOrbit(ket.p)
                ol = modelspace.GetOrbit(ket.q)
                ket_quanta = 2 * ok.n + ok.l + 2 * ol.n + ol.l
                if bra_quanta == ket_quanta:
                    continue
                forward = (
                    (1.0 - oi.occ)
                    * (1.0 - oj.occ)
                    * ok.occ
                    * ol.occ
                )
                reverse = (
                    (1.0 - ok.occ)
                    * (1.0 - ol.occ)
                    * oi.occ
                    * oj.occ
                )
                add_entry(
                    entries,
                    kind="two_body",
                    eta=eta.TwoBody.GetTBME_norm_chij(
                        channel_index, channel_index, ibra, iket
                    ),
                    hamiltonian=hamiltonian.TwoBody.GetTBME_norm_chij(
                        channel_index, channel_index, ibra, iket
                    ),
                    forward_weight=forward,
                    reverse_weight=reverse,
                    forward_denominator=generator.Get2bDenominatorWhiteNCSM(
                        channel_index, ibra, iket
                    ),
                    reverse_denominator=generator.Get2bDenominatorWhiteNCSM(
                        channel_index, iket, ibra
                    ),
                    norm_squared_weight=(2 * channel.J + 1) ** 2,
                    quantum_numbers={
                        "channel": channel_index,
                        "J": channel.J,
                        "parity": channel.parity,
                        "Tz": channel.Tz,
                    },
                    bra=[
                        orbit_record(modelspace, bra.p),
                        orbit_record(modelspace, bra.q),
                    ],
                    ket=[
                        orbit_record(modelspace, ket.p),
                        orbit_record(modelspace, ket.q),
                    ],
                )

    eta_norm = eta.Norm()
    eta_reconstructed = math.sqrt(
        sum(float(entry["eta_norm_squared_contribution"]) for entry in entries)
    )
    if abs(eta_reconstructed - eta_norm) > 1.0e-12 * max(1.0, eta_norm):
        raise RuntimeError(
            "selected-channel eta norm does not reproduce Operator::Norm(): "
            f"{eta_reconstructed:.17g} != {eta_norm:.17g}"
        )

    numerator_norm = math.sqrt(
        sum(
            float(entry["numerator_norm_squared_contribution"])
            for entry in entries
        )
    )
    hamiltonian_mask_norm = math.sqrt(
        sum(
            float(entry["hamiltonian_norm_squared_contribution"])
            for entry in entries
        )
    )
    active_entries = [entry for entry in entries if abs(float(entry["eta_mev"])) > 0.0]
    ranked = sorted(
        active_entries,
        key=lambda entry: float(entry["eta_norm_squared_contribution"]),
        reverse=True,
    )
    top_entries = ranked[:top]
    top_fraction = (
        sum(float(entry["eta_norm_squared_contribution"]) for entry in top_entries)
        / eta_norm**2
        if eta_norm > 0.0
        else 0.0
    )
    result = {
        "label": label,
        "path": str(path),
        "sha256": sha256(path),
        "mr_zero_body_mev": hamiltonian.ZeroBody,
        "hamiltonian_norm": hamiltonian.Norm(),
        "selected_unweighted_hamiltonian_norm_mev": hamiltonian_mask_norm,
        "white_ncsm_numerator_norm_mev": numerator_norm,
        "eta_one_body_norm": eta.OneBodyNorm(),
        "eta_two_body_norm": eta.TwoBodyNorm(),
        "eta_norm": eta_norm,
        "eta_norm_reconstructed": eta_reconstructed,
        "selected_channel_count": len(entries),
        "nonzero_eta_channel_count": len(active_entries),
        "absolute_directional_weight_quantiles": quantiles(
            [
                max(
                    abs(float(entry["forward_weight"])),
                    abs(float(entry["reverse_weight"])),
                )
                for entry in active_entries
            ]
        ),
        "absolute_denominator_quantiles_mev": quantiles(
            [
                abs(value)
                for entry in active_entries
                for value in (
                    float(entry["forward_denominator_mev"]),
                    float(entry["reverse_denominator_mev"]),
                )
            ]
        ),
        "inverse_eta_over_h_quantiles_mev": quantiles(
            [
                float(entry["inverse_eta_over_h_mev"])
                for entry in active_entries
                if math.isfinite(float(entry["inverse_eta_over_h_mev"]))
            ]
        ),
        "top_channel_eta_norm_squared_fraction": top_fraction,
        "top_channels": top_entries,
    }
    if strict_context is not None:
        source_reference, natural_basis = strict_context
        result["strict_decoupling_residual_norm_mev"] = strict_residual_norm(
            vacuum, source_reference, natural_basis
        )
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument(
        "--source-reference",
        type=Path,
        help=(
            "original m-scheme NCSM reference bundle; when supplied, also "
            "evaluate the lambda2-dependent strict D-D^dagger diagnostic"
        ),
    )
    parser.add_argument("--nucleus", required=True)
    parser.add_argument("--emax", type=int, required=True)
    parser.add_argument("--hw", type=float, default=20.0)
    parser.add_argument(
        "--hamiltonian",
        type=parse_labeled_path,
        action="append",
        required=True,
        metavar="LABEL=PATH",
    )
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--pyimsrg-dir", type=Path)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    reference_path = args.reference.expanduser().resolve()
    if not reference_path.is_file():
        raise SystemExit(f"reference does not exist: {reference_path}")
    if args.top < 0:
        raise SystemExit("--top must be nonnegative")
    pyimsrg = import_pyimsrg(args.pyimsrg_dir)
    modelspace = pyimsrg.ModelSpace(args.emax, args.nucleus, args.nucleus)
    modelspace.SetHbarOmega(args.hw)
    modelspace.SetReferenceOcc(
        pyimsrg.MRReference.ReadOccupationMap(modelspace, str(reference_path))
    )
    reference = pyimsrg.MRReference.ReadBinary(modelspace, str(reference_path))
    strict_context = None
    source_reference_path = None
    if args.source_reference is not None:
        source_reference_path = args.source_reference.expanduser().resolve()
        source_reference = load_reference(source_reference_path)
        if int(source_reference.metadata["emax"]) != args.emax:
            raise SystemExit("source reference and diagnostic emax differ")
        if float(source_reference.metadata["hw"]) != args.hw:
            raise SystemExit("source reference and diagnostic hw differ")
        densities = compute_densities(
            source_reference.determinants, source_reference.coefficients
        )
        validate_densities(densities, int(source_reference.metadata["A"]))
        strict_context = (
            source_reference,
            prepare_natural_basis(densities),
        )
    stages = [
        diagnose_stage(
            pyimsrg,
            modelspace,
            reference,
            label,
            path,
            args.top,
            strict_context=strict_context,
        )
        for label, path in args.hamiltonian
    ]
    initial = stages[0]
    for stage in stages:
        for field in (
            "selected_unweighted_hamiltonian_norm_mev",
            "white_ncsm_numerator_norm_mev",
            "eta_norm",
        ):
            denominator = float(initial[field])
            stage[f"{field}_ratio_to_first"] = (
                float(stage[field]) / denominator if denominator > 0.0 else 0.0
            )
        if strict_context is not None:
            field = "strict_decoupling_residual_norm_mev"
            denominator = float(initial[field])
            stage[f"{field}_ratio_to_first"] = (
                float(stage[field]) / denominator if denominator > 0.0 else 0.0
            )
    report = {
        "schema": "mrimsrg_white_ncsm_tail_diagnostic_v1",
        "formula": {
            "source": "Vobig 2020 Sec. 6.5.4, Eqs. (6.5.28)-(6.5.34)",
            "implementation": "src/Generator.cc ConstructGenerator_WhiteNCSM",
            "lambda2_in_generator": False,
            "lambda2_in_mr_commutator": True,
            "lambda3": 0,
            "mask": "Delta e != 0 in the connected natural-orbit blocks",
            "denominator": "leading Epstein-Nesbet",
            "denominator_cutoff_mev": WHITE_DENOMINATOR_CUTOFF_MEV,
            "strict_diagnostic": (
                "masked lambda2-dependent D-D^dagger with lambda3=0"
                if strict_context is not None
                else "not requested; pass --source-reference"
            ),
        },
        "reference": {
            "path": str(reference_path),
            "sha256": sha256(reference_path),
            "nucleus": args.nucleus,
            "emax": args.emax,
            "hw_mev": args.hw,
            "occupations": [
                orbit_record(modelspace, index)
                for index in range(modelspace.GetNumberOrbits())
            ],
            "source_reference": (
                str(source_reference_path)
                if source_reference_path is not None
                else None
            ),
        },
        "stages": stages,
    }
    encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is None:
        print(encoded, end="")
    else:
        output = args.output.expanduser().resolve()
        if output.exists():
            raise SystemExit(f"refusing to overwrite existing report: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded, encoding="utf-8")
        print(output)


if __name__ == "__main__":
    main()
