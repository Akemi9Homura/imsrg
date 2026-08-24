"""Run and materialize one fixed-input MR-IMSRG calculation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

try:
    from .densities import compute_densities, reference_energy, validate_densities
    from .flow import FlowPoint, FlowSettings, integrate_flow
    from .generator import oscillator_quanta_from_orbits
    from .normal_order import (
        MRHamiltonian,
        VacuumHamiltonian,
        normal_order,
        validate_hermitian,
    )
    from .output import save_flow_output
    from .reference_io import load_reference
except ImportError:
    from densities import compute_densities, reference_energy, validate_densities
    from flow import FlowPoint, FlowSettings, integrate_flow
    from generator import oscillator_quanta_from_orbits
    from normal_order import (
        MRHamiltonian,
        VacuumHamiltonian,
        normal_order,
        validate_hermitian,
    )
    from output import save_flow_output
    from reference_io import load_reference


def _print_point(point: FlowPoint) -> None:
    print(
        f"step={point.step:3d} s={point.s:.10g} E={point.zero_body:.12f} "
        f"strict_residual={point.residual:.8e} strict_ratio={point.residual_ratio:.8e} "
        f"generator_ratio={point.generator_residual_ratio:.8e} "
        f"numerator_ratio={point.generator_numerator_residual_ratio:.8e} "
        f"sym={max(point.one_body_hermiticity_error, point.two_body_hermiticity_error, point.two_body_antisymmetry_error):.2e}",
        flush=True,
    )


def _load_resume_state(
    path: Path,
    reference_metadata: dict[str, object],
    densities: Densities,
) -> tuple[MRHamiltonian, float, tuple[float, float, float]]:
    """Load a prior materialized MR state without renormalizing its target."""
    metadata_path = path / "metadata.json"
    if not metadata_path.is_file():
        raise FileNotFoundError(f"resume metadata is unavailable: {metadata_path}")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("schema") != "mrimsrg_flow_v1":
        raise ValueError("resume input is not an MR-IMSRG flow output")
    if metadata.get("reference_metadata") != reference_metadata:
        raise ValueError("resume input was produced from a different reference")
    for name, expected in (
        ("gamma1", densities.gamma1),
        ("gamma2", densities.gamma2),
        ("lambda2", densities.lambda2),
    ):
        actual = np.load(path / f"{name}.npy", allow_pickle=False)
        if not np.array_equal(actual, expected):
            raise ValueError(f"resume input {name} differs from the current reference")
    trajectory = metadata.get("trajectory")
    if not isinstance(trajectory, list) or not trajectory:
        raise ValueError("resume input has no recorded trajectory")
    final_point = trajectory[-1]
    start_s = float(final_point["s"])

    saved = metadata.get("residual_normalization")
    if isinstance(saved, dict):
        normalization = (
            float(saved["strict"]),
            float(saved["generator"]),
            float(saved["generator_numerator"]),
        )
    else:
        first = trajectory[0]
        required = {
            "residual",
            "residual_ratio",
            "generator_residual",
            "generator_residual_ratio",
            "generator_numerator_residual",
            "generator_numerator_residual_ratio",
        }
        if not required.issubset(first):
            raise ValueError(
                "resume input predates the separately recorded generator norms"
            )

        def recover(value_key: str, ratio_key: str) -> float:
            value = float(first[value_key])
            ratio = float(first[ratio_key])
            return value / ratio if ratio > 0.0 else value

        normalization = (
            recover("residual", "residual_ratio"),
            recover("generator_residual", "generator_residual_ratio"),
            recover(
                "generator_numerator_residual",
                "generator_numerator_residual_ratio",
            ),
        )

    one_body = np.load(path / "final_mr_one_body.npy", allow_pickle=False)
    two_body = np.load(path / "final_mr_two_body.npy", allow_pickle=False)
    hamiltonian = MRHamiltonian(
        float(metadata["final_mr_zero_body"]), one_body, two_body
    )
    if one_body.shape != densities.gamma1.shape:
        raise ValueError("resume one-body tensor has an incompatible shape")
    if two_body.shape != densities.gamma2.shape:
        raise ValueError("resume two-body tensor has an incompatible shape")
    if not (
        np.isfinite(hamiltonian.zero_body)
        and np.all(np.isfinite(one_body))
        and np.all(np.isfinite(two_body))
    ):
        raise ValueError("resume Hamiltonian contains non-finite values")
    return hamiltonian, start_s, normalization


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--interaction",
        type=Path,
        help="relocated copy of the fixed interaction, still verified by SHA-256",
    )
    parser.add_argument("--smax", type=float, default=2.0)
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--atol", type=float, default=1e-8)
    parser.add_argument("--max-step", type=float, default=10.0)
    parser.add_argument("--residual-ratio", type=float, default=1e-6)
    parser.add_argument(
        "--resume-from",
        type=Path,
        help="continue from a prior materialized flow without resetting residual ratios",
    )
    parser.add_argument(
        "--checkpoint-s",
        type=float,
        default=None,
        help="materialize one intermediate Hamiltonian at this flow parameter",
    )
    args = parser.parse_args()

    reference = load_reference(args.reference, interaction_path=args.interaction)
    densities = compute_densities(reference.determinants, reference.coefficients)
    validate_densities(densities, int(reference.metadata["A"]))
    contracted = reference_energy(reference, densities)
    if abs(contracted - float(reference.metadata["reference_energy"])) > 1e-8:
        raise ValueError("reference energy does not agree with density contraction")
    vacuum = VacuumHamiltonian(
        float(reference.metadata["zero_body"]),
        reference.one_body,
        reference.two_body,
    )
    initial = normal_order(vacuum, densities)
    start_s = 0.0
    residual_normalization = None
    if args.resume_from is not None:
        initial, start_s, residual_normalization = _load_resume_state(
            args.resume_from.resolve(), reference.metadata, densities
        )
        validate_hermitian(initial)
    settings = FlowSettings(
        smax=args.smax,
        relative_tolerance=args.rtol,
        absolute_tolerance=args.atol,
        max_step=args.max_step,
        residual_ratio=args.residual_ratio,
        checkpoint_s=args.checkpoint_s,
    )
    result = integrate_flow(
        initial,
        densities,
        oscillator_quanta_from_orbits(reference.orbits),
        settings,
        observer=_print_point,
        start_s=start_s,
        residual_normalization=residual_normalization,
    )
    final_vacuum = save_flow_output(
        args.output,
        args.reference,
        reference,
        densities,
        initial,
        result,
        settings,
        resumed_from=args.resume_from,
    )
    validate_hermitian(final_vacuum, tolerance=settings.symmetry_tolerance)
    print(
        f"saved {args.output}: converged={result.converged} "
        f"nfev={result.function_evaluations} {result.message}",
        flush=True,
    )
    return 0 if result.converged else 2


if __name__ == "__main__":
    raise SystemExit(main())
