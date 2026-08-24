"""Run and materialize one fixed-input MR-IMSRG calculation."""

from __future__ import annotations

import argparse
from pathlib import Path

try:
    from .densities import compute_densities, reference_energy, validate_densities
    from .flow import FlowPoint, FlowSettings, integrate_flow
    from .generator import oscillator_quanta_from_orbits
    from .normal_order import VacuumHamiltonian, normal_order, validate_hermitian
    from .output import save_flow_output
    from .reference_io import load_reference
except ImportError:
    from densities import compute_densities, reference_energy, validate_densities
    from flow import FlowPoint, FlowSettings, integrate_flow
    from generator import oscillator_quanta_from_orbits
    from normal_order import VacuumHamiltonian, normal_order, validate_hermitian
    from output import save_flow_output
    from reference_io import load_reference


def _print_point(point: FlowPoint) -> None:
    print(
        f"step={point.step:3d} s={point.s:.10g} E={point.zero_body:.12f} "
        f"residual={point.residual:.8e} ratio={point.residual_ratio:.8e} "
        f"generator_ratio={point.generator_numerator_residual_ratio:.8e} "
        f"sym={max(point.one_body_hermiticity_error, point.two_body_hermiticity_error, point.two_body_antisymmetry_error):.2e}",
        flush=True,
    )


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
    )
    final_vacuum = save_flow_output(
        args.output, args.reference, reference, densities, initial, result, settings
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
