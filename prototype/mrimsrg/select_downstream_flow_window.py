#!/usr/bin/env python3

"""Apply the prototype finite-flow downstream acceptance policy.

This is deliberately separate from the strict White-NCSM generator test.
Vobig Sec. 6.5.5 diagnoses useful or unstable IM-NCSM flow regions from the
post-diagonalized Nmax sequence.  A point accepted here is therefore only a
finite-flow convergence accelerator; it is not a converged MR-IMSRG result.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def flow_parameter(label: str) -> float:
    if not label.startswith("s"):
        raise ValueError(f"invalid flow-point label: {label}")
    return float(label[1:])


def evaluate_window(window: dict[str, object]) -> dict[str, object]:
    policy = window["finite_s_operational_gate"]
    thresholds = policy["thresholds"]
    energies = window["ground_state_mev"]
    bare = energies["s0"]
    full_key = f"Nmax{window['full_space_nmax']}"
    truncated_keys = sorted(
        (key for key in bare if key != full_key),
        key=lambda key: int(key.removeprefix("Nmax")),
    )

    points: list[dict[str, object]] = []
    for label, curve in energies.items():
        if label == "s0":
            continue
        drift_kev = 1000.0 * (curve[full_key] - bare[full_key])
        improvements = {
            key: 1000.0
            * (
                abs(bare[key] - bare[full_key])
                - abs(curve[key] - curve[full_key])
            )
            for key in truncated_keys
        }
        physics_passed = (
            abs(drift_kev) <= thresholds["maximum_full_space_drift_kev"]
            and all(
                value >= thresholds["minimum_each_nmax_improvement_kev"]
                for value in improvements.values()
            )
        )

        diagnostics = policy["point_diagnostics"].get(label)
        numerical_passed = False
        residual_decreased = False
        if diagnostics is not None:
            numerical_passed = (
                diagnostics["maximum_packing_error_kev"]
                <= thresholds["maximum_packing_error_kev"]
                and diagnostics["maximum_ode_error_kev"]
                <= thresholds["maximum_ode_error_kev"]
            )
            residual_decreased = (
                diagnostics["eta_norm_ratio_to_initial"]
                < thresholds["maximum_eta_norm_ratio_to_initial"]
            )

        points.append(
            {
                "label": label,
                "s": flow_parameter(label),
                "full_space_drift_kev": drift_kev,
                "nmax_truncation_improvement_kev": improvements,
                "physics_passed": physics_passed,
                "numerical_passed": numerical_passed,
                "residual_decreased": residual_decreased,
                "passed": physics_passed
                and numerical_passed
                and residual_decreased,
            }
        )

    accepted = [point for point in points if point["passed"]]
    selected = max(accepted, key=lambda point: point["s"]) if accepted else None
    strict = window["strict_decoupling_gate"]
    strict_passed = (
        strict["best_eta_norm_ratio_to_initial"]
        <= strict["maximum_eta_norm_ratio_to_initial"]
    )
    report = {
        "policy": policy["name"],
        "policy_is_literature_threshold": False,
        "points": points,
        "selected_candidate": None if selected is None else selected["label"],
        "finite_s_gate_passed": selected is not None,
        "strict_decoupling_gate_passed": strict_passed,
        "claims_are_independent": True,
    }
    expected = policy.get("expected_selected_candidate")
    if expected is not None and report["selected_candidate"] != expected:
        raise ValueError(
            "stored finite-s decision is stale: expected "
            f"{expected}, recomputed {report['selected_candidate']}"
        )
    if strict_passed and selected is not None:
        # This is allowed in principle, but make sure malformed NaNs never pass.
        if not math.isfinite(strict["best_eta_norm_ratio_to_initial"]):
            raise ValueError("non-finite strict residual")
    return report


def main() -> None:
    repository = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=repository / "docs/MR-IMSRG-Jscheme-large-space.json",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    document = json.loads(args.input.read_text(encoding="utf-8"))
    report = evaluate_window(document["emax4_downstream_flow_window"])
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
