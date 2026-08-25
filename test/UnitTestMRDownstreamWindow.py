#!/usr/bin/env python3

"""Check the finite-flow policy against the recorded emax4 spectra."""

from copy import deepcopy
import json
from pathlib import Path
import sys


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))

from prototype.mrimsrg.select_downstream_flow_window import evaluate_window  # noqa: E402


def main() -> None:
    document = json.loads(
        (REPOSITORY / "docs/MR-IMSRG-Jscheme-large-space.json").read_text(
            encoding="utf-8"
        )
    )
    window = document["emax4_downstream_flow_window"]
    report = evaluate_window(window)
    assert report["selected_candidate"] == "s0.02"
    assert report["finite_s_gate_passed"]
    assert not report["strict_decoupling_gate_passed"]
    by_label = {point["label"]: point for point in report["points"]}
    assert by_label["s0.02"]["passed"]
    assert all(
        improvement > 0.0
        for improvement in by_label["s0.02"][
            "nmax_truncation_improvement_kev"
        ].values()
    )
    assert not by_label["s0.1"]["physics_passed"]

    too_strict = deepcopy(window)
    too_strict["finite_s_operational_gate"]["thresholds"][
        "maximum_full_space_drift_kev"
    ] = 40.0
    too_strict["finite_s_operational_gate"].pop("expected_selected_candidate")
    rejected = evaluate_window(too_strict)
    assert rejected["selected_candidate"] is None
    assert not rejected["finite_s_gate_passed"]

    print("MR downstream finite-flow window test passed")


if __name__ == "__main__":
    main()
