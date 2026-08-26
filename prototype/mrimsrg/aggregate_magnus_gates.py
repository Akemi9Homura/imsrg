#!/usr/bin/env python3
"""Aggregate the fixed four-nucleus matched-endpoint MR-Magnus gates."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from typing import Optional


EXPECTED_NREFMAX = {"He4": 2, "Be8": 0, "C12": 0, "O16": 2}
EXPECTED_NMAX = {"He4": 8, "Be8": 0, "C12": 0, "O16": 2}
EXPECTED_SINGLE_GATE_SCHEMA = "mrimsrg_magnus_gate_v2"
EXPECTED_SINGLE_GATES = {
    "clean_exit",
    "default_residual_ratio_below_limit",
    "default_uses_solver_ode_defaults",
    "matched_endpoint",
    "metadata_matches_requested_system",
    "omega_antihermiticity",
    "omega_segments_materialized",
    "packing_readback",
    "solver_is_adaptive_magnus",
    "spectral_stability",
    "tight_configured_for_default_endpoint",
    "tight_early_stop_disabled",
    "tight_uses_only_tenfold_ode_tolerance",
}
EXPECTED_INTERACTION_SHA256 = (
    "76b7243ef53d30955c0293d29da73688dc3839942143ccf147739108bb58ff84"
)
EXPECTED_PRODUCTION_EXECUTABLE_SHA256 = (
    "a70a28582d232d675c41a8db26a327b6d29bde8e5cf0317695d7d26745218de5"
)
EXPECTED_PRODUCTION_LIBRARY_SHA256 = (
    "c7eb668df585e96c29b0c7e659288ef16cbeecefcf358d9d390399488261927b"
)
EXPECTED_THRESHOLDS = {
    "endpoint_tolerance": 1e-4,
    "packing_tolerance_mev": 5e-6,
    "residual_ratio": 1e-6,
    "spectral_tolerance_mev": 1e-3,
}
EXPECTED_REFERENCE_SHA256 = {
    "He4": "1c9934aa92b50cc55122f15434b0093a860d99ae50a745d89848f76f1f608898",
    "Be8": "bbc3eb70c62e6022bb35d42e05582edd429d6031da1662ff19ab6b9c67416e17",
    "C12": "83b5d9ef6317a9290d37575f0af9819ca5f10babd62c62d95ed7045998192cda",
    "O16": "9c54eb88f9d57afdc09dde8da316ebaac818008fc67e05d092fe7d162907b21c",
}
EXPECTED_BARE_J64_SHA256 = {
    "He4": "55a3c161f78bde8586c997f332f89f0f1043a467af3c5ae90c1f98ed8fd344e2",
    "Be8": "75595b6b5e34f2f723f779497d81d1e0304c67e90bce9c8426d96a13e99249c0",
    "C12": "c669fcde829c3b12a0bcd77edcc60c4c9742e3c28455b5d9af4c306e2205bb46",
    "O16": "fe2027ffda381820a40f8617493152bfa0a017a6fc2e77a6f252b1758fac356e",
}
JOB_PATTERN = re.compile(r"_(\d+)\.txt$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def unique(values, description: str):
    distinct = set(values)
    if len(distinct) != 1:
        raise ValueError(f"four-nucleus gates use inconsistent {description}: {distinct}")
    return distinct.pop()


def job_id(log_path: str) -> Optional[int]:
    match = JOB_PATTERN.search(log_path)
    return None if match is None else int(match.group(1))


def aggregate(entries: list[tuple[Path, dict]]) -> dict[str, object]:
    by_nucleus: dict[str, tuple[Path, dict]] = {}
    for path, report in entries:
        if report.get("schema") != EXPECTED_SINGLE_GATE_SCHEMA:
            raise ValueError(f"unsupported gate schema in {path}")
        gates = report.get("gates")
        if not isinstance(gates, dict) or set(gates) != EXPECTED_SINGLE_GATES:
            raise ValueError(f"incomplete single-nucleus gates in {path}")
        if report.get("passed") is not True or not all(
            value is True for value in gates.values()
        ):
            raise ValueError(f"failed single-nucleus gate in {path}")
        nucleus = report.get("nucleus")
        if not isinstance(nucleus, str):
            raise ValueError(f"gate does not name a nucleus: {path}")
        if nucleus in by_nucleus:
            raise ValueError(f"duplicate gate for {nucleus}")
        by_nucleus[nucleus] = (path, report)
    if set(by_nucleus) != set(EXPECTED_NREFMAX):
        raise ValueError(
            f"expected gates for {sorted(EXPECTED_NREFMAX)}, got {sorted(by_nucleus)}"
        )

    for nucleus, expected_nrefmax in EXPECTED_NREFMAX.items():
        actual = by_nucleus[nucleus][1].get("nrefmax")
        if actual != expected_nrefmax:
            raise ValueError(
                f"{nucleus} gate uses Nrefmax={actual}, expected {expected_nrefmax}"
            )
        actual_nmax = by_nucleus[nucleus][1].get("nmax")
        if actual_nmax != EXPECTED_NMAX[nucleus]:
            raise ValueError(
                f"{nucleus} gate uses Nmax={actual_nmax}, expected {EXPECTED_NMAX[nucleus]}"
            )
        if by_nucleus[nucleus][1].get("states") != 3:
            raise ValueError(f"{nucleus} gate does not contain exactly three states")

    reports = [by_nucleus[nucleus][1] for nucleus in EXPECTED_NREFMAX]
    interaction_sha256 = unique(
        (report["interaction_sha256"] for report in reports), "interaction SHA-256"
    )
    downstream_validator_sha256 = unique(
        (report["downstream_validator"]["sha256"] for report in reports),
        "downstream validator SHA-256",
    )
    production_executable_sha256 = unique(
        (
            report[profile]["metadata"]["executable_sha256"]
            for report in reports
            for profile in ("default", "tight")
        ),
        "production executable SHA-256",
    )
    production_library_sha256 = unique(
        (
            report[profile]["metadata"]["shared_library_sha256"]
            for report in reports
            for profile in ("default", "tight")
        ),
        "production shared-library SHA-256",
    )
    omega_validator_sha256 = unique(
        (
            report[profile]["omega_validation"]["executable_sha256"]
            for report in reports
            for profile in ("default", "tight")
        ),
        "Omega validator SHA-256",
    )
    thresholds = unique(
        (json.dumps(report["thresholds"], sort_keys=True) for report in reports),
        "gate thresholds",
    )

    nuclei: dict[str, object] = {}
    for nucleus in EXPECTED_NREFMAX:
        path, report = by_nucleus[nucleus]
        omega_validations = [
            report[profile]["omega_validation"]
            for profile in ("default", "tight")
        ]
        nuclei[nucleus] = {
            "source_gate": str(path.resolve()),
            "source_gate_sha256": sha256(path),
            "nrefmax": report["nrefmax"],
            "nmax": report["nmax"],
            "reference_sha256": report["default"]["metadata"][
                "reference_sha256"
            ],
            "bare_jcoupled64_sha256": report["default"]["metadata"][
                "interaction_sha256"
            ],
            "default_job_id": job_id(report["default"]["log"]),
            "tight_job_id": job_id(report["tight"]["log"]),
            "default_stop_s": report["default"]["flow"]["final"]["s"],
            "default_residual_ratio": report["default"]["flow"]["final"][
                "eta_ratio_to_initial"
            ],
            "spectral_max_abs_keV": 1e3 * report["spectral_max_abs_mev"],
            "packing_max_abs_eV": 1e6 * max(
                report["default"]["packing_max_abs_mev"],
                report["tight"]["packing_max_abs_mev"],
            ),
            "default_omega_files": report["default"]["omega_files"],
            "tight_omega_files": report["tight"]["omega_files"],
            "omega_zero_body_max_abs": max(
                validation["zero_body_max_abs"] for validation in omega_validations
            ),
            "omega_one_body_antihermiticity_max_abs": max(
                validation["one_body_antihermiticity_max_abs"]
                for validation in omega_validations
            ),
            "omega_two_body_antihermiticity_max_abs": max(
                validation["two_body_antihermiticity_max_abs"]
                for validation in omega_validations
            ),
            "magnus_series_rejections": sum(
                report[profile]["magnus_series_rejections"]
                for profile in ("default", "tight")
            ),
            "passed": report["passed"],
        }

    consistency_gates = {
        "all_single_nucleus_gates_passed": all(
            report["passed"] and all(report["gates"].values())
            for report in reports
        ),
        "interaction_is_frozen": (
            interaction_sha256 == EXPECTED_INTERACTION_SHA256
        ),
        "downstream_validator_identical": bool(downstream_validator_sha256),
        "production_executable_is_frozen": (
            production_executable_sha256
            == EXPECTED_PRODUCTION_EXECUTABLE_SHA256
        ),
        "production_library_is_frozen": (
            production_library_sha256 == EXPECTED_PRODUCTION_LIBRARY_SHA256
        ),
        "omega_validator_identical": bool(omega_validator_sha256),
        "thresholds_are_frozen": json.loads(thresholds) == EXPECTED_THRESHOLDS,
        "input_artifacts_are_frozen": all(
            report[profile]["metadata"]["reference_sha256"]
            == EXPECTED_REFERENCE_SHA256[nucleus]
            and report[profile]["metadata"]["interaction_sha256"]
            == EXPECTED_BARE_J64_SHA256[nucleus]
            for nucleus, (_, report) in by_nucleus.items()
            for profile in ("default", "tight")
        ),
    }
    return {
        "schema": "mrimsrg_four_nucleus_magnus_gate_v1",
        "interaction_sha256": interaction_sha256,
        "downstream_validator_sha256": downstream_validator_sha256,
        "production_executable_sha256": production_executable_sha256,
        "production_library_sha256": production_library_sha256,
        "omega_validator_sha256": omega_validator_sha256,
        "thresholds": json.loads(thresholds),
        "nuclei": nuclei,
        "max_default_residual_ratio": max(
            entry["default_residual_ratio"] for entry in nuclei.values()
        ),
        "max_spectral_abs_keV": max(
            entry["spectral_max_abs_keV"] for entry in nuclei.values()
        ),
        "max_packing_abs_eV": max(
            entry["packing_max_abs_eV"] for entry in nuclei.values()
        ),
        "max_omega_antihermiticity_abs": max(
            max(
                entry["omega_one_body_antihermiticity_max_abs"],
                entry["omega_two_body_antihermiticity_max_abs"],
            )
            for entry in nuclei.values()
        ),
        "consistency_gates": consistency_gates,
        "passed": all(consistency_gates.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gates", nargs=4, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    args = parser.parse_args()
    entries = [(path, json.loads(path.read_text(encoding="utf-8"))) for path in args.gates]
    report = aggregate(entries)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    args.json.write_text(rendered, encoding="utf-8")
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
