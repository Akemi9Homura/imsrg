#!/usr/bin/env python3
"""Unit checks for the MR-Magnus gate report parsers."""

import json
from pathlib import Path
import sys
import tempfile


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))

from prototype.mrimsrg.summarize_magnus_gate import (  # noqa: E402
    omega_segments_are_materialized,
    parse_energies,
    parse_flow,
    parse_resource_usage,
)
from prototype.mrimsrg.aggregate_magnus_gates import (  # noqa: E402
    EXPECTED_BARE_J64_SHA256,
    EXPECTED_INTERACTION_SHA256,
    EXPECTED_NREFMAX,
    EXPECTED_NMAX,
    EXPECTED_REFERENCE_SHA256,
    EXPECTED_THRESHOLDS,
    aggregate,
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    with tempfile.TemporaryDirectory(prefix="mr-magnus-gate-") as temporary:
        root = Path(temporary)
        flow = root / "flow.dat"
        flow.write_text(
            "0 0 -10 20 0 0 0 0 3 4 0 0 0 0 0 0 0\n"
            "2 7 -11 19 0 0.1 0.2 0 3e-6 4e-6 0 42 0 0 0 0 0\n",
            encoding="utf-8",
        )
        parsed = parse_flow(flow)
        require(parsed["initial"]["eta_norm"] == 5.0,
                "initial eta norm was parsed incorrectly")
        require(abs(parsed["final"]["eta_norm"] - 5e-6) < 1e-18,
                "final eta norm was parsed incorrectly")
        require(abs(parsed["final"]["eta_ratio_to_initial"] - 1e-6) < 1e-18,
                "eta ratio was parsed incorrectly")
        require(parsed["final"]["scalar_commutators"] == 42,
                "commutator count was parsed incorrectly")

        resources = root / "resource_usage.txt"
        resources.write_text(
            "wall_seconds=73\nmaximum_rss_kib=14152\nexit_status=0\n",
            encoding="utf-8",
        )
        usage = parse_resource_usage(resources)
        require(usage["wall_seconds"] == 73 and usage["exit_status"] == 0,
                "resource usage was parsed incorrectly")

        energies = parse_energies(
            "state=2 E=-1.5 Ex=2 twoJ=0\n"
            "state=0 E=-3.5 Ex=0 twoJ=0\n"
            "state=1 E=-2.0 Ex=1.5 twoJ=2\n"
        )
        require(energies == [-3.5, -2.0, -1.5],
                "NCSM energies were not sorted by state")

        require(
            omega_segments_are_materialized(
                {"omega_files": 20, "empty_omega_files": 0},
                {"omega_files": 21, "empty_omega_files": 0},
            ),
            "nonempty Omega segments were rejected",
        )
        require(
            not omega_segments_are_materialized(
                {"omega_files": 0, "empty_omega_files": 0},
                {"omega_files": 1, "empty_omega_files": 0},
            ),
            "missing Omega segments were accepted",
        )
        require(
            not omega_segments_are_materialized(
                {"omega_files": 1, "empty_omega_files": 1},
                {"omega_files": 1, "empty_omega_files": 0},
            ),
            "an empty Omega segment was accepted",
        )

        entries = []
        for index, (nucleus, nrefmax) in enumerate(EXPECTED_NREFMAX.items()):
            profile = {
                "metadata": {
                    "executable_sha256": "exe-sha",
                    "shared_library_sha256": "library-sha",
                    "reference_sha256": EXPECTED_REFERENCE_SHA256[nucleus],
                    "interaction_sha256": EXPECTED_BARE_J64_SHA256[nucleus],
                },
                "log": f"/result/log_{nucleus}_{1000 + index}.txt",
                "flow": {"final": {"s": 10.0 + index,
                                     "eta_ratio_to_initial": 9e-7}},
                "packing_max_abs_mev": 2e-6,
                "omega_files": index + 1,
                "omega_validation": {
                    "executable_sha256": "omega-validator-sha",
                    "zero_body_max_abs": 1e-22,
                    "one_body_antihermiticity_max_abs": 0.0,
                    "two_body_antihermiticity_max_abs": 0.0,
                },
                "magnus_series_rejections": 0,
            }
            report = {
                "schema": "mrimsrg_magnus_gate_v1",
                "nucleus": nucleus,
                "nrefmax": nrefmax,
                "nmax": EXPECTED_NMAX[nucleus],
                "states": 3,
                "interaction_sha256": EXPECTED_INTERACTION_SHA256,
                "downstream_validator": {"sha256": "validator-sha"},
                "thresholds": EXPECTED_THRESHOLDS,
                "spectral_max_abs_mev": (index + 1) * 1e-4,
                "default": dict(profile),
                "tight": dict(profile),
                "passed": True,
            }
            path = root / f"{nucleus}.json"
            path.write_text(json.dumps(report), encoding="utf-8")
            entries.append((path, report))

        combined = aggregate(entries)
        require(combined["passed"], "consistent four-nucleus gates were rejected")
        require(abs(combined["max_spectral_abs_keV"] - 0.4) < 1e-15,
                "four-nucleus spectral maximum was aggregated incorrectly")
        broken_entries = list(entries)
        broken_report = dict(broken_entries[0][1])
        broken_default = dict(broken_report["default"])
        broken_metadata = dict(broken_default["metadata"])
        broken_metadata["executable_sha256"] = "different-executable"
        broken_default["metadata"] = broken_metadata
        broken_report["default"] = broken_default
        broken_entries[0] = (broken_entries[0][0], broken_report)
        try:
            aggregate(broken_entries)
        except ValueError:
            pass
        else:
            raise AssertionError("inconsistent production executables were accepted")

    print("MR Magnus gate parser regression passed")


if __name__ == "__main__":
    main()
