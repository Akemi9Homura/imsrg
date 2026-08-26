#!/usr/bin/env python3
"""Unit checks for the MR-Magnus gate report parsers."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))

from prototype.mrimsrg.summarize_magnus_gate import (  # noqa: E402
    omega_segments_are_materialized,
    parse_energies,
    parse_flow,
    profile_metadata_gates,
    parse_resource_usage,
)
from prototype.mrimsrg.aggregate_magnus_gates import (  # noqa: E402
    EXPECTED_BARE_J64_SHA256,
    EXPECTED_INTERACTION_SHA256,
    EXPECTED_NREFMAX,
    EXPECTED_NMAX,
    EXPECTED_PRODUCTION_EXECUTABLE_SHA256,
    EXPECTED_PRODUCTION_LIBRARY_SHA256,
    EXPECTED_REFERENCE_SHA256,
    EXPECTED_SINGLE_GATES,
    EXPECTED_SINGLE_GATE_SCHEMA,
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

        default_profile = {
            "metadata": {
                "settings": {
                    "nucleus": "He4", "nrefmax": 2, "emax": 2,
                    "start_s": 0.0, "target_s": 1e7,
                    "eta_criterion": 7e-7, "tight_validation": False,
                },
                "solver_method": "magnus_adaptive",
                "ode_parameter_source": "imsrg++ adaptive Magnus runtime defaults",
                "ode_parameter_overrides": [],
                "validation_profile": "production_default",
            },
            "flow": {"final": {"s": 12.5}},
        }
        tight_profile = {
            "metadata": {
                "settings": {
                    "nucleus": "He4", "nrefmax": 2, "emax": 2,
                    "start_s": 0.0, "target_s": 12.5,
                    "eta_criterion": 1e-30, "tight_validation": True,
                },
                "solver_method": "magnus_adaptive",
                "ode_parameter_source": (
                    "explicit tenfold validation override of the imsrg++ default"
                ),
                "ode_parameter_overrides": ["ode_tolerance=1e-7"],
                "validation_profile": "tight_ode",
                "cumulative_target_s": 12.5,
                "omega_segment_target_s": 12.5,
            },
            "flow": {"final": {"s": 12.5}},
        }
        metadata_gates = profile_metadata_gates(
            default_profile, tight_profile,
            nucleus="He4", nrefmax=2, endpoint_tolerance=1e-4,
        )
        require(all(metadata_gates.values()), "valid flow profiles were rejected")
        broken_tight_profile = {
            **tight_profile,
            "metadata": {
                **tight_profile["metadata"],
                "ode_parameter_overrides": [
                    "ode_tolerance=1e-7", "dsmax=0.05"
                ],
            },
        }
        broken_gates = profile_metadata_gates(
            default_profile, broken_tight_profile,
            nucleus="He4", nrefmax=2, endpoint_tolerance=1e-4,
        )
        require(
            not broken_gates["tight_uses_only_tenfold_ode_tolerance"],
            "a tight profile with an extra step-size override was accepted",
        )
        broken_tight_profile = {
            **tight_profile,
            "metadata": {
                **tight_profile["metadata"],
                "settings": {
                    **tight_profile["metadata"]["settings"],
                    "target_s": 12.6,
                },
                "cumulative_target_s": 12.6,
                "omega_segment_target_s": 12.6,
            },
        }
        broken_gates = profile_metadata_gates(
            default_profile, broken_tight_profile,
            nucleus="He4", nrefmax=2, endpoint_tolerance=1e-4,
        )
        require(
            not broken_gates["tight_configured_for_default_endpoint"],
            "a tight profile configured for a different endpoint was accepted",
        )

        entries = []
        for index, (nucleus, nrefmax) in enumerate(EXPECTED_NREFMAX.items()):
            profile = {
                "metadata": {
                    "executable_sha256": EXPECTED_PRODUCTION_EXECUTABLE_SHA256,
                    "shared_library_sha256": EXPECTED_PRODUCTION_LIBRARY_SHA256,
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
                "schema": EXPECTED_SINGLE_GATE_SCHEMA,
                "nucleus": nucleus,
                "nrefmax": nrefmax,
                "nmax": EXPECTED_NMAX[nucleus],
                "states": 3,
                "interaction_sha256": EXPECTED_INTERACTION_SHA256,
                "downstream_validator": {"sha256": "validator-sha"},
                "thresholds": EXPECTED_THRESHOLDS,
                "spectral_max_abs_mev": (index + 1) * 1e-4,
                "gates": {gate: True for gate in EXPECTED_SINGLE_GATES},
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
        combined_path = root / "four_nucleus_gate.json"
        subprocess.run(
            [
                sys.executable,
                str(REPOSITORY / "prototype/mrimsrg/aggregate_magnus_gates.py"),
                *(str(path) for path, _ in entries),
                "--json", str(combined_path),
            ],
            check=True,
            text=True,
            capture_output=True,
        )
        require(
            json.loads(combined_path.read_text(encoding="utf-8"))["passed"],
            "four-nucleus gate CLI did not write a passing report",
        )
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

        stale_entries = list(entries)
        stale_report = dict(stale_entries[0][1])
        stale_report["schema"] = "mrimsrg_magnus_gate_v1"
        stale_entries[0] = (stale_entries[0][0], stale_report)
        try:
            aggregate(stale_entries)
        except ValueError:
            pass
        else:
            raise AssertionError("a stale single-nucleus gate schema was accepted")

        unfrozen_entries = []
        for path, report in entries:
            unfrozen_report = dict(report)
            for profile_name in ("default", "tight"):
                unfrozen_profile = dict(unfrozen_report[profile_name])
                unfrozen_metadata = dict(unfrozen_profile["metadata"])
                unfrozen_metadata["executable_sha256"] = "consistent-but-unfrozen"
                unfrozen_profile["metadata"] = unfrozen_metadata
                unfrozen_report[profile_name] = unfrozen_profile
            unfrozen_entries.append((path, unfrozen_report))
        require(
            not aggregate(unfrozen_entries)["passed"],
            "a consistent but unfrozen production executable was accepted",
        )

        unfrozen_entries = []
        for path, report in entries:
            unfrozen_report = dict(report)
            for profile_name in ("default", "tight"):
                unfrozen_profile = dict(unfrozen_report[profile_name])
                unfrozen_metadata = dict(unfrozen_profile["metadata"])
                unfrozen_metadata["shared_library_sha256"] = (
                    "consistent-but-unfrozen"
                )
                unfrozen_profile["metadata"] = unfrozen_metadata
                unfrozen_report[profile_name] = unfrozen_profile
            unfrozen_entries.append((path, unfrozen_report))
        require(
            not aggregate(unfrozen_entries)["passed"],
            "a consistent but unfrozen production library was accepted",
        )

    print("MR Magnus gate parser regression passed")


if __name__ == "__main__":
    main()
