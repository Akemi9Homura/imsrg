#!/usr/bin/env python3
"""Unit checks for the MR-Magnus gate report parsers."""

from pathlib import Path
import sys
import tempfile


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))

from prototype.mrimsrg.summarize_magnus_gate import (  # noqa: E402
    parse_energies,
    parse_flow,
    parse_resource_usage,
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

    print("MR Magnus gate parser regression passed")


if __name__ == "__main__":
    main()
