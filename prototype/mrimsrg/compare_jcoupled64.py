#!/usr/bin/env python3
"""Compare two lossless scalar J-coupled Hamiltonian checkpoints."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

try:
    from .sr_imsrgpp_check import read_jcoupled64
except ImportError:
    from sr_imsrgpp_check import read_jcoupled64


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def compare(first_path: Path, second_path: Path) -> dict:
    first = read_jcoupled64(first_path)
    second = read_jcoupled64(second_path)
    if first.hw != second.hw or first.emax != second.emax:
        raise ValueError("jcoupled64 files use different hw/emax")
    if not np.array_equal(first.orbits, second.orbits):
        raise ValueError("jcoupled64 files use different orbit tables")
    first_keys = [record[:5] for record in first.records]
    second_keys = [record[:5] for record in second.records]
    if first_keys != second_keys:
        raise ValueError("jcoupled64 files use different two-body record order")

    one_difference = first.one_body - second.one_body
    first_two = np.asarray([record[5] for record in first.records])
    second_two = np.asarray([record[5] for record in second.records])
    two_difference = first_two - second_two
    worst_one_flat = int(np.argmax(np.abs(one_difference)))
    worst_one = np.unravel_index(worst_one_flat, one_difference.shape)
    worst_two = int(np.argmax(np.abs(two_difference)))
    return {
        "schema": "mrimsrg_jcoupled64_comparison_v1",
        "first": str(first_path.resolve()),
        "second": str(second_path.resolve()),
        "first_sha256": sha256(first_path),
        "second_sha256": sha256(second_path),
        "hw": first.hw,
        "emax": first.emax,
        "norbits": int(len(first.orbits)),
        "ntbme": int(len(first.records)),
        "zero_body_max_abs_mev": abs(first.zero_body - second.zero_body),
        "one_body_max_abs_mev": float(np.max(np.abs(one_difference))),
        "one_body_relative_frobenius": float(
            np.linalg.norm(one_difference) / max(np.linalg.norm(first.one_body), 1e-300)
        ),
        "one_body_worst_indices": [int(value) for value in worst_one],
        "two_body_max_abs_mev": float(np.max(np.abs(two_difference))),
        "two_body_relative_frobenius": float(
            np.linalg.norm(two_difference) / max(np.linalg.norm(first_two), 1e-300)
        ),
        "two_body_worst_record": [int(value) for value in first_keys[worst_two]],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("first", type=Path)
    parser.add_argument("second", type=Path)
    parser.add_argument("--absolute-tolerance", type=float, default=1e-10)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    report = compare(args.first, args.second)
    maximum = max(
        report["zero_body_max_abs_mev"],
        report["one_body_max_abs_mev"],
        report["two_body_max_abs_mev"],
    )
    report["absolute_tolerance_mev"] = args.absolute_tolerance
    report["passed"] = maximum <= args.absolute_tolerance
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if args.json is not None:
        args.json.write_text(rendered, encoding="utf-8")
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
