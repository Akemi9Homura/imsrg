#!/usr/bin/env python3
"""Summarize one matched-endpoint MR-Magnus production/tight gate."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess

try:
    from .compare_jcoupled64 import compare as compare_jcoupled64
except ImportError:
    from compare_jcoupled64 import compare as compare_jcoupled64


ENERGY_PATTERN = re.compile(r"^state=(\d+) E=([-+0-9.eE]+)", re.MULTILINE)
GATE_SCHEMA = "mrimsrg_magnus_gate_v2"
DEFAULT_ODE_PARAMETER_SOURCE = "imsrg++ adaptive Magnus runtime defaults"
TIGHT_ODE_PARAMETER_SOURCE = (
    "explicit tenfold validation override of the imsrg++ default"
)
FLOW_PREFIX_WIDTHS = (5, 12, *(16 for _ in range(9)), 7)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def unique(root: Path, pattern: str) -> Path:
    matches = sorted(root.glob(pattern))
    if len(matches) != 1:
        raise ValueError(
            f"expected one {pattern!r} below {root}, found {len(matches)}"
        )
    return matches[0]


def parse_flow_prefix(line: str) -> list[float] | None:
    """Parse through Ncomm using the fixed widths emitted by IMSRGSolver."""
    fixed_width = sum(FLOW_PREFIX_WIDTHS)
    if len(line) >= fixed_width:
        fields: list[float] = []
        start = 0
        try:
            for width in FLOW_PREFIX_WIDTHS:
                fields.append(float(line[start:start + width]))
                start += width
            return fields
        except ValueError:
            pass

    tokens = line.split()
    if len(tokens) < len(FLOW_PREFIX_WIDTHS):
        return None
    try:
        return [float(value) for value in tokens[:len(FLOW_PREFIX_WIDTHS)]]
    except ValueError:
        return None


def parse_flow(path: Path) -> dict[str, object]:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = parse_flow_prefix(line)
        if fields is not None:
            rows.append(fields)
    if not rows:
        raise ValueError(f"flow file contains no status rows: {path}")

    def status(row: list[float]) -> dict[str, object]:
        eta_components = row[8:11]
        omega_components = row[5:8]
        return {
            "step": int(row[0]),
            "s": row[1],
            "zero_body_mev": row[2],
            "hamiltonian_norm": row[3],
            "omega_components": omega_components,
            "omega_norm": math.sqrt(sum(value * value for value in omega_components)),
            "eta_components": eta_components,
            "eta_norm": math.sqrt(sum(value * value for value in eta_components)),
            "scalar_commutators": int(row[11]),
        }

    initial = status(rows[0])
    final = status(rows[-1])
    final["eta_ratio_to_initial"] = (
        final["eta_norm"] / initial["eta_norm"]
    )
    return {
        "path": str(path.resolve()),
        "status_rows": len(rows),
        "initial": initial,
        "final": final,
    }


def parse_resource_usage(path: Path) -> dict[str, object]:
    result: dict[str, object] = {"path": str(path.resolve())}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        try:
            result[key] = int(value)
        except ValueError:
            result[key] = value
    return result


def parse_energies(output: str) -> list[float]:
    states = sorted(
        (int(index), float(energy))
        for index, energy in ENERGY_PATTERN.findall(output)
    )
    if not states or [index for index, _ in states] != list(range(len(states))):
        raise ValueError("NCSM output does not contain a contiguous state list")
    return [energy for _, energy in states]


def diagonalize(
    validator: Path,
    hamiltonian: Path,
    *,
    interaction: Path,
    proton_number: int,
    neutron_number: int,
    nmax: int,
    states: int,
) -> list[float]:
    if hamiltonian.suffix == ".jcoupled64":
        command = [
            str(validator), "--interaction", str(interaction),
            "--jcoupled64", str(hamiltonian),
        ]
    else:
        command = [str(validator), "--no2bpack", str(hamiltonian)]
    command.extend([
        "--Z", str(proton_number), "--N", str(neutron_number),
        "--nmax", str(nmax), "--states", str(states),
    ])
    completed = subprocess.run(
        command, check=True, text=True, capture_output=True
    )
    energies = parse_energies(completed.stdout)
    if len(energies) != states:
        raise ValueError(
            f"NCSM returned {len(energies)} states instead of {states}"
        )
    return energies


def validate_omega_files(
    executable: Path,
    root: Path,
    *,
    nucleus: str,
    emax: int,
    hw: float,
    tolerance: float,
) -> dict[str, object]:
    omega_paths = sorted(root.glob("*_Omega_*"))
    if not omega_paths:
        return {
            "schema": "imsrg_scalar_omega_validation_v1",
            "files": 0,
            "passed": False,
            "reason": "no Omega segments were materialized",
            "executable": str(executable.resolve()),
            "executable_sha256": sha256(executable),
            "returncode": None,
        }
    command = [
        str(executable),
        "--emax", str(emax),
        "--hw", str(hw),
        "--reference", nucleus,
        "--tolerance", str(tolerance),
        *(str(path) for path in omega_paths),
    ]
    completed = subprocess.run(command, text=True, capture_output=True)
    json_lines = [
        line for line in completed.stdout.splitlines() if line.startswith("{")
    ]
    if not json_lines:
        raise RuntimeError(
            "Omega validator returned no JSON report:\n" + completed.stderr
        )
    report = json.loads(json_lines[-1])
    report["executable"] = str(executable.resolve())
    report["executable_sha256"] = sha256(executable)
    report["returncode"] = completed.returncode
    return report


def profile(root: Path) -> dict[str, object]:
    flow = parse_flow(unique(root, "flow_*.dat"))
    log_path = unique(root, "log_*.txt")
    log = log_path.read_text(encoding="utf-8")
    metadata_path = unique(root, "metadata.json")
    omega_paths = sorted(root.glob("*_Omega_*"))
    return {
        "root": str(root.resolve()),
        "metadata": json.loads(metadata_path.read_text(encoding="utf-8")),
        "flow": flow,
        "resource_usage": parse_resource_usage(unique(root, "resource_usage.txt")),
        "jcoupled64": unique(root, "H_*.jcoupled64"),
        "no2bpack": unique(root, "H_*.no2bpack"),
        "omega_files": len(omega_paths),
        "empty_omega_files": sum(path.stat().st_size == 0 for path in omega_paths),
        "magnus_series_rejections": log.count("Magnus series rejected"),
        "magnus_series_recovery_segments": log.count(
            "materialized the accepted segment"
        ),
        "log": str(log_path.resolve()),
    }


def omega_segments_are_materialized(*entries: dict[str, object]) -> bool:
    return all(
        entry["omega_files"] > 0 and entry["empty_omega_files"] == 0
        for entry in entries
    )


def profile_metadata_gates(
    default: dict[str, object],
    tight: dict[str, object],
    *,
    nucleus: str,
    nrefmax: int,
    endpoint_tolerance: float,
) -> dict[str, bool]:
    """Check that the numerical comparison really uses the frozen profiles."""
    default_metadata = default["metadata"]
    tight_metadata = tight["metadata"]
    default_settings = default_metadata["settings"]
    tight_settings = tight_metadata["settings"]

    def system_matches(metadata: dict[str, object]) -> bool:
        settings = metadata["settings"]
        return (
            settings.get("nucleus") == nucleus
            and settings.get("nrefmax") == nrefmax
            and settings.get("emax") == 2
            and settings.get("start_s") == 0.0
        )

    default_stop_s = default["flow"]["final"]["s"]
    tight_target_s = tight_settings.get("target_s")
    configured_endpoint_matches = (
        isinstance(tight_target_s, (int, float))
        and abs(tight_target_s - default_stop_s) <= endpoint_tolerance
        and tight_metadata.get("cumulative_target_s") == tight_target_s
        and tight_metadata.get("omega_segment_target_s") == tight_target_s
    )
    return {
        "metadata_matches_requested_system": (
            system_matches(default_metadata) and system_matches(tight_metadata)
        ),
        "solver_is_adaptive_magnus": (
            default_metadata.get("solver_method") == "magnus_adaptive"
            and tight_metadata.get("solver_method") == "magnus_adaptive"
        ),
        "default_uses_solver_ode_defaults": (
            default_metadata.get("ode_parameter_source")
            == DEFAULT_ODE_PARAMETER_SOURCE
            and default_metadata.get("ode_parameter_overrides") == []
            and default_settings.get("tight_validation") in (None, False)
            and default_metadata.get("validation_profile")
            in (None, "production_default")
        ),
        "tight_uses_only_tenfold_ode_tolerance": (
            tight_metadata.get("ode_parameter_source")
            == TIGHT_ODE_PARAMETER_SOURCE
            and tight_metadata.get("ode_parameter_overrides")
            == ["ode_tolerance=1e-7"]
            and tight_settings.get("tight_validation") is True
            and tight_metadata.get("validation_profile") == "tight_ode"
        ),
        "tight_configured_for_default_endpoint": configured_endpoint_matches,
        "tight_early_stop_disabled": (
            isinstance(tight_settings.get("eta_criterion"), (int, float))
            and tight_settings["eta_criterion"] <= 1e-20
        ),
    }


def summarize(args: argparse.Namespace) -> dict[str, object]:
    default = profile(args.default)
    tight = profile(args.tight)
    matrix = compare_jcoupled64(default["jcoupled64"], tight["jcoupled64"])
    for entry in (default, tight):
        entry["omega_validation"] = validate_omega_files(
            args.omega_validator,
            Path(entry["root"]),
            nucleus=args.nucleus,
            emax=matrix["emax"],
            hw=matrix["hw"],
            tolerance=args.omega_tolerance,
        )
        entry["jcoupled64_sha256"] = sha256(entry["jcoupled64"])
        entry["no2bpack_sha256"] = sha256(entry["no2bpack"])
        entry["jcoupled64_energies_mev"] = diagonalize(
            args.validator, entry["jcoupled64"], interaction=args.interaction,
            proton_number=args.proton_number,
            neutron_number=args.neutron_number,
            nmax=args.nmax, states=args.states,
        )
        entry["no2bpack_energies_mev"] = diagonalize(
            args.validator, entry["no2bpack"], interaction=args.interaction,
            proton_number=args.proton_number,
            neutron_number=args.neutron_number,
            nmax=args.nmax, states=args.states,
        )
        entry["packing_max_abs_mev"] = max(
            abs(left - right) for left, right in zip(
                entry["jcoupled64_energies_mev"],
                entry["no2bpack_energies_mev"],
            )
        )
        entry["jcoupled64"] = str(entry["jcoupled64"].resolve())
        entry["no2bpack"] = str(entry["no2bpack"].resolve())

    spectral_differences = [
        tight_energy - default_energy
        for default_energy, tight_energy in zip(
            default["jcoupled64_energies_mev"],
            tight["jcoupled64_energies_mev"],
        )
    ]
    maximum_spectral_difference = max(abs(value) for value in spectral_differences)
    endpoint_difference = abs(
        default["flow"]["final"]["s"] - tight["flow"]["final"]["s"]
    )
    gates = profile_metadata_gates(
        default,
        tight,
        nucleus=args.nucleus,
        nrefmax=args.nrefmax,
        endpoint_tolerance=args.endpoint_tolerance,
    )
    gates.update({
        "default_residual_ratio_below_limit": (
            default["flow"]["final"]["eta_ratio_to_initial"]
            < args.residual_ratio
        ),
        "matched_endpoint": endpoint_difference <= args.endpoint_tolerance,
        "spectral_stability": (
            maximum_spectral_difference < args.spectral_tolerance_mev
        ),
        "packing_readback": max(
            default["packing_max_abs_mev"], tight["packing_max_abs_mev"]
        ) < args.packing_tolerance_mev,
        "clean_exit": (
            default["resource_usage"].get("exit_status") == 0
            and tight["resource_usage"].get("exit_status") == 0
        ),
        "omega_segments_materialized": omega_segments_are_materialized(
            default, tight
        ),
        "omega_antihermiticity": (
            default["omega_validation"]["passed"]
            and tight["omega_validation"]["passed"]
            and default["omega_validation"]["returncode"] == 0
            and tight["omega_validation"]["returncode"] == 0
        ),
    })
    return {
        "schema": GATE_SCHEMA,
        "gate_implementation": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "nucleus": args.nucleus,
        "nrefmax": args.nrefmax,
        "nmax": args.nmax,
        "states": args.states,
        "interaction": str(args.interaction.resolve()),
        "interaction_sha256": sha256(args.interaction),
        "downstream_validator": {
            "path": str(args.validator.resolve()),
            "sha256": sha256(args.validator),
        },
        "default": default,
        "tight": tight,
        "matrix_difference": matrix,
        "endpoint_abs_difference": endpoint_difference,
        "tight_minus_default_spectrum_mev": spectral_differences,
        "spectral_max_abs_mev": maximum_spectral_difference,
        "thresholds": {
            "residual_ratio": args.residual_ratio,
            "endpoint_tolerance": args.endpoint_tolerance,
            "spectral_tolerance_mev": args.spectral_tolerance_mev,
            "packing_tolerance_mev": args.packing_tolerance_mev,
        },
        "gates": gates,
        "passed": all(gates.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nucleus", required=True)
    parser.add_argument("--nrefmax", required=True, type=int)
    parser.add_argument("--default", required=True, type=Path)
    parser.add_argument("--tight", required=True, type=Path)
    parser.add_argument("--interaction", required=True, type=Path)
    parser.add_argument("--validator", required=True, type=Path)
    parser.add_argument("--omega-validator", required=True, type=Path)
    parser.add_argument("--Z", dest="proton_number", required=True, type=int)
    parser.add_argument("--N", dest="neutron_number", required=True, type=int)
    parser.add_argument("--nmax", required=True, type=int)
    parser.add_argument("--states", type=int, default=3)
    parser.add_argument("--residual-ratio", type=float, default=1e-6)
    parser.add_argument("--endpoint-tolerance", type=float, default=1e-4)
    parser.add_argument("--spectral-tolerance-mev", type=float, default=1e-3)
    parser.add_argument("--packing-tolerance-mev", type=float, default=5e-6)
    parser.add_argument("--omega-tolerance", type=float, default=1e-10)
    parser.add_argument("--json", required=True, type=Path)
    args = parser.parse_args()
    report = summarize(args)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    args.json.write_text(rendered, encoding="utf-8")
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
