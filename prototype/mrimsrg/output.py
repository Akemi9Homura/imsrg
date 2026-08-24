"""Self-describing materialized output for the rapid MR-IMSRG prototype."""

from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import json
from pathlib import Path
import struct

import numpy as np

try:
    from .densities import Densities
    from .flow import FlowResult, FlowSettings
    from .normal_order import MRHamiltonian, VacuumHamiltonian, to_vacuum
    from .reference_io import ReferenceData
except ImportError:
    from densities import Densities
    from flow import FlowResult, FlowSettings
    from normal_order import MRHamiltonian, VacuumHamiltonian, to_vacuum
    from reference_io import ReferenceData


_BRIDGE_MAGIC = b"mrimsrg_m_v1\0\0\0\0"


def _write_bridge_payload(
    path: Path, norb: int, hamiltonian: VacuumHamiltonian
) -> None:
    with path.open("wb") as stream:
        stream.write(_BRIDGE_MAGIC)
        stream.write(struct.pack("<Qd", norb, hamiltonian.zero_body))
        stream.write(
            np.asarray(hamiltonian.one_body, dtype="<f8").tobytes(order="C")
        )
        stream.write(
            np.asarray(hamiltonian.two_body, dtype="<f8").tobytes(order="C")
        )


def _checkpoint_name(s: float) -> str:
    return "s" + format(s, ".10g").replace("-", "m").replace(".", "p")


def _write_vacuum_checkpoint(
    root: Path,
    s: float,
    point: dict[str, object],
    hamiltonian: VacuumHamiltonian,
    norb: int,
) -> dict[str, object]:
    name = _checkpoint_name(s)
    path = root / "checkpoints" / name
    path.mkdir(parents=True)
    np.save(path / "vacuum_one_body.npy", hamiltonian.one_body, allow_pickle=False)
    np.save(path / "vacuum_two_body.npy", hamiltonian.two_body, allow_pickle=False)
    _write_bridge_payload(path / "vacuum_mscheme.bin", norb, hamiltonian)
    metadata = {
        "schema": "mrimsrg_vacuum_checkpoint_v1",
        "s": s,
        "vacuum_zero_body": hamiltonian.zero_body,
        "flow_point": point,
        "bridge_payload": "vacuum_mscheme.bin",
    }
    with (path / "metadata.json").open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return {"s": s, "path": f"checkpoints/{name}", **metadata}


def save_flow_output(
    path: str | Path,
    reference_path: str | Path,
    reference: ReferenceData,
    densities: Densities,
    initial: MRHamiltonian,
    result: FlowResult,
    settings: FlowSettings,
) -> VacuumHamiltonian:
    """Save all inputs and both MR/vacuum final representations.

    The directory is created atomically enough for the intended single writer:
    an existing path is always rejected so a prior result cannot be silently
    overwritten.
    """
    root = Path(path)
    if root.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {root}")
    root.mkdir(parents=True)
    final_vacuum = to_vacuum(result.hamiltonian, densities)

    arrays = {
        "orbits": reference.orbits,
        "gamma1": densities.gamma1,
        "gamma2": densities.gamma2,
        "lambda2": densities.lambda2,
        "initial_mr_one_body": initial.one_body,
        "initial_mr_two_body": initial.two_body,
        "final_mr_one_body": result.hamiltonian.one_body,
        "final_mr_two_body": result.hamiltonian.two_body,
        "final_vacuum_one_body": final_vacuum.one_body,
        "final_vacuum_two_body": final_vacuum.two_body,
    }
    for name, values in arrays.items():
        np.save(root / f"{name}.npy", values, allow_pickle=False)

    # Compact bridge payload: 16-byte magic, uint64 norb, float64 E0, then
    # row-major float64 one- and two-body arrays.  The NPY files remain the
    # human-inspectable canonical output; this file only avoids embedding a
    # general NPY parser in the validated C++ NCSM bridge.
    _write_bridge_payload(root / "vacuum_mscheme.bin", reference.norb, final_vacuum)

    initial_vacuum = to_vacuum(initial, densities)
    saved_checkpoints = [
        _write_vacuum_checkpoint(
            root,
            0.0,
            asdict(result.trajectory[0]),
            initial_vacuum,
            reference.norb,
        )
    ]
    for checkpoint in result.checkpoints:
        saved_checkpoints.append(
            _write_vacuum_checkpoint(
                root,
                checkpoint.point.s,
                asdict(checkpoint.point),
                to_vacuum(checkpoint.hamiltonian, densities),
                reference.norb,
            )
        )

    metadata = {
        "schema": "mrimsrg_flow_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "reference_path": str(Path(reference_path).resolve()),
        "reference_metadata": reference.metadata,
        "generator": "brillouin_delta_e_masked",
        "commutator": "MR-IMSRG(2), lambda3=0",
        "density_approximation": "lambda3=0",
        "decoupling_mask": {
            "one_body": "2*n(p)+l(p) != 2*n(q)+l(q)",
            "two_body": "e(p)+e(q) != e(r)+e(s)",
        },
        "ode_method": "DOP853 direct flow",
        "flow_settings": asdict(settings),
        "flow_converged": result.converged,
        "flow_message": result.message,
        "function_evaluations": result.function_evaluations,
        "initial_mr_zero_body": initial.zero_body,
        "initial_vacuum_zero_body": initial_vacuum.zero_body,
        "final_mr_zero_body": result.hamiltonian.zero_body,
        "final_vacuum_zero_body": final_vacuum.zero_body,
        "bridge_payload": "vacuum_mscheme.bin",
        "bridge_payload_layout": "magic[16], uint64 norb, float64 E0, float64 t[norb,norb], float64 V[norb,norb,norb,norb], little-endian C-order",
        "one_body_convention": "t[p,q] a^dagger_p a_q",
        "two_body_convention": "(1/4) V[p,q,r,s] a^dagger_p a^dagger_q a_s a_r",
        "trajectory": [asdict(point) for point in result.trajectory],
        "vacuum_checkpoints": saved_checkpoints,
        "final_vacuum_location": ".",
    }
    with (root / "metadata.json").open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return final_vacuum
