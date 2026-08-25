#!/usr/bin/env python3

"""Embed a validated fixed-Nrefmax MR reference in a larger HO model space."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import struct

try:
    from .export_jref import ENDIAN_MARKER, MAGIC, SCHEMA_VERSION
    from .sr_imsrgpp_check import import_pyimsrg
except ImportError:
    from export_jref import ENDIAN_MARKER, MAGIC, SCHEMA_VERSION
    from sr_imsrgpp_check import import_pyimsrg


_JREF_HEADER = struct.Struct("<16sII7idQQ")
_MINIPACK_HEADER = struct.Struct("<8sdiI")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_jref_metadata(path: Path) -> dict[str, int | float]:
    with path.open("rb") as stream:
        payload = stream.read(_JREF_HEADER.size)
    if len(payload) != _JREF_HEADER.size:
        raise ValueError(f"truncated MR reference header: {path}")
    (
        magic,
        endian,
        version,
        mass,
        charge,
        nrefmax,
        j2,
        parity,
        emax,
        e2max,
        hw,
        norbits,
        nchannels,
    ) = _JREF_HEADER.unpack(payload)
    if magic != MAGIC or endian != ENDIAN_MARKER or version != SCHEMA_VERSION:
        raise ValueError(f"unsupported MR reference header: {path}")
    return {
        "A": mass,
        "Z": charge,
        "Nrefmax": nrefmax,
        "J2": j2,
        "parity": parity,
        "emax": emax,
        "e2max": e2max,
        "hw": hw,
        "norbits": norbits,
        "nchannels": nchannels,
    }


def _read_minipack_metadata(path: Path) -> dict[str, int | float]:
    with path.open("rb") as stream:
        payload = stream.read(_MINIPACK_HEADER.size)
    if len(payload) != _MINIPACK_HEADER.size:
        raise ValueError(f"truncated minipack header: {path}")
    magic, hw, emax, tbme_count = _MINIPACK_HEADER.unpack(payload)
    if magic != b"minipack":
        raise ValueError(f"unsupported interaction header: {path}")
    if emax < 0 or tbme_count < 0:
        raise ValueError(f"invalid minipack dimensions: {path}")
    return {"hw": hw, "emax": emax, "e2max": 2 * emax, "tbme_count": tbme_count}


def embed_reference(
    source: Path,
    interaction: Path,
    output: Path,
    pyimsrg_dir: Path,
    *,
    tolerance: float = 1e-10,
) -> dict[str, object]:
    """Embed, write, and independently read back one larger-space reference."""

    source_metadata = _read_jref_metadata(source)
    interaction_metadata = _read_minipack_metadata(interaction)
    if interaction_metadata["emax"] < source_metadata["emax"]:
        raise ValueError("target minipack emax is smaller than the source reference")
    if interaction_metadata["e2max"] < source_metadata["e2max"]:
        raise ValueError("target minipack e2max is smaller than the source reference")
    if abs(float(interaction_metadata["hw"]) - float(source_metadata["hw"])) > tolerance:
        raise ValueError("target minipack hw differs from the source reference")
    if output.exists():
        raise FileExistsError(f"refusing to overwrite embedded MR reference: {output}")

    pyimsrg = import_pyimsrg(pyimsrg_dir)
    source_modelspace = pyimsrg.ModelSpace(int(source_metadata["emax"]), "He4", "He4")
    source_modelspace.SetE2max(int(source_metadata["e2max"]))
    source_modelspace.SetHbarOmega(float(source_metadata["hw"]))
    source_modelspace.SetReferenceOcc(
        pyimsrg.MRReference.ReadOccupationMap(source_modelspace, str(source))
    )
    source_reference = pyimsrg.MRReference.ReadBinary(
        source_modelspace, str(source), tolerance
    )

    target_modelspace = pyimsrg.ModelSpace(
        int(interaction_metadata["emax"]), "He4", "He4"
    )
    target_modelspace.SetE2max(int(interaction_metadata["e2max"]))
    target_modelspace.SetHbarOmega(float(interaction_metadata["hw"]))
    interaction_sha256 = _sha256(interaction)
    embedded = source_reference.EmbedInModelSpace(
        target_modelspace, interaction_sha256, tolerance
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    embedded.WriteBinary(str(output))

    roundtrip_modelspace = pyimsrg.ModelSpace(
        int(interaction_metadata["emax"]), "He4", "He4"
    )
    roundtrip_modelspace.SetE2max(int(interaction_metadata["e2max"]))
    roundtrip_modelspace.SetHbarOmega(float(interaction_metadata["hw"]))
    roundtrip_modelspace.SetReferenceOcc(
        pyimsrg.MRReference.ReadOccupationMap(roundtrip_modelspace, str(output))
    )
    roundtrip = pyimsrg.MRReference.ReadBinary(
        roundtrip_modelspace, str(output), tolerance
    )
    return {
        "schema": "mrimsrg_embedded_reference_v1",
        "source": str(source.resolve()),
        "output": str(output.resolve()),
        "A": roundtrip.A,
        "Z": roundtrip.Z,
        "Nrefmax": roundtrip.Nrefmax,
        "source_emax": source_metadata["emax"],
        "target_emax": roundtrip.emax,
        "target_e2max": roundtrip.e2max,
        "hw": roundtrip.hw,
        "target_interaction": str(interaction.resolve()),
        "target_interaction_sha256": roundtrip.interaction_sha256,
        "source_rdm_sha256": roundtrip.rdm_sha256,
        "source_wavefunction_sha256": roundtrip.wavefunction_sha256,
        "maximum_hermiticity_violation": roundtrip.MaximumHermiticityViolation(),
        "maximum_contraction_violation": roundtrip.MaximumContractionViolation(),
        "embedding_definition": (
            "source blocks unchanged; added occupations/lambda2 zero; "
            "added natural-orbit block identity"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--interaction", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--pyimsrg-dir", required=True, type=Path)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    if args.tolerance <= 0.0:
        parser.error("--tolerance must be positive")
    report = embed_reference(
        args.source,
        args.interaction,
        args.output,
        args.pyimsrg_dir,
        tolerance=args.tolerance,
    )
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if args.json is not None:
        if args.json.exists():
            raise FileExistsError(f"refusing to overwrite report: {args.json}")
        args.json.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
