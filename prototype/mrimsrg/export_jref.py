"""Export a validated NCSM reference as compact J-scheme MR input.

The source remains the existing shell-model-obs/simpleFCI bridge.  This module
only performs the independently tested density construction, temporary natural
basis rotation and scalar J coupling before writing the production C++ reader
format.  No m-scheme density is needed by the C++ flow.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import struct

import numpy as np

try:
    from .basis import prepare_natural_basis
    from .densities import compute_densities, validate_densities
    from .jcoupling import couple_scalar_two_body, extract_j_orbits
    from .reference_io import load_reference
except ImportError:
    from basis import prepare_natural_basis
    from densities import compute_densities, validate_densities
    from jcoupling import couple_scalar_two_body, extract_j_orbits
    from reference_io import load_reference


MAGIC = b"mrimsrg_jref1\0\0\0"
ENDIAN_MARKER = 0x01020304
SCHEMA_VERSION = 1


def _digest_files(paths: tuple[Path, ...]) -> str:
    digest = hashlib.sha256()
    for path in paths:
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
    return digest.hexdigest()


def _digest_arrays(arrays: tuple[np.ndarray, ...]) -> str:
    digest = hashlib.sha256()
    for values in arrays:
        digest.update(np.asarray(values, dtype="<f8").tobytes(order="C"))
    return digest.hexdigest()


def _extract_spherical_transformation(
    vectors: np.ndarray, orbits: np.ndarray, tolerance: float
) -> np.ndarray:
    """Collapse a rotationally scalar m-scheme transformation to J-orbits."""

    j_orbits = extract_j_orbits(orbits)
    result = np.zeros((len(j_orbits), len(j_orbits)), dtype=np.float64)
    for p in range(orbits.shape[0]):
        for q in range(orbits.shape[0]):
            same_channel_and_m = (
                orbits[p, 2] == orbits[q, 2]
                and orbits[p, 3] == orbits[q, 3]
                and orbits[p, 4] == orbits[q, 4]
                and orbits[p, 5] == orbits[q, 5]
            )
            if not same_channel_and_m and abs(float(vectors[p, q])) > tolerance:
                raise ValueError("natural-orbit transformation is not rotationally scalar")
    for a in j_orbits:
        for b in j_orbits:
            if (a.l, a.j2, a.tz2) != (b.l, b.j2, b.tz2):
                continue
            values = [vectors[p, q] for p, q in zip(a.substates, b.substates)]
            if np.max(np.abs(np.asarray(values) - values[0])) > tolerance:
                raise ValueError("natural-orbit transformation depends on magnetic substate")
            result[a.index, b.index] = float(values[0])
    if np.linalg.norm(result.T @ result - np.eye(len(j_orbits)), ord=np.inf) > tolerance:
        raise ValueError("collapsed natural-orbit transformation is not orthogonal")
    return result


def export_reference(
    reference_path: str | Path,
    output_path: str | Path,
    *,
    tolerance: float = 1e-10,
) -> dict[str, object]:
    """Build densities from the saved wavefunction and write mrimsrg_jref_v1."""

    root = Path(reference_path)
    output = Path(output_path)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing MR reference: {output}")
    reference = load_reference(root)
    densities = compute_densities(reference.determinants, reference.coefficients)
    validate_densities(densities, reference.metadata["A"], tolerance=tolerance)
    natural = prepare_natural_basis(densities, tolerance=tolerance)
    validate_densities(natural.densities, reference.metadata["A"], tolerance=5 * tolerance)
    j_orbits = extract_j_orbits(reference.orbits)
    transformation = _extract_spherical_transformation(
        natural.vectors, reference.orbits, tolerance
    )

    occupations = np.zeros(len(j_orbits), dtype=np.float64)
    for orbit in j_orbits:
        diagonal = np.diag(natural.densities.gamma1)[list(orbit.substates)]
        if np.max(np.abs(diagonal - diagonal[0])) > tolerance:
            raise ValueError("natural occupations depend on magnetic substate")
        occupations[orbit.index] = float(np.mean(diagonal))
    lambda_blocks = couple_scalar_two_body(
        natural.densities.lambda2, reference.orbits
    )

    interaction_sha256 = reference.metadata["interaction_sha256"]
    rdm_sha256 = _digest_arrays(
        (densities.gamma1, densities.gamma2, densities.lambda2)
    )
    wavefunction_sha256 = _digest_files(
        (root / "orbits.npy", root / "determinants.npy", root / "coefficients.npy")
    )
    for name, value in (
        ("interaction", interaction_sha256),
        ("RDM", rdm_sha256),
        ("wavefunction", wavefunction_sha256),
    ):
        if len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
            raise ValueError(f"invalid lowercase SHA-256 for {name}")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("xb") as stream:
        stream.write(MAGIC)
        stream.write(struct.pack("<II", ENDIAN_MARKER, SCHEMA_VERSION))
        metadata = reference.metadata
        stream.write(
            struct.pack(
                "<7idQQ",
                metadata["A"],
                metadata["Z"],
                metadata["Nrefmax"],
                metadata["J2"],
                metadata["parity"],
                metadata["emax"],
                metadata["e2max"],
                float(metadata["hw"]),
                len(j_orbits),
                len(lambda_blocks),
            )
        )
        stream.write(interaction_sha256.encode("ascii"))
        stream.write(rdm_sha256.encode("ascii"))
        stream.write(wavefunction_sha256.encode("ascii"))
        for orbit in j_orbits:
            stream.write(
                struct.pack(
                    "<5id",
                    orbit.index,
                    orbit.n,
                    orbit.l,
                    orbit.j2,
                    orbit.tz2,
                    occupations[orbit.index],
                )
            )
        stream.write(np.asarray(transformation, dtype="<f8").tobytes(order="C"))
        for block in lambda_blocks:
            stream.write(struct.pack("<3iIQ", block.J, block.parity, block.Tz, 0, len(block.pairs)))
            for a, b in block.pairs:
                stream.write(struct.pack("<II", a, b))
            stream.write(np.asarray(block.matrix, dtype="<f8").tobytes(order="C"))

    return {
        "schema": "mrimsrg_jref_v1",
        "path": str(output.resolve()),
        "A": reference.metadata["A"],
        "Z": reference.metadata["Z"],
        "Nrefmax": reference.metadata["Nrefmax"],
        "norbits": len(j_orbits),
        "nchannels": len(lambda_blocks),
        "interaction_sha256": interaction_sha256,
        "rdm_sha256": rdm_sha256,
        "wavefunction_sha256": wavefunction_sha256,
        "natural_basis_is_identity": natural.is_identity,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    args = parser.parse_args()
    summary = export_reference(args.reference, args.output, tolerance=args.tolerance)
    for key, value in summary.items():
        print(f"{key}={value}")


if __name__ == "__main__":
    main()
