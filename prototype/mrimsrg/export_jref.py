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
    from .jref_format import (ENDIAN_MARKER, MAGIC, REFERENCE_ENERGY_LABEL_BYTES,
                              SCHEMA_VERSION, SCHEMA_VERSION_LAMBDA3)
    from .jcoupling import couple_scalar_two_body, extract_j_orbits
    from .reference_io import load_reference
    from .lambda3 import coupled_lambda3
except ImportError:
    from basis import prepare_natural_basis
    from densities import compute_densities, validate_densities
    from jref_format import (ENDIAN_MARKER, MAGIC, REFERENCE_ENERGY_LABEL_BYTES,
                             SCHEMA_VERSION, SCHEMA_VERSION_LAMBDA3)
    from jcoupling import couple_scalar_two_body, extract_j_orbits
    from reference_io import load_reference
    from lambda3 import coupled_lambda3


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


def write_jref(
    output: Path,
    *,
    header: tuple,
    orbits: list[tuple[int, int, int, int, int]],
    occupations: np.ndarray,
    transformation: np.ndarray,
    lambda_blocks: list[tuple],
    digests: tuple[str, str, str],
    lambda3_blocks: list[tuple] | None = None,
    reference_energy: tuple[float, str] | None = None,
) -> None:
    """Write mrimsrg_jref (version 1, or 2 with lambda3/reference energy).

    header: (A, Z, Nrefmax, J2, parity, emax, e2max, hw); orbits: rows
    (index, n, l, 2j, 2tz) with contiguous indices; transformation[old, new]
    in that orbit order; lambda_blocks: (J, parity, Tz, pairs, matrix) with
    normalized fixed-M pair matrix elements; digests: interaction, RDM and
    wavefunction SHA-256 (lowercase hex).

    Version 2 appends: reference-energy flag, energy and a 64-byte ASCII label
    (which Hamiltonian produced it); lambda3 flag and channel count; per
    channel (2J, parity, 2Tz, reserved, nstates), the states (a, b, c, J12)
    with a <= b and the row-major matrix (lambda3.py conventions).  An empty
    lambda3_blocks list means "lambda3 present and zero"; None means absent.
    """
    for name, value in zip(("interaction", "RDM", "wavefunction"), digests):
        if len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
            raise ValueError(f"invalid lowercase SHA-256 for {name}")
    version = SCHEMA_VERSION if lambda3_blocks is None and reference_energy is None \
        else SCHEMA_VERSION_LAMBDA3
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("xb") as stream:
        stream.write(MAGIC)
        stream.write(struct.pack("<II", ENDIAN_MARKER, version))
        A, Z, nrefmax, J2, parity, emax, e2max, hw = header
        stream.write(struct.pack("<7idQQ", A, Z, nrefmax, J2, parity, emax, e2max, float(hw),
                                 len(orbits), len(lambda_blocks)))
        for value in digests:
            stream.write(value.encode("ascii"))
        for index, n, l, j2, tz2 in orbits:
            stream.write(struct.pack("<5id", index, n, l, j2, tz2, float(occupations[index])))
        stream.write(np.asarray(transformation, dtype="<f8").tobytes(order="C"))
        for J, block_parity, Tz, pairs, matrix in lambda_blocks:
            stream.write(struct.pack("<3iIQ", J, block_parity, Tz, 0, len(pairs)))
            for a, b in pairs:
                stream.write(struct.pack("<II", a, b))
            stream.write(np.asarray(matrix, dtype="<f8").tobytes(order="C"))
        if version == SCHEMA_VERSION:
            return
        energy, label = reference_energy if reference_energy is not None else (0.0, "")
        encoded = label.encode("ascii")
        if len(encoded) > REFERENCE_ENERGY_LABEL_BYTES:
            raise ValueError("reference energy label is longer than 64 bytes")
        stream.write(struct.pack("<Id", int(reference_energy is not None), float(energy)))
        stream.write(encoded.ljust(REFERENCE_ENERGY_LABEL_BYTES, b"\0"))
        blocks = lambda3_blocks if lambda3_blocks is not None else []
        stream.write(struct.pack("<IQ", int(lambda3_blocks is not None), len(blocks)))
        for twoJ, block_parity, twoTz, states, matrix in blocks:
            stream.write(struct.pack("<3iIQ", twoJ, block_parity, twoTz, 0, len(states)))
            for a, b, c, J12 in states:
                if a > b:
                    raise ValueError("lambda3 states must have a <= b")
                stream.write(struct.pack("<3Ii", a, b, c, J12))
            stream.write(np.asarray(matrix, dtype="<f8").tobytes(order="C"))


def export_reference(
    reference_path: str | Path,
    output_path: str | Path,
    *,
    tolerance: float = 1e-10,
    interaction_path: str | Path | None = None,
    lambda3: bool = False,
) -> dict[str, object]:
    """Build densities from the saved wavefunction and write mrimsrg_jref.

    With lambda3=True the three-body cumulant and the reference energy are
    added (version 2); otherwise the version-1 file is unchanged.
    """

    root = Path(reference_path)
    output = Path(output_path)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing MR reference: {output}")
    reference = load_reference(root, interaction_path=interaction_path)
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
    metadata = reference.metadata
    lambda3_blocks = None
    reference_energy = None
    lambda3_report = {}
    if lambda3:
        mstates = [(int(row[0]), int(row[4])) for row in reference.orbits]
        lambda3_blocks, lambda3_report = coupled_lambda3(
            [(orbit.n, orbit.l, orbit.j2, orbit.tz2) for orbit in j_orbits],
            mstates,
            np.asarray(reference.determinants, dtype=np.uint8),
            np.asarray(reference.coefficients, dtype=np.float64),
            transformation,
            occupations,
        )
        reference_energy = (float(metadata["reference_energy"]), "simpleFCI NCSM eigenvalue")
    write_jref(
        output,
        header=(metadata["A"], metadata["Z"], metadata["Nrefmax"], metadata["J2"],
                metadata["parity"], metadata["emax"], metadata["e2max"], float(metadata["hw"])),
        orbits=[(orbit.index, orbit.n, orbit.l, orbit.j2, orbit.tz2) for orbit in j_orbits],
        occupations=occupations,
        transformation=transformation,
        lambda_blocks=[(block.J, block.parity, block.Tz, block.pairs, block.matrix)
                       for block in lambda_blocks],
        digests=(interaction_sha256, rdm_sha256, wavefunction_sha256),
        lambda3_blocks=lambda3_blocks,
        reference_energy=reference_energy,
    )

    return {
        **lambda3_report,
        "schema": "mrimsrg_jref_v2" if lambda3 else "mrimsrg_jref_v1",
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
    parser.add_argument(
        "--interaction",
        type=Path,
        help="relocated byte-identical interaction to verify instead of metadata path",
    )
    parser.add_argument("--tolerance", type=float, default=1e-10)
    parser.add_argument("--lambda3", action="store_true",
                        help="also write the three-body cumulant and reference energy (jref v2)")
    args = parser.parse_args()
    summary = export_reference(
        args.reference,
        args.output,
        tolerance=args.tolerance,
        interaction_path=args.interaction,
        lambda3=args.lambda3,
    )
    for key, value in summary.items():
        print(f"{key}={value}")


if __name__ == "__main__":
    main()
