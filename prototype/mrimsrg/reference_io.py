"""Read and validate the shell-model-obs input-bridge bundle."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np


FIXED_INTERACTION_SHA256 = (
    "76b7243ef53d30955c0293d29da73688dc3839942143ccf147739108bb58ff84"
)


@dataclass(frozen=True)
class ReferenceData:
    metadata: dict[str, Any]
    orbits: np.ndarray
    one_body: np.ndarray
    two_body: np.ndarray
    determinants: np.ndarray
    coefficients: np.ndarray

    @property
    def norb(self) -> int:
        return int(self.orbits.shape[0])


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load_reference(path: str | Path, verify_interaction: bool = True) -> ReferenceData:
    root = Path(path)
    with (root / "metadata.json").open(encoding="utf-8") as stream:
        metadata = json.load(stream)
    data = ReferenceData(
        metadata=metadata,
        orbits=np.load(root / "orbits.npy"),
        one_body=np.load(root / "one_body.npy"),
        two_body=np.load(root / "two_body.npy"),
        determinants=np.load(root / "determinants.npy"),
        coefficients=np.load(root / "coefficients.npy"),
    )
    validate_reference(data, verify_interaction=verify_interaction)
    return data


def validate_reference(data: ReferenceData, verify_interaction: bool = True) -> None:
    metadata = data.metadata
    _require(metadata.get("schema") == "mrimsrg_reference_v1", "unsupported reference schema")
    _require(metadata.get("hw") == 20, "reference must use hw=20 MeV")
    _require(metadata.get("emax") == 2, "reference must use emax=2")
    _require(metadata.get("e2max") == 4, "reference must use e2max=4")
    _require(metadata.get("J2") == 0 and metadata.get("parity") == 1, "reference must be J=0 positive parity")
    _require(metadata.get("interaction_sha256") == FIXED_INTERACTION_SHA256, "metadata has an unexpected interaction SHA-256")
    _require(metadata.get("one_body_convention") == "t[p,q] a^dagger_p a_q", "unexpected one-body convention")
    _require(
        metadata.get("two_body_convention")
        == "(1/4) V[p,q,r,s] a^dagger_p a^dagger_q a_s a_r",
        "unexpected two-body convention",
    )

    norb = data.norb
    nconfig = int(data.determinants.shape[0])
    _require(data.orbits.shape == (norb, 6), "orbit table must have columns (jindex,n,l,2j,2m,2tz)")
    _require(data.one_body.shape == (norb, norb), "invalid one-body tensor shape")
    _require(data.two_body.shape == (norb, norb, norb, norb), "invalid two-body tensor shape")
    _require(data.determinants.shape == (nconfig, norb), "invalid determinant array shape")
    _require(data.coefficients.shape == (nconfig,), "invalid coefficient array shape")
    _require(np.all((data.determinants == 0) | (data.determinants == 1)), "determinants are not binary")
    _require(abs(float(data.coefficients @ data.coefficients) - 1.0) < 1e-10, "reference wavefunction is not normalized")

    particle_counts = data.determinants.sum(axis=1)
    _require(np.all(particle_counts == metadata["A"]), "determinant particle number does not equal A")
    proton_mask = data.orbits[:, 5] == -1
    neutron_mask = data.orbits[:, 5] == 1
    _require(np.all(data.determinants[:, proton_mask].sum(axis=1) == metadata["Z"]), "determinant proton number does not equal Z")
    _require(np.all(data.determinants[:, neutron_mask].sum(axis=1) == metadata["N"]), "determinant neutron number does not equal N")

    tolerance = 2e-12
    _require(np.max(np.abs(data.one_body - data.one_body.T)) < tolerance, "one-body Hamiltonian is not Hermitian")
    _require(np.max(np.abs(data.two_body + data.two_body.swapaxes(0, 1))) < tolerance, "two-body tensor is not antisymmetric in bra indices")
    _require(np.max(np.abs(data.two_body + data.two_body.swapaxes(2, 3))) < tolerance, "two-body tensor is not antisymmetric in ket indices")
    _require(np.max(np.abs(data.two_body - data.two_body.transpose(2, 3, 0, 1))) < tolerance, "two-body Hamiltonian is not Hermitian")

    if verify_interaction:
        interaction = Path(metadata["interaction"])
        _require(interaction.is_file(), f"interaction file is unavailable: {interaction}")
        actual_sha256 = sha256_file(interaction)
        _require(actual_sha256 == FIXED_INTERACTION_SHA256, f"unexpected interaction SHA-256: {actual_sha256}")
