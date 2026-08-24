"""Single Brillouin generator with the relaxed IM-NCSM decoupling mask.

The matrix elements are Eqs. (58)--(59) of Hergert, Phys. Scripta 92,
023002 (2017), with ``lambda3=0``.  They are the irreducible Brillouin
residuals ``<Psi|[H,{A}]|Psi>`` and therefore form an anti-Hermitian
generator without energy denominators.  The same expressions are generated
by ``refs/qcombo/examples/Brillouin.ipynb``.

Only matrix elements that change the sum of HO single-particle quanta are
retained.  This is the relaxed IM-NCSM pattern of Vobig Eqs. (6.5.13)--
(6.5.16) and Mongelli Eqs. (5.218)--(5.221).
"""

from __future__ import annotations

import numpy as np

try:
    from .commutator import _natural_occupations
    from .densities import Densities
    from .normal_order import MRHamiltonian
except ImportError:
    from commutator import _natural_occupations
    from densities import Densities
    from normal_order import MRHamiltonian


def oscillator_quanta_from_orbits(orbits: np.ndarray) -> np.ndarray:
    """Return ``e=2n+l`` from bridge columns ``(jindex,n,l,2j,2m,2tz)``."""
    if orbits.ndim != 2 or orbits.shape[1] < 3:
        raise ValueError("orbit table must contain n and l columns")
    return 2 * orbits[:, 1].astype(np.int64) + orbits[:, 2].astype(np.int64)


def decoupling_masks(oscillator_quanta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Construct the relaxed one- and two-body IM-NCSM masks."""
    quanta = np.asarray(oscillator_quanta)
    if quanta.ndim != 1:
        raise ValueError("oscillator quanta must be a one-dimensional array")
    one_body = quanta[:, None] != quanta[None, :]
    pair_quanta = quanta[:, None] + quanta[None, :]
    two_body = pair_quanta[:, :, None, None] != pair_quanta[None, None, :, :]
    return one_body, two_body


def brillouin_generator(
    hamiltonian: MRHamiltonian,
    densities: Densities,
    oscillator_quanta: np.ndarray,
    *,
    natural_tolerance: float = 1e-10,
) -> MRHamiltonian:
    """Build the masked, lambda3=0 Brillouin generator."""
    n = _natural_occupations(densities, natural_tolerance)
    nbar = 1.0 - n
    f = hamiltonian.one_body
    gamma = hamiltonian.two_body
    lambda2 = densities.lambda2
    norb = len(n)
    if f.shape != (norb, norb) or gamma.shape != (norb,) * 4:
        raise ValueError("Hamiltonian and density shapes are incompatible")
    mask1, mask2 = decoupling_masks(oscillator_quanta)
    if mask1.shape != f.shape:
        raise ValueError("oscillator-quanta array has an incompatible length")

    # QCombo form of Eq. (58), with output indices (i,j).
    eta1 = (n[None, :] - n[:, None]) * f.T
    if np.any(lambda2):
        eta1 += 0.5 * (
            np.einsum("abci,abcj->ij", gamma, lambda2, optimize=True)
            - np.einsum("ajcd,aicd->ij", gamma, lambda2, optimize=True)
        )

    # Eq. (59), first line and the f-lambda2 terms.
    forward_weight = (
        nbar[:, None, None, None]
        * nbar[None, :, None, None]
        * n[None, None, :, None]
        * n[None, None, None, :]
    )
    reverse_weight = (
        n[:, None, None, None]
        * n[None, :, None, None]
        * nbar[None, None, :, None]
        * nbar[None, None, None, :]
    )
    eta2 = gamma.transpose(2, 3, 0, 1) * (forward_weight - reverse_weight)

    if np.any(lambda2):
        eta2 += (
            np.einsum("ai,ajkl->ijkl", f, lambda2, optimize=True)
            - np.einsum("aj,aikl->ijkl", f, lambda2, optimize=True)
            - np.einsum("ka,ijal->ijkl", f, lambda2, optimize=True)
            + np.einsum("la,ijak->ijkl", f, lambda2, optimize=True)
        )
        eta2 += 0.5 * (
            np.einsum("abij,abkl,ij->ijkl", gamma, lambda2, 1.0 - n[:, None] - n[None, :], optimize=True)
            - np.einsum("klab,ijab,kl->ijkl", gamma, lambda2, 1.0 - n[:, None] - n[None, :], optimize=True)
        )

        # (1-Pij)(1-Pkl) sum_ac (n_j-n_k) Gamma^{ak}_{cj} lambda^{ai}_{cl}.
        # Write all four permutations explicitly so every exchanged index is
        # visible and independently testable.
        eta2 += np.einsum(
            "akcj,aicl,jk->ijkl", gamma, lambda2, n[:, None] - n[None, :], optimize=True
        )
        eta2 -= np.einsum(
            "akci,ajcl,ik->ijkl", gamma, lambda2, n[:, None] - n[None, :], optimize=True
        )
        eta2 -= np.einsum(
            "alcj,aick,jl->ijkl", gamma, lambda2, n[:, None] - n[None, :], optimize=True
        )
        eta2 += np.einsum(
            "alci,ajck,il->ijkl", gamma, lambda2, n[:, None] - n[None, :], optimize=True
        )

    eta1 *= mask1
    eta2 *= mask2
    return MRHamiltonian(0.0, eta1, eta2)


def masked_residual_norm(generator: MRHamiltonian) -> float:
    """Frobenius norm used to monitor the selected Brillouin residuals."""
    return float(
        np.sqrt(
            np.vdot(generator.one_body, generator.one_body).real
            + 0.25 * np.vdot(generator.two_body, generator.two_body).real
        )
    )
