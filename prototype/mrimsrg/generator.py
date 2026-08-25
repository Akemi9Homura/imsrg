"""Single modified White generator for the rapid IM-NCSM flow.

The directional decoupling matrix elements ``D1`` and ``D2`` are Mongelli
Eqs. (5.216)--(5.217), with ``lambda3=0`` and terms quadratic in ``lambda2``
omitted. Their anti-Hermitian combination reproduces the Brillouin residual
of Hergert Eqs. (58)--(59); the positive products were also regenerated with
QCombo's generalized Wick expansion.  They are retained as an independent,
strict decoupling diagnostic.

The one generator used by the flow is the practical White-NCSM generator of
Vobig Sec. 6.5.4.  In addition to omitting the explicitly indicated
``O(lambda2)`` denominator corrections, that implementation defines
"White-NCSM" by neglecting *all* irreducible-density terms in the generator
matrix elements.  The MR-IMSRG(2) commutator itself still retains
``lambda2``.  A sign-preserving cutoff follows the established ``imsrg++``
handling of small denominators.

Only natural-orbital matrix elements whose inherited labels change the sum of
HO single-particle quanta are retained, i.e. the relaxed IM-NCSM pattern of
Vobig Eqs. (6.5.13)--(6.5.27) and Mongelli Eqs. (5.218)--(5.221).
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


WHITE_DENOMINATOR_CUTOFF = 1e-6
GENERATOR_IMPLEMENTATION = "white_ncsm_spherical_monopole_v1"


def oscillator_quanta_from_orbits(orbits: np.ndarray) -> np.ndarray:
    """Return ``e=2n+l`` from bridge columns ``(jindex,n,l,2j,2m,2tz)``."""
    if orbits.ndim != 2 or orbits.shape[1] < 3:
        raise ValueError("orbit table must contain n and l columns")
    return 2 * orbits[:, 1].astype(np.int64) + orbits[:, 2].astype(np.int64)


def spherical_orbit_groups_from_orbits(orbits: np.ndarray) -> np.ndarray:
    """Return and validate spherical-orbit labels for an m-scheme orbit table.

    The bridge columns are ``(jindex,n,l,2j,2m,2tz)``.  Every ``jindex``
    must contain one complete magnetic multiplet with common ``n,l,j,tz``.
    Natural orbitals inherit these output-slot labels after radial mixing.
    """
    values = np.asarray(orbits)
    if values.ndim != 2 or values.shape[1] < 6:
        raise ValueError("orbit table must contain (jindex,n,l,2j,2m,2tz)")
    groups = values[:, 0].astype(np.int64)
    if not np.array_equal(groups, values[:, 0]):
        raise ValueError("spherical j-orbit indices must be integers")
    for group in np.unique(groups):
        members = values[groups == group]
        if not np.all(members[:, [1, 2, 3, 5]] == members[0, [1, 2, 3, 5]]):
            raise ValueError(f"j-orbit group {group} mixes n,l,j, or tz")
        j2 = int(members[0, 3])
        expected_m2 = np.arange(-j2, j2 + 1, 2, dtype=np.int64)
        if not np.array_equal(np.sort(members[:, 4].astype(np.int64)), expected_m2):
            raise ValueError(
                f"j-orbit group {group} is not a complete magnetic multiplet"
            )
    return groups


def _validate_spherical_groups(
    spherical_orbit_groups: np.ndarray | None, norb: int
) -> np.ndarray:
    if spherical_orbit_groups is None:
        return np.arange(norb, dtype=np.int64)
    groups = np.asarray(spherical_orbit_groups)
    if groups.shape != (norb,):
        raise ValueError("spherical-orbit group array has an incompatible length")
    _, inverse = np.unique(groups, return_inverse=True)
    return inverse.astype(np.int64)


def _spherical_scalar_diagonals(
    hamiltonian: MRHamiltonian,
    occupations: np.ndarray,
    spherical_orbit_groups: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return m-independent occupations, 1B diagonals, and 2B monopoles.

    For spherical orbits ``a,b``, averaging ``Gamma[p,q,p,q]`` over every
    ``m_p in a, m_q in b`` is exactly the unnormalized monopole convention
    used by ``TwoBodyME::GetTBMEmonopole``:

    ``sqrt((1+d_ab)^2) sum_J (2J+1) Gamma^J_abab / (d_a d_b)``.

    The ordered-m average automatically supplies the factor two for
    identical spherical orbits because antisymmetric ``p=q`` elements vanish.
    """
    norb = len(occupations)
    groups = _validate_spherical_groups(spherical_orbit_groups, norb)
    ngroup = int(groups.max()) + 1 if norb else 0
    counts = np.bincount(groups, minlength=ngroup).astype(np.float64)

    occupation_mean = np.bincount(
        groups, weights=np.asarray(occupations, dtype=np.float64), minlength=ngroup
    ) / counts
    one_diagonal = np.diag(hamiltonian.one_body)
    one_mean = np.bincount(
        groups, weights=np.asarray(one_diagonal, dtype=np.float64), minlength=ngroup
    ) / counts

    m_diagonal = np.einsum("ijij->ij", hamiltonian.two_body)
    pair_groups = groups[:, None] * ngroup + groups[None, :]
    pair_sum = np.bincount(
        pair_groups.ravel(),
        weights=np.asarray(m_diagonal, dtype=np.float64).ravel(),
        minlength=ngroup * ngroup,
    ).reshape(ngroup, ngroup)
    pair_count = counts[:, None] * counts[None, :]
    monopole = pair_sum / pair_count
    return (
        occupation_mean[groups],
        one_mean[groups],
        monopole[groups[:, None], groups[None, :]],
    )


def decoupling_masks(oscillator_quanta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Construct the relaxed one- and two-body IM-NCSM masks."""
    quanta = np.asarray(oscillator_quanta)
    if quanta.ndim != 1:
        raise ValueError("oscillator quanta must be a one-dimensional array")
    one_body = quanta[:, None] != quanta[None, :]
    pair_quanta = quanta[:, None] + quanta[None, :]
    two_body = pair_quanta[:, :, None, None] != pair_quanta[None, None, :, :]
    return one_body, two_body


def decoupling_matrix_elements(
    hamiltonian: MRHamiltonian,
    densities: Densities,
    *,
    natural_tolerance: float = 1e-10,
) -> tuple[np.ndarray, np.ndarray]:
    """Return directional ``D1=<H A1>`` and ``D2=<H A2>`` elements.

    Array indices follow this prototype's coefficient convention. In
    particular, the leading one-body term is
    ``D1[i,j] = n[j] (1-n[i]) f[j,i]``. This apparent transpose relative to
    some printed formulas is fixed by the convention that ``f[p,q]``
    multiplies ``a^dagger_p a_q``.
    """
    n = _natural_occupations(densities, natural_tolerance)
    nbar = 1.0 - n
    f = hamiltonian.one_body
    gamma = hamiltonian.two_body
    lambda2 = densities.lambda2
    norb = len(n)
    if f.shape != (norb, norb) or gamma.shape != (norb,) * 4:
        raise ValueError("Hamiltonian and density shapes are incompatible")

    # QCombo positive products H(1B) A(1B) and H(2B) A(1B), after setting
    # lambda3=0. The f-lambda2 term is Hermitian for Hermitian inputs, but it
    # must be retained before applying the directional White denominators.
    d1 = n[None, :] * nbar[:, None] * f.T
    if np.any(lambda2):
        d1 += np.einsum("pq,piqj->ij", f, lambda2, optimize=True)
        d1 += 0.5 * (
            np.einsum("i,abci,abcj->ij", nbar, gamma, lambda2, optimize=True)
            - np.einsum("j,ajcd,aicd->ij", n, gamma, lambda2, optimize=True)
        )

    # Mongelli Eq. (5.217). Each permutation is explicit so the two
    # antisymmetric index pairs can be checked independently.
    d2 = (
        nbar[:, None, None, None]
        * nbar[None, :, None, None]
        * n[None, None, :, None]
        * n[None, None, None, :]
        * gamma.transpose(2, 3, 0, 1)
    )
    if np.any(lambda2):
        d2 += np.einsum("l,lp,ijpk->ijkl", n, f, lambda2, optimize=True)
        d2 -= np.einsum("k,kp,ijpl->ijkl", n, f, lambda2, optimize=True)
        d2 -= np.einsum("j,pj,pikl->ijkl", nbar, f, lambda2, optimize=True)
        d2 += np.einsum("i,pi,pjkl->ijkl", nbar, f, lambda2, optimize=True)
        d2 += 0.5 * np.einsum(
            "i,j,pqij,pqkl->ijkl", nbar, nbar, gamma, lambda2, optimize=True
        )
        d2 += 0.5 * np.einsum(
            "k,l,klpq,ijpq->ijkl", n, n, gamma, lambda2, optimize=True
        )
        mixed = np.einsum(
            "k,i,pkiq,pjql->ijkl", n, nbar, gamma, lambda2, optimize=True
        )
        d2 += (
            -mixed
            + mixed.swapaxes(0, 1)
            + mixed.swapaxes(2, 3)
            - mixed.swapaxes(0, 1).swapaxes(2, 3)
        )
    return d1, d2


def white_ncsm_matrix_elements(
    hamiltonian: MRHamiltonian,
    densities: Densities,
    *,
    natural_tolerance: float = 1e-10,
    spherical_orbit_groups: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the lambda-free numerator used by White-NCSM.

    Vobig Sec. 6.5.4 defines the production ``White-NCSM`` choice by
    neglecting all terms involving irreducible densities in the generator.
    These are therefore the leading terms of Eqs. (6.5.9)--(6.5.10), not the
    strict ``D1/D2`` diagnostic returned by :func:`decoupling_matrix_elements`.
    """
    n = _natural_occupations(densities, natural_tolerance)
    n, _, _ = _spherical_scalar_diagonals(
        hamiltonian, n, spherical_orbit_groups
    )
    nbar = 1.0 - n
    norb = len(n)
    if hamiltonian.one_body.shape != (norb, norb):
        raise ValueError("Hamiltonian and density shapes are incompatible")
    if hamiltonian.two_body.shape != (norb,) * 4:
        raise ValueError("Hamiltonian and density shapes are incompatible")

    d1 = n[None, :] * nbar[:, None] * hamiltonian.one_body.T
    d2 = (
        nbar[:, None, None, None]
        * nbar[None, :, None, None]
        * n[None, None, :, None]
        * n[None, None, None, :]
        * hamiltonian.two_body.transpose(2, 3, 0, 1)
    )
    return d1, d2


def epstein_nesbet_denominators(
    hamiltonian: MRHamiltonian,
    densities: Densities,
    *,
    natural_tolerance: float = 1e-10,
    spherical_orbit_groups: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return leading MR Epstein--Nesbet denominators.

    These are Mongelli Eqs. (5.213)--(5.214), without their
    ``O(lambda2)`` terms. In the Slater limit the two-body expression reduces
    exactly to the Epstein--Nesbet denominator in ``src/Generator.cc``.
    """
    n = _natural_occupations(densities, natural_tolerance)
    n, f_diagonal, gamma_diagonal = _spherical_scalar_diagonals(
        hamiltonian, n, spherical_orbit_groups
    )
    nbar = 1.0 - n
    energy = hamiltonian.zero_body

    delta1 = (
        -(nbar[:, None] ** 2)
        * (n[None, :] ** 2)
        * gamma_diagonal
        + (nbar[:, None] ** 2) * n[None, :] * f_diagonal[:, None]
        - nbar[:, None] * (n[None, :] ** 2) * f_diagonal[None, :]
        + energy * (nbar[:, None] * n[None, :] - 1.0)
    )

    weight = (
        nbar[:, None, None, None]
        * nbar[None, :, None, None]
        * n[None, None, :, None]
        * n[None, None, None, :]
    )
    first_ordering = weight * (
        0.5
        * nbar[:, None, None, None]
        * nbar[None, :, None, None]
        * gamma_diagonal[:, :, None, None]
        + 0.5
        * n[None, None, :, None]
        * n[None, None, None, :]
        * gamma_diagonal[None, None, :, :]
        - nbar[:, None, None, None]
        * n[None, None, None, :]
        * gamma_diagonal[:, None, None, :]
        - nbar[:, None, None, None]
        * n[None, None, :, None]
        * gamma_diagonal[:, None, :, None]
        + nbar[:, None, None, None] * f_diagonal[:, None, None, None]
        - n[None, None, :, None] * f_diagonal[None, None, :, None]
    )
    first_ordering += 0.5 * energy * (weight - 1.0)
    delta2 = first_ordering + first_ordering.swapaxes(0, 1).swapaxes(2, 3)
    return delta1, delta2


def _safe_denominator(values: np.ndarray, cutoff: float) -> np.ndarray:
    if cutoff <= 0.0:
        raise ValueError("denominator cutoff must be positive")
    signs = np.where(values < 0.0, -1.0, 1.0)
    return np.where(np.abs(values) < cutoff, signs * cutoff, values)


def masked_decoupling_residual(
    hamiltonian: MRHamiltonian,
    densities: Densities,
    oscillator_quanta: np.ndarray,
    *,
    natural_tolerance: float = 1e-10,
) -> MRHamiltonian:
    """Build the masked anti-Hermitian ``D-D^dagger`` residual."""
    residual = decoupling_residual(
        hamiltonian, densities, natural_tolerance=natural_tolerance
    )
    mask1, mask2 = decoupling_masks(oscillator_quanta)
    if mask1.shape != residual.one_body.shape:
        raise ValueError("oscillator-quanta array has an incompatible length")
    return MRHamiltonian(
        0.0, residual.one_body * mask1, residual.two_body * mask2
    )


def decoupling_residual(
    hamiltonian: MRHamiltonian,
    densities: Densities,
    *,
    natural_tolerance: float = 1e-10,
) -> MRHamiltonian:
    """Build the unmasked anti-Hermitian ``D-D^dagger`` residual."""
    d1, d2 = decoupling_matrix_elements(
        hamiltonian, densities, natural_tolerance=natural_tolerance
    )
    return MRHamiltonian(
        0.0,
        d1 - d1.T,
        d2 - d2.transpose(2, 3, 0, 1),
    )


def white_ncsm_numerator_residual(
    hamiltonian: MRHamiltonian,
    densities: Densities,
    *,
    natural_tolerance: float = 1e-10,
    spherical_orbit_groups: np.ndarray | None = None,
) -> MRHamiltonian:
    """Return the unmasked anti-Hermitian White-NCSM numerator.

    This is recorded separately from the strict, lambda2-dependent
    decoupling residual and from the denominator-weighted anti-Hermitian
    generator.  Its masked norm is a diagnostic, not the formal acceptance
    condition: denominator weighting and anti-Hermitization do not commute
    for fractional occupations.
    """
    d1, d2 = white_ncsm_matrix_elements(
        hamiltonian,
        densities,
        natural_tolerance=natural_tolerance,
        spherical_orbit_groups=spherical_orbit_groups,
    )
    return MRHamiltonian(
        0.0,
        d1 - d1.T,
        d2 - d2.transpose(2, 3, 0, 1),
    )


def white_generator(
    hamiltonian: MRHamiltonian,
    densities: Densities,
    oscillator_quanta: np.ndarray,
    *,
    denominator_cutoff: float = WHITE_DENOMINATOR_CUTOFF,
    natural_tolerance: float = 1e-10,
    spherical_orbit_groups: np.ndarray | None = None,
) -> MRHamiltonian:
    """Build the single masked White generator used by the prototype."""
    eta = white_generator_unmasked(
        hamiltonian,
        densities,
        denominator_cutoff=denominator_cutoff,
        natural_tolerance=natural_tolerance,
        spherical_orbit_groups=spherical_orbit_groups,
    )
    mask1, mask2 = decoupling_masks(oscillator_quanta)
    if mask1.shape != eta.one_body.shape:
        raise ValueError("oscillator-quanta array has an incompatible length")
    return MRHamiltonian(0.0, eta.one_body * mask1, eta.two_body * mask2)


def white_generator_unmasked(
    hamiltonian: MRHamiltonian,
    densities: Densities,
    *,
    denominator_cutoff: float = WHITE_DENOMINATOR_CUTOFF,
    natural_tolerance: float = 1e-10,
    spherical_orbit_groups: np.ndarray | None = None,
) -> MRHamiltonian:
    """Build the unmasked White generator in a natural-orbital basis."""
    d1, d2 = white_ncsm_matrix_elements(
        hamiltonian,
        densities,
        natural_tolerance=natural_tolerance,
        spherical_orbit_groups=spherical_orbit_groups,
    )
    delta1, delta2 = epstein_nesbet_denominators(
        hamiltonian,
        densities,
        natural_tolerance=natural_tolerance,
        spherical_orbit_groups=spherical_orbit_groups,
    )
    weighted1 = d1 / _safe_denominator(delta1, denominator_cutoff)
    weighted2 = d2 / _safe_denominator(delta2, denominator_cutoff)
    eta1 = weighted1 - weighted1.T
    eta2 = weighted2 - weighted2.transpose(2, 3, 0, 1)
    return MRHamiltonian(0.0, eta1, eta2)


def masked_residual_norm(residual: MRHamiltonian) -> float:
    """Frobenius norm used to monitor the selected decoupling residuals."""
    return float(
        np.sqrt(
            np.vdot(residual.one_body, residual.one_body).real
            + 0.25 * np.vdot(residual.two_body, residual.two_body).real
        )
    )
