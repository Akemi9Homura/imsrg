"""Explicit natural-orbital MR-IMSRG(2) commutator.

This implements Eqs. (49)--(51) of Hergert, Phys. Scripta 92, 023002
(2017), arXiv:1607.06882, with the three-body irreducible density set to
zero.  The index order is the same as :mod:`densities` and
:mod:`normal_order`::

    one_body[p,q]       <-> a^dagger_p a_q
    two_body[p,q,r,s]   <-> a^dagger_p a^dagger_q a_s a_r / 4

The equations used here assume a natural-orbital basis, so ``gamma1`` must be
diagonal.  The fixed Nrefmax=0 benchmarks satisfy this condition directly in
the HO m-scheme basis.  Refusing a non-diagonal density is intentional: using
its diagonal alone would silently evaluate the wrong generalized-Wick
contractions.

The individual rank contributions were also regenerated with
``refs/qcombo/examples/MR_IMSRG2.ipynb`` (QCombo 0.2.0).  Runtime code does not
depend on the ignored local reference tree.
"""

from __future__ import annotations

import numpy as np

try:
    from .densities import Densities
    from .normal_order import MRHamiltonian
except ImportError:  # Direct use from the prototype directory.
    from densities import Densities
    from normal_order import MRHamiltonian


def _validate_shapes(
    left: MRHamiltonian, right: MRHamiltonian, densities: Densities
) -> int:
    norb = int(left.one_body.shape[0])
    matrix_shape = (norb, norb)
    tensor_shape = (norb, norb, norb, norb)
    if left.one_body.shape != matrix_shape or right.one_body.shape != matrix_shape:
        raise ValueError("one-body tensors have incompatible shapes")
    if left.two_body.shape != tensor_shape or right.two_body.shape != tensor_shape:
        raise ValueError("two-body tensors have incompatible shapes")
    if densities.gamma1.shape != matrix_shape:
        raise ValueError("gamma1 has an incompatible shape")
    if densities.lambda2.shape != tensor_shape:
        raise ValueError("lambda2 has an incompatible shape")
    return norb


def _natural_occupations(densities: Densities, tolerance: float) -> np.ndarray:
    gamma1 = densities.gamma1
    diagonal = np.diag(np.diag(gamma1))
    off_diagonal_error = float(np.max(np.abs(gamma1 - diagonal)))
    if off_diagonal_error > tolerance:
        raise ValueError(
            "MR-IMSRG(2) equations require a natural-orbital basis: "
            f"max off-diagonal gamma1 = {off_diagonal_error:.3e}"
        )
    occupations = np.diag(gamma1).copy()
    if np.min(occupations) < -tolerance or np.max(occupations) > 1.0 + tolerance:
        raise ValueError("natural occupations lie outside [0,1]")
    return occupations


def commutator(
    left: MRHamiltonian,
    right: MRHamiltonian,
    densities: Densities,
    *,
    natural_tolerance: float = 1e-10,
) -> MRHamiltonian:
    """Return the lambda3=0 MR-IMSRG(2) truncation of ``[left, right]``.

    Zero-body inputs commute with every operator and therefore do not enter.
    No symmetry projection is applied to the result; the tests intentionally
    check that the equations themselves preserve the expected symmetries.
    """
    contributions = commutator_contributions(
        left, right, densities, natural_tolerance=natural_tolerance
    )
    return MRHamiltonian(
        sum(term.zero_body for term in contributions.values()),
        sum(term.one_body for term in contributions.values()),
        sum(term.two_body for term in contributions.values()),
    )


def commutator_contributions(
    left: MRHamiltonian,
    right: MRHamiltonian,
    densities: Densities,
    *,
    natural_tolerance: float = 1e-10,
) -> dict[str, MRHamiltonian]:
    """Return the named contractions used by :func:`commutator`.

    The eight ``comm...ss`` keys are the scalar IMSRG(2) decomposition used
    by the current ``src/Commutator.cc``.  The two ``mr_lambda2`` keys are the
    additional contractions in Hergert Eqs. (49)--(50); they vanish exactly
    in the single-reference limit.  This is a diagnostic view of the actual
    production calculation: :func:`commutator` sums this mapping rather than
    evaluating a second implementation.
    """
    norb = _validate_shapes(left, right, densities)
    n = _natural_occupations(densities, natural_tolerance)
    nbar = 1.0 - n
    lambda2 = densities.lambda2
    x1, x2 = left.one_body, left.two_body
    y1, y2 = right.one_body, right.two_body

    # Eq. (51): two-body part.  This is evaluated first because Eq. (49)
    # contains 1/4 C2^{ab}_{cd} lambda2^{ab}_{cd}.
    c2_122 = (
        np.einsum("ia,ajkl->ijkl", x1, y2, optimize=True)
        + np.einsum("ja,iakl->ijkl", x1, y2, optimize=True)
        - np.einsum("ak,ijal->ijkl", x1, y2, optimize=True)
        - np.einsum("al,ijka->ijkl", x1, y2, optimize=True)
        - np.einsum("ia,ajkl->ijkl", y1, x2, optimize=True)
        - np.einsum("ja,iakl->ijkl", y1, x2, optimize=True)
        + np.einsum("ak,ijal->ijkl", y1, x2, optimize=True)
        + np.einsum("al,ijka->ijkl", y1, x2, optimize=True)
    )
    pair_weight = 1.0 - n[:, None] - n[None, :]
    c2_222_pp_hh = 0.5 * (
        np.einsum("ijab,abkl,ab->ijkl", x2, y2, pair_weight, optimize=True)
        - np.einsum("ijab,abkl,ab->ijkl", y2, x2, pair_weight, optimize=True)
    )
    occupation_difference = n[:, None] - n[None, :]
    c2_222_ph = (
        np.einsum(
            "iakb,jbla,ab->ijkl", x2, y2, occupation_difference, optimize=True
        )
        - np.einsum(
            "iakb,jbla,ab->ijkl", y2, x2, occupation_difference, optimize=True
        )
        - np.einsum(
            "jakb,ibla,ab->ijkl", x2, y2, occupation_difference, optimize=True
        )
        + np.einsum(
            "jakb,ibla,ab->ijkl", y2, x2, occupation_difference, optimize=True
        )
    )

    # Eq. (50): one-body part.
    c1_111 = np.einsum("ia,aj->ij", x1, y1, optimize=True) - np.einsum(
        "ia,aj->ij", y1, x1, optimize=True
    )
    c1_121 = np.einsum(
        "ab,biaj,ab->ij", x1, y2, occupation_difference, optimize=True
    ) - np.einsum(
        "ab,biaj,ab->ij", y1, x2, occupation_difference, optimize=True
    )
    three_index_weight = (
        n[:, None, None] * nbar[None, :, None] * nbar[None, None, :]
        + nbar[:, None, None] * n[None, :, None] * n[None, None, :]
    )
    c1_221 = 0.5 * (
        np.einsum("iabc,bcja,abc->ij", x2, y2, three_index_weight, optimize=True)
        - np.einsum("iabc,bcja,abc->ij", y2, x2, three_index_weight, optimize=True)
    )

    c1_mr_lambda2 = np.zeros((norb, norb), dtype=np.result_type(x1, x2, y1, y2))
    if np.any(lambda2):
        c1_mr_lambda2 += 0.25 * (
            np.einsum("iabc,deja,debc->ij", x2, y2, lambda2, optimize=True)
            - np.einsum("iabc,deja,debc->ij", y2, x2, lambda2, optimize=True)
        )
        c1_mr_lambda2 += np.einsum(
            "iabc,bejd,aecd->ij", x2, y2, lambda2, optimize=True
        ) - np.einsum(
            "iabc,bejd,aecd->ij", y2, x2, lambda2, optimize=True
        )
        c1_mr_lambda2 -= 0.5 * (
            np.einsum("iajb,cdae,cdbe->ij", x2, y2, lambda2, optimize=True)
            - np.einsum("iajb,cdae,cdbe->ij", y2, x2, lambda2, optimize=True)
        )
        c1_mr_lambda2 += 0.5 * (
            np.einsum("iajb,bcde,acde->ij", x2, y2, lambda2, optimize=True)
            - np.einsum("iajb,bcde,acde->ij", y2, x2, lambda2, optimize=True)
        )

    # Eq. (49), with lambda3=0.
    c0_110 = np.einsum(
        "ab,ba,ab->", x1, y1, occupation_difference, optimize=True
    )
    c0_220 = 0.25 * (
        np.einsum("abcd,cdab,a,b,c,d->", x2, y2, n, n, nbar, nbar, optimize=True)
        - np.einsum(
            "abcd,cdab,a,b,c,d->", y2, x2, n, n, nbar, nbar, optimize=True
        )
    )
    c2 = c2_122 + c2_222_pp_hh + c2_222_ph
    c0_mr_lambda2 = 0.25 * np.einsum(
        "abcd,abcd->", c2, lambda2, optimize=True
    )

    zero1 = np.zeros((norb, norb), dtype=np.result_type(x1, x2, y1, y2))
    zero2 = np.zeros((norb, norb, norb, norb), dtype=zero1.dtype)

    def term(
        *,
        zero_body: complex | float = 0.0,
        one_body: np.ndarray = zero1,
        two_body: np.ndarray = zero2,
    ) -> MRHamiltonian:
        return MRHamiltonian(np.real_if_close(zero_body).item(), one_body, two_body)

    return {
        "comm110ss": term(zero_body=c0_110),
        "comm220ss": term(zero_body=c0_220),
        "comm111ss": term(one_body=c1_111),
        "comm121ss": term(one_body=c1_121),
        "comm221ss": term(one_body=c1_221),
        "comm122ss": term(two_body=c2_122),
        "comm222_pp_hhss": term(two_body=c2_222_pp_hh),
        "comm222_phss": term(two_body=c2_222_ph),
        "mr_lambda2_one_body": term(one_body=c1_mr_lambda2),
        "mr_lambda2_zero_body": term(zero_body=c0_mr_lambda2),
    }
