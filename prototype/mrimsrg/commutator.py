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
    _validate_shapes(left, right, densities)
    n = _natural_occupations(densities, natural_tolerance)
    nbar = 1.0 - n
    lambda2 = densities.lambda2
    x1, x2 = left.one_body, left.two_body
    y1, y2 = right.one_body, right.two_body

    # Eq. (51): two-body part.  This is evaluated first because Eq. (49)
    # contains 1/4 C2^{ab}_{cd} lambda2^{ab}_{cd}.
    c2 = (
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
    c2 += 0.5 * (
        np.einsum("ijab,abkl,ab->ijkl", x2, y2, pair_weight, optimize=True)
        - np.einsum("ijab,abkl,ab->ijkl", y2, x2, pair_weight, optimize=True)
    )
    occupation_difference = n[:, None] - n[None, :]
    c2 += (
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
    c1 = np.einsum("ia,aj->ij", x1, y1, optimize=True) - np.einsum(
        "ia,aj->ij", y1, x1, optimize=True
    )
    c1 += np.einsum(
        "ab,biaj,ab->ij", x1, y2, occupation_difference, optimize=True
    ) - np.einsum(
        "ab,biaj,ab->ij", y1, x2, occupation_difference, optimize=True
    )
    three_index_weight = (
        n[:, None, None] * nbar[None, :, None] * nbar[None, None, :]
        + nbar[:, None, None] * n[None, :, None] * n[None, None, :]
    )
    c1 += 0.5 * (
        np.einsum("iabc,bcja,abc->ij", x2, y2, three_index_weight, optimize=True)
        - np.einsum("iabc,bcja,abc->ij", y2, x2, three_index_weight, optimize=True)
    )

    if np.any(lambda2):
        c1 += 0.25 * (
            np.einsum("iabc,deja,debc->ij", x2, y2, lambda2, optimize=True)
            - np.einsum("iabc,deja,debc->ij", y2, x2, lambda2, optimize=True)
        )
        c1 += np.einsum(
            "iabc,bejd,aecd->ij", x2, y2, lambda2, optimize=True
        ) - np.einsum(
            "iabc,bejd,aecd->ij", y2, x2, lambda2, optimize=True
        )
        c1 -= 0.5 * (
            np.einsum("iajb,cdae,cdbe->ij", x2, y2, lambda2, optimize=True)
            - np.einsum("iajb,cdae,cdbe->ij", y2, x2, lambda2, optimize=True)
        )
        c1 += 0.5 * (
            np.einsum("iajb,bcde,acde->ij", x2, y2, lambda2, optimize=True)
            - np.einsum("iajb,bcde,acde->ij", y2, x2, lambda2, optimize=True)
        )

    # Eq. (49), with lambda3=0.
    c0 = np.einsum(
        "ab,ba,ab->", x1, y1, occupation_difference, optimize=True
    )
    c0 += 0.25 * (
        np.einsum("abcd,cdab,a,b,c,d->", x2, y2, n, n, nbar, nbar, optimize=True)
        - np.einsum(
            "abcd,cdab,a,b,c,d->", y2, x2, n, n, nbar, nbar, optimize=True
        )
    )
    c0 += 0.25 * np.einsum("abcd,abcd->", c2, lambda2, optimize=True)

    scalar = np.real_if_close(c0).item()
    return MRHamiltonian(scalar, c1, c2)

