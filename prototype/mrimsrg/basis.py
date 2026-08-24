"""Small basis-covariance helpers for non-diagonal one-body densities.

The published explicit MR-IMSRG(2) equations used by this prototype are
written in a natural-orbital basis.  A correlated Nrefmax=2 NCSM reference
need not have diagonal ``gamma1`` in the original HO basis.  We therefore
diagonalize only connected off-diagonal blocks of ``gamma1``, evaluate the
natural-basis equations, and transform operators back before applying the HO
``Delta e`` mask.  The externally visible Hamiltonian always remains in the
original HO orbit ordering.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

try:
    from .densities import Densities
    from .normal_order import MRHamiltonian
except ImportError:
    from densities import Densities
    from normal_order import MRHamiltonian


@dataclass(frozen=True)
class NaturalBasis:
    """Transformation from the original basis to a temporary natural basis.

    Columns of ``vectors`` are natural orbitals expressed in the original
    basis, so a one-body matrix transforms as ``U.T @ O @ U``.
    """

    vectors: np.ndarray
    densities: Densities
    is_identity: bool


def _connected_components(adjacency: np.ndarray) -> list[list[int]]:
    remaining = set(range(adjacency.shape[0]))
    components: list[list[int]] = []
    while remaining:
        start = min(remaining)
        stack = [start]
        remaining.remove(start)
        component: list[int] = []
        while stack:
            node = stack.pop()
            component.append(node)
            neighbors = [
                int(index)
                for index in np.flatnonzero(adjacency[node])
                if int(index) in remaining
            ]
            for neighbor in neighbors:
                remaining.remove(neighbor)
                stack.append(neighbor)
        components.append(sorted(component))
    return components


def transform_array(
    values: np.ndarray,
    vectors: np.ndarray,
    *,
    to_natural: bool,
    coefficient_tolerance: float = 1e-14,
) -> np.ndarray:
    """Transform every tensor index with a sparse orthogonal matrix.

    The Nrefmax=2 references mix only small radial blocks.  Applying the
    nonzero matrix entries to one axis at a time keeps this transparent and
    avoids treating the mostly identity transformation as a dense n^5
    contraction.
    """
    if values.ndim not in (2, 4):
        raise ValueError("only one- and two-body coefficient arrays are supported")
    norb = vectors.shape[0]
    if vectors.shape != (norb, norb) or values.shape != (norb,) * values.ndim:
        raise ValueError("basis transformation and tensor shapes are incompatible")
    matrix = vectors.T if to_natural else vectors
    nonzero = np.argwhere(np.abs(matrix) > coefficient_tolerance)
    result = np.asarray(values)
    for axis in range(values.ndim):
        transformed = np.zeros_like(result)
        for output_index, input_index in nonzero:
            output_slice = [slice(None)] * values.ndim
            input_slice = [slice(None)] * values.ndim
            output_slice[axis] = int(output_index)
            input_slice[axis] = int(input_index)
            transformed[tuple(output_slice)] += (
                matrix[output_index, input_index] * result[tuple(input_slice)]
            )
        result = transformed
    return result


def transform_hamiltonian(
    hamiltonian: MRHamiltonian, vectors: np.ndarray, *, to_natural: bool
) -> MRHamiltonian:
    """Transform an MR-normal-ordered operator without changing its scalar."""
    return MRHamiltonian(
        hamiltonian.zero_body,
        transform_array(hamiltonian.one_body, vectors, to_natural=to_natural),
        transform_array(hamiltonian.two_body, vectors, to_natural=to_natural),
    )


def prepare_natural_basis(
    densities: Densities, tolerance: float = 1e-10
) -> NaturalBasis:
    """Diagonalize only connected non-diagonal blocks of ``gamma1``."""
    gamma1 = np.asarray(densities.gamma1)
    if gamma1.ndim != 2 or gamma1.shape[0] != gamma1.shape[1]:
        raise ValueError("gamma1 must be a square matrix")
    if np.max(np.abs(gamma1 - gamma1.T)) > tolerance:
        raise ValueError("gamma1 must be Hermitian")
    norb = gamma1.shape[0]
    off_diagonal = gamma1 - np.diag(np.diag(gamma1))
    adjacency = np.abs(off_diagonal) > tolerance
    if not np.any(adjacency):
        return NaturalBasis(np.eye(norb), densities, True)

    vectors = np.eye(norb)
    for component in _connected_components(adjacency):
        if len(component) == 1:
            continue
        block = gamma1[np.ix_(component, component)]
        occupations, block_vectors = np.linalg.eigh(block)
        order = np.argsort(occupations)[::-1]
        block_vectors = block_vectors[:, order]
        # Remove the arbitrary column sign so output is deterministic.
        for column in range(block_vectors.shape[1]):
            pivot = int(np.argmax(np.abs(block_vectors[:, column])))
            if block_vectors[pivot, column] < 0.0:
                block_vectors[:, column] *= -1.0
        vectors[np.ix_(component, component)] = block_vectors

    gamma1_natural = vectors.T @ gamma1 @ vectors
    gamma2_natural = transform_array(
        densities.gamma2, vectors, to_natural=True
    )
    lambda2_natural = transform_array(
        densities.lambda2, vectors, to_natural=True
    )
    natural_densities = Densities(
        gamma1=gamma1_natural,
        gamma2=gamma2_natural,
        lambda2=lambda2_natural,
    )
    off_diagonal_error = np.max(
        np.abs(gamma1_natural - np.diag(np.diag(gamma1_natural)))
    )
    if off_diagonal_error > 5.0 * tolerance:
        raise ValueError("failed to diagonalize gamma1 in connected blocks")
    return NaturalBasis(vectors, natural_densities, False)
