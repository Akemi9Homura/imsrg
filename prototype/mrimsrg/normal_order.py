"""MR normal ordering for a vacuum 0B+1B+2B Hamiltonian.

Array conventions match :mod:`densities` and the C++ input bridge.  In
particular, ``V[p,q,r,s]`` multiplies
``a^dagger_p a^dagger_q a_s a_r / 4``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

try:
    from .densities import Densities
except ImportError:
    from densities import Densities


@dataclass(frozen=True)
class VacuumHamiltonian:
    zero_body: float
    one_body: np.ndarray
    two_body: np.ndarray


@dataclass(frozen=True)
class MRHamiltonian:
    zero_body: float
    one_body: np.ndarray
    two_body: np.ndarray


def normal_order(hamiltonian: VacuumHamiltonian, densities: Densities) -> MRHamiltonian:
    """Normal order ``hamiltonian`` with respect to the correlated reference."""
    gamma1 = densities.gamma1
    gamma2 = densities.gamma2
    interaction_contraction = np.einsum(
        "prqs,rs->pq", hamiltonian.two_body, gamma1, optimize=True
    )
    one_body = hamiltonian.one_body + interaction_contraction
    zero_body = (
        hamiltonian.zero_body
        + np.einsum("pq,pq->", hamiltonian.one_body, gamma1, optimize=True)
        + 0.25
        * np.einsum("pqrs,pqrs->", hamiltonian.two_body, gamma2, optimize=True)
    )
    return MRHamiltonian(float(zero_body), one_body, hamiltonian.two_body.copy())


def to_vacuum(hamiltonian: MRHamiltonian, densities: Densities) -> VacuumHamiltonian:
    """Convert MR-normal-ordered 0B/1B/2B pieces back to vacuum ordering."""
    gamma1 = densities.gamma1
    gamma2 = densities.gamma2
    two_body = hamiltonian.two_body.copy()
    one_body = hamiltonian.one_body - np.einsum(
        "prqs,rs->pq", two_body, gamma1, optimize=True
    )
    zero_body = (
        hamiltonian.zero_body
        - np.einsum("pq,pq->", one_body, gamma1, optimize=True)
        - 0.25 * np.einsum("pqrs,pqrs->", two_body, gamma2, optimize=True)
    )
    return VacuumHamiltonian(float(zero_body), one_body, two_body)


def validate_hermitian(
    hamiltonian: VacuumHamiltonian | MRHamiltonian, tolerance: float = 2e-10
) -> None:
    if np.max(np.abs(hamiltonian.one_body - hamiltonian.one_body.T)) > tolerance:
        raise ValueError("one-body tensor is not Hermitian")
    two_body = hamiltonian.two_body
    if np.max(np.abs(two_body + two_body.swapaxes(0, 1))) > tolerance:
        raise ValueError("two-body tensor is not antisymmetric in bra indices")
    if np.max(np.abs(two_body + two_body.swapaxes(2, 3))) > tolerance:
        raise ValueError("two-body tensor is not antisymmetric in ket indices")
    if np.max(np.abs(two_body - two_body.transpose(2, 3, 0, 1))) > tolerance:
        raise ValueError("two-body tensor is not Hermitian")

