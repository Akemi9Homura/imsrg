"""Independent m-scheme oracle for scalar MR-Magnus(2).

The production C++ path accumulates a finite Magnus operator with
``exp(Omega_new) = exp(dOmega) exp(Omega_old)`` and reconstructs operators as
``exp(Omega) O exp(-Omega)``.  This module mirrors those two algebraic steps
using the explicit :mod:`commutator` implementation.  It is deliberately not
an ODE solver and is never used as a production runtime backend.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

try:
    from .commutator import commutator
    from .densities import Densities
    from .normal_order import MRHamiltonian
except ImportError:  # Direct use from the prototype directory.
    from commutator import commutator
    from densities import Densities
    from normal_order import MRHamiltonian


BERNOULLI = (1.0, -0.5, 1.0 / 6.0, 0.0, -1.0 / 30.0, 0.0, 1.0 / 42.0, 0.0, -1.0 / 30.0)
FACTORIAL = (1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0)


def scale(operator: MRHamiltonian, coefficient: float) -> MRHamiltonian:
    return MRHamiltonian(
        coefficient * operator.zero_body,
        coefficient * operator.one_body,
        coefficient * operator.two_body,
    )


def add(*operators: MRHamiltonian) -> MRHamiltonian:
    if not operators:
        raise ValueError("at least one operator is required")
    return MRHamiltonian(
        sum(operator.zero_body for operator in operators),
        sum((operator.one_body for operator in operators), np.zeros_like(operators[0].one_body)),
        sum((operator.two_body for operator in operators), np.zeros_like(operators[0].two_body)),
    )


def operator_norm(operator: MRHamiltonian) -> float:
    """Return the explicit coefficient two-norm used for oracle convergence."""
    return float(
        np.sqrt(
            abs(operator.zero_body) ** 2
            + np.vdot(operator.one_body, operator.one_body).real
            + np.vdot(operator.two_body, operator.two_body).real
        )
    )


def max_abs(operator: MRHamiltonian) -> float:
    return max(
        abs(operator.zero_body),
        float(np.max(np.abs(operator.one_body))),
        float(np.max(np.abs(operator.two_body))),
    )


@dataclass(frozen=True)
class BCHTransformResult:
    operator: MRHamiltonian
    # terms[k] = ad_Omega^k(initial) / k!
    terms: tuple[MRHamiltonian, ...]


def bch_transform(
    initial: MRHamiltonian,
    omega: MRHamiltonian,
    densities: Densities,
    *,
    max_order: int = 40,
    relative_threshold: float | None = None,
) -> BCHTransformResult:
    """Return ``exp(Omega) initial exp(-Omega)`` in MR-Magnus(2).

    When ``relative_threshold`` is supplied, the iterative BCH construction
    stops once the newest factorial-scaled term is small relative to the
    accumulated operator, matching Vobig Eq. (4.6.11).  Tests that compare
    individual orders leave it unset and specify an exact ``max_order``.
    """
    if max_order < 0:
        raise ValueError("max_order must be non-negative")
    if relative_threshold is not None and relative_threshold <= 0.0:
        raise ValueError("relative_threshold must be positive")

    terms = [initial]
    transformed = initial
    for order in range(1, max_order + 1):
        # If terms[k-1] = ad^(k-1)/(k-1)!, then this is ad^k/k!.
        term = scale(commutator(omega, terms[-1], densities), 1.0 / order)
        terms.append(term)
        transformed = add(transformed, term)
        if relative_threshold is not None:
            denominator = operator_norm(transformed)
            if denominator == 0.0:
                if operator_norm(term) == 0.0:
                    break
            elif operator_norm(term) / denominator < relative_threshold:
                break
    return BCHTransformResult(transformed, tuple(terms))


@dataclass(frozen=True)
class BCHProductResult:
    omega: MRHamiltonian
    # Named additive pieces in the precise order used by imsrg++ BCH_Product.
    contributions: tuple[tuple[str, MRHamiltonian], ...]


def bch_product(
    d_omega: MRHamiltonian,
    omega: MRHamiltonian,
    densities: Densities,
    *,
    threshold: float = 1e-4,
) -> BCHProductResult:
    """Mirror ``imsrg++``'s production ``BCH_Product(dOmega, Omega)``.

    The result satisfies ``exp(result) = exp(dOmega) exp(Omega)`` up to the
    IMSRG(2) and finite BCH-product truncations.  In particular the first
    nested object is ``[Omega,dOmega]`` and its Bernoulli-1 coefficient is
    ``-1/2``, yielding the conventional ``+[dOmega,Omega]/2`` term.
    """
    if threshold < 0.0:
        raise ValueError("threshold must be non-negative")

    contributions: list[tuple[str, MRHamiltonian]] = [
        ("domega", d_omega),
        ("omega", omega),
    ]
    result = add(d_omega, omega)
    nested = commutator(omega, d_omega, densities)
    nested_norm = operator_norm(nested)

    # Existing imsrg++ separately retains 1/12 [dOmega,[dOmega,Omega]],
    # written there as 1/12 [[Omega,dOmega],dOmega].
    if nested_norm * operator_norm(d_omega) > threshold:
        x_nested = scale(commutator(nested, d_omega, densities), 1.0 / 12.0)
        result = add(result, x_nested)
        contributions.append(("domega_domega_omega", x_nested))

    order = 1
    while nested_norm > threshold:
        if order < 2 or order % 2 == 0:
            term = scale(nested, BERNOULLI[order] / FACTORIAL[order])
            result = add(result, term)
            contributions.append((f"bernoulli_{order}", term))

        order += 1
        if order >= len(BERNOULLI):
            break
        if 2.0 * operator_norm(omega) * nested_norm < threshold:
            break
        nested = commutator(omega, nested, densities)
        nested_norm = operator_norm(nested)

    return BCHProductResult(result, tuple(contributions))


def magnus_derivative(
    omega: MRHamiltonian,
    eta: MRHamiltonian,
    densities: Densities,
    *,
    relative_threshold: float = 1e-2,
) -> MRHamiltonian:
    """Return the production Bernoulli series for ``dOmega/ds``.

    Set ``relative_threshold=0`` to retain every term through k=8.  Positive
    values implement Vobig Eq. (4.6.7) together with the monotonic nested-norm
    guard of Eq. (4.6.8).
    """
    if relative_threshold < 0.0:
        raise ValueError("relative_threshold must be non-negative")
    derivative = eta
    nested = eta
    previous_nested_norm = operator_norm(nested)
    monotonically_decreasing = True
    for order in range(1, len(BERNOULLI)):
        nested = commutator(omega, nested, densities)
        nested_norm = operator_norm(nested)
        if previous_nested_norm > 0.0 and nested_norm >= previous_nested_norm:
            monotonically_decreasing = False
        if BERNOULLI[order] != 0.0:
            term = scale(nested, BERNOULLI[order] / FACTORIAL[order])
            derivative = add(derivative, term)
            if (
                relative_threshold > 0.0
                and order >= 2
                and monotonically_decreasing
                and operator_norm(term)
                < relative_threshold * operator_norm(derivative)
            ):
                return derivative
        previous_nested_norm = nested_norm
    return derivative
