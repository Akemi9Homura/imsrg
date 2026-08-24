"""Direct adaptive integration of the rapid MR-IMSRG(2) flow.

The structure follows the existing ``IMSRGSolver`` direct-flow path: update
the generator from the current Hamiltonian, evaluate ``[eta,H]``, use an
adaptive explicit Runge--Kutta method, and stop on the generator norm.  SciPy's
DOP853 is used here so accepted states can be inspected without retaining a
copy of every 40^4 tensor in memory.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from typing import Callable

import numpy as np
from scipy.integrate import DOP853

try:
    from .commutator import commutator
    from .densities import Densities
    from .generator import brillouin_generator, masked_residual_norm
    from .normal_order import MRHamiltonian
except ImportError:
    from commutator import commutator
    from densities import Densities
    from generator import brillouin_generator, masked_residual_norm
    from normal_order import MRHamiltonian


@dataclass(frozen=True)
class FlowSettings:
    smax: float = 2.0
    relative_tolerance: float = 1e-6
    absolute_tolerance: float = 1e-8
    max_step: float = 0.05
    initial_step: float | None = None
    residual_ratio: float = 1e-6
    max_accepted_steps: int = 400
    symmetry_tolerance: float = 2e-8


@dataclass(frozen=True)
class FlowPoint:
    step: int
    s: float
    zero_body: float
    residual: float
    residual_ratio: float
    one_body_hermiticity_error: float
    two_body_hermiticity_error: float
    two_body_antisymmetry_error: float


@dataclass(frozen=True)
class FlowResult:
    hamiltonian: MRHamiltonian
    trajectory: tuple[FlowPoint, ...]
    function_evaluations: int
    converged: bool
    message: str


@lru_cache(maxsize=None)
def _layout(norb: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    one_upper = np.triu_indices(norb)
    pairs = np.asarray(list(combinations(range(norb), 2)), dtype=np.intp)
    pair_upper = np.triu_indices(len(pairs))
    return one_upper[0], one_upper[1], pairs, np.asarray(pair_upper)


def _pack(hamiltonian: MRHamiltonian) -> np.ndarray:
    norb = int(hamiltonian.one_body.shape[0])
    one_i, one_j, pairs, pair_upper = _layout(norb)
    p, q = pairs[:, 0], pairs[:, 1]
    pair_matrix = hamiltonian.two_body[
        p[:, None], q[:, None], p[None, :], q[None, :]
    ]
    return np.concatenate(
        (
            np.asarray([hamiltonian.zero_body], dtype=np.float64),
            np.asarray(hamiltonian.one_body[one_i, one_j], dtype=np.float64),
            np.asarray(
                pair_matrix[pair_upper[0], pair_upper[1]], dtype=np.float64
            ),
        )
    )


def _unpack(values: np.ndarray, norb: int) -> MRHamiltonian:
    one_i, one_j, pairs, pair_upper = _layout(norb)
    one_end = 1 + len(one_i)
    expected = one_end + len(pair_upper[0])
    if values.size != expected:
        raise ValueError(f"packed Hamiltonian has size {values.size}, expected {expected}")
    one_body = np.zeros((norb, norb), dtype=values.dtype)
    one_values = values[1:one_end]
    one_body[one_i, one_j] = one_values
    one_body[one_j, one_i] = one_values

    pair_matrix = np.zeros((len(pairs), len(pairs)), dtype=values.dtype)
    pair_values = values[one_end:]
    pair_matrix[pair_upper[0], pair_upper[1]] = pair_values
    pair_matrix[pair_upper[1], pair_upper[0]] = pair_values
    p, q = pairs[:, 0], pairs[:, 1]
    two_body = np.zeros((norb,) * 4, dtype=values.dtype)
    two_body[p[:, None], q[:, None], p[None, :], q[None, :]] = pair_matrix
    two_body[q[:, None], p[:, None], p[None, :], q[None, :]] = -pair_matrix
    two_body[p[:, None], q[:, None], q[None, :], p[None, :]] = -pair_matrix
    two_body[q[:, None], p[:, None], q[None, :], p[None, :]] = pair_matrix
    return MRHamiltonian(
        float(values[0]),
        one_body,
        two_body,
    )


def _symmetry_errors(hamiltonian: MRHamiltonian) -> tuple[float, float, float]:
    one_error = float(np.max(np.abs(hamiltonian.one_body - hamiltonian.one_body.T)))
    two_body = hamiltonian.two_body
    hermitian_error = float(
        np.max(np.abs(two_body - two_body.transpose(2, 3, 0, 1)))
    )
    antisymmetric_error = max(
        float(np.max(np.abs(two_body + two_body.swapaxes(0, 1)))),
        float(np.max(np.abs(two_body + two_body.swapaxes(2, 3)))),
    )
    return one_error, hermitian_error, antisymmetric_error


def _flow_point(
    step: int,
    s: float,
    hamiltonian: MRHamiltonian,
    residual: float,
    initial_residual: float,
) -> FlowPoint:
    one_error, two_error, antisymmetry_error = _symmetry_errors(hamiltonian)
    ratio = residual / initial_residual if initial_residual > 0.0 else 0.0
    return FlowPoint(
        step=step,
        s=float(s),
        zero_body=float(hamiltonian.zero_body),
        residual=float(residual),
        residual_ratio=float(ratio),
        one_body_hermiticity_error=one_error,
        two_body_hermiticity_error=two_error,
        two_body_antisymmetry_error=antisymmetry_error,
    )


def integrate_flow(
    initial_hamiltonian: MRHamiltonian,
    densities: Densities,
    oscillator_quanta: np.ndarray,
    settings: FlowSettings = FlowSettings(),
    observer: Callable[[FlowPoint], None] | None = None,
) -> FlowResult:
    """Integrate ``dH/ds=[eta(H),H]`` until convergence or ``smax``."""
    if settings.smax <= 0.0:
        raise ValueError("smax must be positive")
    if settings.relative_tolerance <= 0.0 or settings.absolute_tolerance <= 0.0:
        raise ValueError("ODE tolerances must be positive")
    if not 0.0 < settings.residual_ratio < 1.0:
        raise ValueError("residual_ratio must lie between zero and one")

    norb = int(initial_hamiltonian.one_body.shape[0])
    initial_eta = brillouin_generator(
        initial_hamiltonian, densities, oscillator_quanta
    )
    initial_residual = masked_residual_norm(initial_eta)
    trajectory = [
        _flow_point(0, 0.0, initial_hamiltonian, initial_residual, initial_residual)
    ]
    if observer is not None:
        observer(trajectory[0])
    if initial_residual == 0.0:
        return FlowResult(
            hamiltonian=initial_hamiltonian,
            trajectory=tuple(trajectory),
            function_evaluations=0,
            converged=True,
            message="initial Hamiltonian already satisfies the masked decoupling condition",
        )

    function_evaluations = 0

    def right_hand_side(s: float, values: np.ndarray) -> np.ndarray:
        del s
        nonlocal function_evaluations
        hamiltonian = _unpack(values, norb)
        eta = brillouin_generator(hamiltonian, densities, oscillator_quanta)
        derivative = commutator(eta, hamiltonian, densities)
        packed = _pack(derivative)
        if not np.all(np.isfinite(packed)):
            raise FloatingPointError("non-finite MR-IMSRG derivative")
        function_evaluations += 1
        return packed

    solver = DOP853(
        right_hand_side,
        0.0,
        _pack(initial_hamiltonian),
        settings.smax,
        rtol=settings.relative_tolerance,
        atol=settings.absolute_tolerance,
        max_step=settings.max_step,
        first_step=settings.initial_step,
    )

    converged = False
    message = "maximum flow parameter reached before residual target"
    for step in range(1, settings.max_accepted_steps + 1):
        if solver.status == "finished":
            break
        previous_s = solver.t
        solver.step()
        if solver.status == "failed":
            message = "ODE solver failed"
            break
        if solver.t == previous_s:
            message = "ODE solver made no progress"
            break
        hamiltonian = _unpack(solver.y, norb)
        eta = brillouin_generator(hamiltonian, densities, oscillator_quanta)
        residual = masked_residual_norm(eta)
        point = _flow_point(step, solver.t, hamiltonian, residual, initial_residual)
        trajectory.append(point)
        if observer is not None:
            observer(point)
        if max(
            point.one_body_hermiticity_error,
            point.two_body_hermiticity_error,
            point.two_body_antisymmetry_error,
        ) > settings.symmetry_tolerance:
            message = "Hamiltonian symmetry error exceeded the configured tolerance"
            break
        if point.residual_ratio <= settings.residual_ratio:
            converged = True
            message = "masked Brillouin residual target reached"
            break
        if solver.status == "finished":
            break
    else:
        message = "maximum number of accepted ODE steps reached"

    final_hamiltonian = _unpack(solver.y.copy(), norb)
    return FlowResult(
        hamiltonian=final_hamiltonian,
        trajectory=tuple(trajectory),
        function_evaluations=function_evaluations,
        converged=converged,
        message=message,
    )
