"""Direct adaptive integration of the rapid MR-IMSRG(2) flow.

The structure follows the existing ``IMSRGSolver`` direct-flow path: update
the generator from the current Hamiltonian, evaluate ``[eta,H]``, use an
adaptive explicit Runge--Kutta method, and stop on the masked decoupling
residual. SciPy's DOP853 is used here so accepted states can be inspected
without retaining a copy of every 40^4 tensor in memory.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from typing import Callable

import numpy as np
from scipy.integrate import DOP853

try:
    from .basis import prepare_natural_basis, transform_hamiltonian
    from .commutator import commutator
    from .densities import Densities
    from .generator import (
        decoupling_masks,
        decoupling_residual,
        masked_residual_norm,
        white_generator_unmasked,
        white_ncsm_numerator_residual,
    )
    from .normal_order import MRHamiltonian
except ImportError:
    from basis import prepare_natural_basis, transform_hamiltonian
    from commutator import commutator
    from densities import Densities
    from generator import (
        decoupling_masks,
        decoupling_residual,
        masked_residual_norm,
        white_generator_unmasked,
        white_ncsm_numerator_residual,
    )
    from normal_order import MRHamiltonian


@dataclass(frozen=True)
class FlowSettings:
    smax: float = 2.0
    relative_tolerance: float = 1e-6
    absolute_tolerance: float = 1e-8
    # DOP853 still selects smaller steps from the error estimate; this only
    # prevents an unnecessary fixed ceiling on the smooth late-time flow.
    max_step: float = 10.0
    initial_step: float | None = None
    residual_ratio: float = 1e-6
    max_accepted_steps: int = 400
    symmetry_tolerance: float = 2e-8
    checkpoint_s: float | None = None


@dataclass(frozen=True)
class FlowPoint:
    step: int
    s: float
    zero_body: float
    residual: float
    residual_ratio: float
    generator_numerator_residual: float
    generator_numerator_residual_ratio: float
    one_body_hermiticity_error: float
    two_body_hermiticity_error: float
    two_body_antisymmetry_error: float


@dataclass(frozen=True)
class FlowCheckpoint:
    point: FlowPoint
    hamiltonian: MRHamiltonian


@dataclass(frozen=True)
class FlowResult:
    hamiltonian: MRHamiltonian
    trajectory: tuple[FlowPoint, ...]
    function_evaluations: int
    converged: bool
    message: str
    checkpoints: tuple[FlowCheckpoint, ...] = ()


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
    generator_residual: float,
    initial_generator_residual: float,
) -> FlowPoint:
    one_error, two_error, antisymmetry_error = _symmetry_errors(hamiltonian)
    ratio = residual / initial_residual if initial_residual > 0.0 else 0.0
    generator_ratio = (
        generator_residual / initial_generator_residual
        if initial_generator_residual > 0.0
        else 0.0
    )
    return FlowPoint(
        step=step,
        s=float(s),
        zero_body=float(hamiltonian.zero_body),
        residual=float(residual),
        residual_ratio=float(ratio),
        generator_numerator_residual=float(generator_residual),
        generator_numerator_residual_ratio=float(generator_ratio),
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
    if settings.checkpoint_s is not None and not 0.0 < settings.checkpoint_s < settings.smax:
        raise ValueError("checkpoint_s must lie strictly between zero and smax")

    norb = int(initial_hamiltonian.one_body.shape[0])
    mask1, mask2 = decoupling_masks(oscillator_quanta)
    if mask1.shape != initial_hamiltonian.one_body.shape:
        raise ValueError("oscillator-quanta array has an incompatible length")

    natural_basis = prepare_natural_basis(densities)
    working_densities = natural_basis.densities
    working_initial = (
        initial_hamiltonian
        if natural_basis.is_identity
        else transform_hamiltonian(
            initial_hamiltonian, natural_basis.vectors, to_natural=True
        )
    )

    def apply_ho_mask(operator: MRHamiltonian) -> MRHamiltonian:
        original_basis_operator = (
            operator
            if natural_basis.is_identity
            else transform_hamiltonian(
                operator, natural_basis.vectors, to_natural=False
            )
        )
        masked_original = MRHamiltonian(
            0.0,
            original_basis_operator.one_body * mask1,
            original_basis_operator.two_body * mask2,
        )
        return (
            masked_original
            if natural_basis.is_identity
            else transform_hamiltonian(
                masked_original, natural_basis.vectors, to_natural=True
            )
        )

    initial_residual_operator = apply_ho_mask(
        decoupling_residual(working_initial, working_densities)
    )
    initial_residual = masked_residual_norm(initial_residual_operator)
    initial_generator_residual = masked_residual_norm(
        apply_ho_mask(
            white_ncsm_numerator_residual(working_initial, working_densities)
        )
    )
    trajectory = [
        _flow_point(
            0,
            0.0,
            working_initial,
            initial_residual,
            initial_residual,
            initial_generator_residual,
            initial_generator_residual,
        )
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
        eta = apply_ho_mask(
            white_generator_unmasked(hamiltonian, working_densities)
        )
        derivative = commutator(eta, hamiltonian, working_densities)
        packed = _pack(derivative)
        if not np.all(np.isfinite(packed)):
            raise FloatingPointError("non-finite MR-IMSRG derivative")
        function_evaluations += 1
        return packed

    solver = DOP853(
        right_hand_side,
        0.0,
        _pack(working_initial),
        settings.smax,
        rtol=settings.relative_tolerance,
        atol=settings.absolute_tolerance,
        max_step=settings.max_step,
        first_step=settings.initial_step,
    )

    converged = False
    message = "maximum flow parameter reached before residual target"
    checkpoints: list[FlowCheckpoint] = []
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
        if (
            settings.checkpoint_s is not None
            and not checkpoints
            and previous_s < settings.checkpoint_s <= solver.t
        ):
            checkpoint_hamiltonian = _unpack(
                solver.dense_output()(settings.checkpoint_s), norb
            )
            checkpoint_residual_operator = apply_ho_mask(
                decoupling_residual(checkpoint_hamiltonian, working_densities)
            )
            checkpoint_residual = masked_residual_norm(checkpoint_residual_operator)
            checkpoint_generator_residual = masked_residual_norm(
                apply_ho_mask(
                    white_ncsm_numerator_residual(
                        checkpoint_hamiltonian, working_densities
                    )
                )
            )
            checkpoints.append(
                FlowCheckpoint(
                    point=_flow_point(
                        step,
                        settings.checkpoint_s,
                        checkpoint_hamiltonian,
                        checkpoint_residual,
                        initial_residual,
                        checkpoint_generator_residual,
                        initial_generator_residual,
                    ),
                    hamiltonian=checkpoint_hamiltonian,
                )
            )
        hamiltonian = _unpack(solver.y, norb)
        residual_operator = apply_ho_mask(
            decoupling_residual(hamiltonian, working_densities)
        )
        residual = masked_residual_norm(residual_operator)
        generator_residual = masked_residual_norm(
            apply_ho_mask(
                white_ncsm_numerator_residual(hamiltonian, working_densities)
            )
        )
        point = _flow_point(
            step,
            solver.t,
            hamiltonian,
            residual,
            initial_residual,
            generator_residual,
            initial_generator_residual,
        )
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
            message = "masked decoupling residual target reached"
            break
        if solver.status == "finished":
            break
    else:
        message = "maximum number of accepted ODE steps reached"

    working_final = _unpack(solver.y.copy(), norb)
    final_hamiltonian = (
        working_final
        if natural_basis.is_identity
        else transform_hamiltonian(
            working_final, natural_basis.vectors, to_natural=False
        )
    )
    if not natural_basis.is_identity:
        checkpoints = [
            FlowCheckpoint(
                checkpoint.point,
                transform_hamiltonian(
                    checkpoint.hamiltonian,
                    natural_basis.vectors,
                    to_natural=False,
                ),
            )
            for checkpoint in checkpoints
        ]
    return FlowResult(
        hamiltonian=final_hamiltonian,
        trajectory=tuple(trajectory),
        function_evaluations=function_evaluations,
        converged=converged,
        message=message,
        checkpoints=tuple(checkpoints),
    )
