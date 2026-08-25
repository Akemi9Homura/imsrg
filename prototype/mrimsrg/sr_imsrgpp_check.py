"""Strict single-reference degeneration check against the actual imsrg++ code.

The input is the float64 ``jcoupled64`` diagnostic generated from the same
ordinary N2LOopt Hamiltonian used by the m-scheme prototype.  It is loaded into
an ``imsrg++`` ``Operator``, normal ordered by ``Operator::DoNormalOrdering``,
and passed through ``Generator::Update`` and ``Commutator::Commutator``.  The
resulting J-coupled tensors are expanded to the prototype's m-scheme ordering
and compared element by element with the *production* MR path.

This file deliberately contains no alternative SR generator or commutator.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import importlib
import json
from pathlib import Path
import struct
import sys
from typing import Any

import numpy as np
from sympy import Rational
from sympy.physics.wigner import clebsch_gordan

try:
    from .commutator import commutator, commutator_contributions
    from .densities import compute_densities, validate_densities
    from .generator import (
        GENERATOR_IMPLEMENTATION,
        WHITE_DENOMINATOR_CUTOFF,
        _safe_denominator,
        decoupling_masks,
        epstein_nesbet_denominators,
        oscillator_quanta_from_orbits,
        spherical_orbit_groups_from_orbits,
        white_generator,
    )
    from .normal_order import MRHamiltonian, VacuumHamiltonian, normal_order, to_vacuum
    from .reference_io import load_reference
except ImportError:
    from commutator import commutator, commutator_contributions
    from densities import compute_densities, validate_densities
    from generator import (
        GENERATOR_IMPLEMENTATION,
        WHITE_DENOMINATOR_CUTOFF,
        _safe_denominator,
        decoupling_masks,
        epstein_nesbet_denominators,
        oscillator_quanta_from_orbits,
        spherical_orbit_groups_from_orbits,
        white_generator,
    )
    from normal_order import MRHamiltonian, VacuumHamiltonian, normal_order, to_vacuum
    from reference_io import load_reference


J64_MAGIC = b"mrimsrg_j64_v1\0\0"
SOURCE_ORACLE_FILES = ("Generator.cc", "Commutator.cc", "IMSRGSolver.cc")
NUCLEI = {"He4": (4, 2), "O16": (16, 8)}


@dataclass(frozen=True)
class JCoupled64:
    hw: float
    emax: int
    orbits: np.ndarray
    zero_body: float
    one_body: np.ndarray
    records: tuple[tuple[int, int, int, int, int, float], ...]


def _read_exact(stream: Any, size: int, field: str) -> bytes:
    value = stream.read(size)
    if len(value) != size:
        raise ValueError(f"truncated jcoupled64 while reading {field}")
    return value


def read_jcoupled64(path: Path) -> JCoupled64:
    with path.open("rb") as stream:
        if _read_exact(stream, 16, "magic") != J64_MAGIC:
            raise ValueError("unsupported jcoupled64 payload")
        hw = struct.unpack("<d", _read_exact(stream, 8, "hw"))[0]
        emax = struct.unpack("<i", _read_exact(stream, 4, "emax"))[0]
        norb, nobme, ntbme = struct.unpack(
            "<QQQ", _read_exact(stream, 24, "counts")
        )
        orbits = np.asarray(
            [
                struct.unpack("<iiii", _read_exact(stream, 16, "orbit"))
                for _ in range(norb)
            ],
            dtype=np.int64,
        )
        zero_body = struct.unpack("<d", _read_exact(stream, 8, "zero body"))[0]
        if nobme != norb * (norb + 1) // 2:
            raise ValueError("invalid jcoupled64 one-body count")
        one_body = np.zeros((norb, norb), dtype=np.float64)
        for a in range(norb):
            for b in range(a, norb):
                value = struct.unpack("<d", _read_exact(stream, 8, "OBME"))[0]
                one_body[a, b] = value
                one_body[b, a] = value
        records = []
        for _ in range(ntbme):
            a, b, c, d, j = struct.unpack(
                "<iiiii", _read_exact(stream, 20, "TBME indices")
            )
            value = struct.unpack("<d", _read_exact(stream, 8, "TBME"))[0]
            records.append((a, b, c, d, j, value))
        if stream.read(1):
            raise ValueError("jcoupled64 payload has trailing bytes")
    return JCoupled64(hw, emax, orbits, zero_body, one_body, tuple(records))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_oracle_sources(repository: Path, oracle_source: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for name in SOURCE_ORACLE_FILES:
        local = repository / "src" / name
        oracle = oracle_source / "src" / name
        if not local.is_file() or not oracle.is_file():
            raise ValueError(f"missing source file needed for oracle check: {name}")
        local_hash = _sha256(local)
        oracle_hash = _sha256(oracle)
        if local_hash != oracle_hash:
            raise ValueError(
                f"oracle {name} differs from the current repository: "
                f"{oracle_hash} != {local_hash}"
            )
        result[name] = local_hash
    return result


def import_pyimsrg(module_dir: Path) -> Any:
    sys.path.insert(0, str(module_dir.resolve()))
    try:
        return importlib.import_module("pyIMSRG")
    finally:
        sys.path.pop(0)


def _modelspace_orbits(modelspace: Any) -> np.ndarray:
    return np.asarray(
        [
            (
                modelspace.GetOrbit(i).n,
                modelspace.GetOrbit(i).l,
                modelspace.GetOrbit(i).j2,
                modelspace.GetOrbit(i).tz2,
            )
            for i in range(modelspace.GetNumberOrbits())
        ],
        dtype=np.int64,
    )


def make_imsrg_operator(pyimsrg: Any, payload: JCoupled64, reference: str) -> tuple[Any, Any]:
    modelspace = pyimsrg.ModelSpace(payload.emax, reference)
    modelspace.SetHbarOmega(payload.hw)
    if not np.array_equal(_modelspace_orbits(modelspace), payload.orbits):
        raise ValueError("imsrg++ and jcoupled64 spherical-orbit tables differ")

    operator = pyimsrg.Operator(modelspace, 0, 0, 0, 2)
    operator.SetHermitian()
    operator.ZeroBody = payload.zero_body
    for a in range(payload.one_body.shape[0]):
        for b in range(payload.one_body.shape[1]):
            operator.SetOneBody(a, b, float(payload.one_body[a, b]))
    for a, b, c, d, j, value in payload.records:
        oa = modelspace.GetOrbit(a)
        ob = modelspace.GetOrbit(b)
        channel_index = modelspace.GetTwoBodyChannelIndex(
            j, (oa.l + ob.l) % 2, (oa.tz2 + ob.tz2) // 2
        )
        channel = modelspace.GetTwoBodyChannel(channel_index)
        bra = channel.GetLocalIndex(a, b)
        ket = channel.GetLocalIndex(c, d)
        if bra < 0 or ket < 0:
            raise ValueError("jcoupled64 TBME is absent from the imsrg++ channel")
        operator.TwoBody.SetTBME_chij(
            channel_index, channel_index, bra, ket, float(value)
        )
    return modelspace, operator


@lru_cache(maxsize=None)
def _cg(j1_2: int, j2_2: int, j_2: int, m1_2: int, m2_2: int) -> float:
    m_2 = m1_2 + m2_2
    return float(
        clebsch_gordan(
            Rational(j1_2, 2),
            Rational(j2_2, 2),
            Rational(j_2, 2),
            Rational(m1_2, 2),
            Rational(m2_2, 2),
            Rational(m_2, 2),
        )
    )


def _phase(exponent: int) -> float:
    return -1.0 if exponent % 2 else 1.0


def operator_to_mscheme(operator: Any, orbits: np.ndarray) -> MRHamiltonian:
    """Expand a scalar imsrg++ operator using the validated bridge convention."""
    norb = orbits.shape[0]
    one_body = np.zeros((norb, norb), dtype=np.float64)
    for p in range(norb):
        a, _, lp, jp, mp, tzp = (int(x) for x in orbits[p])
        for q in range(norb):
            b, _, lq, jq, mq, tzq = (int(x) for x in orbits[q])
            if (lp, jp, mp, tzp) == (lq, jq, mq, tzq):
                one_body[p, q] = operator.GetOneBody(a, b)

    two_body = np.zeros((norb, norb, norb, norb), dtype=np.float64)
    for p in range(norb):
        ap, _, lp, jp, mp, tzp = (int(x) for x in orbits[p])
        for q in range(norb):
            if p == q:
                continue
            aq, _, lq, jq, mq, tzq = (int(x) for x in orbits[q])
            for r in range(norb):
                ar, _, lr, jr, mr, tzr = (int(x) for x in orbits[r])
                for s in range(norb):
                    if r == s:
                        continue
                    ass, _, ls, js, ms, tzs = (int(x) for x in orbits[s])
                    if mp + mq != mr + ms or tzp + tzq != tzr + tzs:
                        continue
                    if (lp + lq - lr - ls) % 2:
                        continue

                    a, b, c, d = ap, aq, ar, ass
                    swap_ab = a > b
                    swap_cd = c > d
                    if swap_ab:
                        a, b = b, a
                    if swap_cd:
                        c, d = d, c
                    jmin = max(abs(jp - jq), abs(jr - js)) // 2
                    jmax = min(jp + jq, jr + js) // 2
                    value = 0.0
                    for j in range(jmin, jmax + 1):
                        if (a == b or c == d) and j % 2:
                            continue
                        cg_bra = _cg(jp, jq, 2 * j, mp, mq)
                        cg_ket = _cg(jr, js, 2 * j, mr, ms)
                        if cg_bra == 0.0 or cg_ket == 0.0:
                            continue
                        pair_phase = 1.0
                        if a == b:
                            pair_phase *= np.sqrt(2.0)
                        if c == d:
                            pair_phase *= np.sqrt(2.0)
                        if swap_ab:
                            pair_phase *= _phase((jp + jq) // 2 - j + 1)
                        if swap_cd:
                            pair_phase *= _phase((jr + js) // 2 - j + 1)
                        value += (
                            pair_phase
                            * operator.TwoBody.GetTBME_J_norm(j, j, a, b, c, d)
                            * cg_bra
                            * cg_ket
                        )
                    two_body[p, q, r, s] = value
    return MRHamiltonian(float(operator.ZeroBody), one_body, two_body)


def _comparison(actual: np.ndarray, expected: np.ndarray) -> dict[str, Any]:
    difference = np.asarray(actual) - np.asarray(expected)
    flat_index = int(np.argmax(np.abs(difference)))
    worst_index = tuple(int(x) for x in np.unravel_index(flat_index, difference.shape))
    expected_norm = float(np.linalg.norm(np.asarray(expected).ravel()))
    difference_norm = float(np.linalg.norm(difference.ravel()))
    return {
        "max_abs": float(np.max(np.abs(difference))),
        "relative_frobenius": (
            difference_norm / expected_norm if expected_norm > 1e-14 else None
        ),
        "worst_index": worst_index,
        "actual_at_worst": float(np.asarray(actual)[worst_index]),
        "expected_at_worst": float(np.asarray(expected)[worst_index]),
    }


def _operator_comparison(actual: MRHamiltonian, expected: MRHamiltonian) -> dict[str, Any]:
    zero_difference = float(actual.zero_body - expected.zero_body)
    return {
        "zero_body": {
            "actual": actual.zero_body,
            "expected": expected.zero_body,
            "difference": zero_difference,
            "max_abs": abs(zero_difference),
        },
        "one_body": _comparison(actual.one_body, expected.one_body),
        "two_body": _comparison(actual.two_body, expected.two_body),
    }


def _commutator_contribution_report(
    pyimsrg: Any,
    modelspace: Any,
    left_j: Any,
    right_j: Any,
    production_terms: dict[str, MRHamiltonian],
    orbits: np.ndarray,
) -> dict[str, Any]:
    """Compare the named production contractions to current C++ routines."""
    report: dict[str, Any] = {}
    cpp_names = (
        "comm110ss",
        "comm220ss",
        "comm111ss",
        "comm121ss",
        "comm221ss",
        "comm122ss",
        "comm222_pp_hhss",
        "comm222_phss",
    )
    for name in cpp_names:
        result_j = pyimsrg.Operator(modelspace, 0, 0, 0, 2)
        result_j.SetHermitian()
        getattr(pyimsrg.Commutator, name)(left_j, right_j, result_j)
        report[name] = _operator_comparison(
            operator_to_mscheme(result_j, orbits), production_terms[name]
        )
    for name in ("mr_lambda2_one_body", "mr_lambda2_zero_body"):
        zero = MRHamiltonian(
            0.0,
            np.zeros_like(production_terms[name].one_body),
            np.zeros_like(production_terms[name].two_body),
        )
        report[name] = _operator_comparison(zero, production_terms[name])
    return report


def _mask_report(densities: Any, oscillator_quanta: np.ndarray) -> dict[str, int]:
    occupations = np.diag(densities.gamma1)
    holes = occupations > 0.5
    particles = ~holes
    sr_one = (particles[:, None] & holes[None, :]) | (
        holes[:, None] & particles[None, :]
    )
    pair_holes = holes[:, None] & holes[None, :]
    pair_particles = particles[:, None] & particles[None, :]
    sr_two = (
        pair_particles[:, :, None, None] & pair_holes[None, None, :, :]
    ) | (pair_holes[:, :, None, None] & pair_particles[None, None, :, :])
    relaxed_one, relaxed_two = decoupling_masks(oscillator_quanta)
    return {
        "sr_one_body_elements": int(np.count_nonzero(sr_one)),
        "sr_two_body_elements": int(np.count_nonzero(sr_two)),
        "sr_one_body_removed_by_delta_e": int(np.count_nonzero(sr_one & ~relaxed_one)),
        "sr_two_body_removed_by_delta_e": int(np.count_nonzero(sr_two & ~relaxed_two)),
    }


def _mr_add_scaled(
    base: MRHamiltonian, increment: MRHamiltonian, coefficient: float
) -> MRHamiltonian:
    return MRHamiltonian(
        base.zero_body + coefficient * increment.zero_body,
        base.one_body + coefficient * increment.one_body,
        base.two_body + coefficient * increment.two_body,
    )


def _mr_rk4_step(
    hamiltonian: MRHamiltonian, step: float, rhs: Any
) -> MRHamiltonian:
    k1 = rhs(hamiltonian)
    k2 = rhs(_mr_add_scaled(hamiltonian, k1, 0.5 * step))
    k3 = rhs(_mr_add_scaled(hamiltonian, k2, 0.5 * step))
    k4 = rhs(_mr_add_scaled(hamiltonian, k3, step))
    return MRHamiltonian(
        hamiltonian.zero_body
        + step * (k1.zero_body + 2 * k2.zero_body + 2 * k3.zero_body + k4.zero_body) / 6,
        hamiltonian.one_body
        + step * (k1.one_body + 2 * k2.one_body + 2 * k3.one_body + k4.one_body) / 6,
        hamiltonian.two_body
        + step * (k1.two_body + 2 * k2.two_body + 2 * k3.two_body + k4.two_body) / 6,
    )


def _imsrg_rk4_step(hamiltonian: Any, step: float, rhs: Any) -> Any:
    k1 = rhs(hamiltonian)
    k2 = rhs(hamiltonian + (0.5 * step) * k1)
    k3 = rhs(hamiltonian + (0.5 * step) * k2)
    k4 = rhs(hamiltonian + step * k3)
    return hamiltonian + (step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def _denominator_report(
    generator: Any,
    modelspace: Any,
    production_h: MRHamiltonian,
    densities: Any,
    orbits: np.ndarray,
) -> dict[str, Any]:
    groups = spherical_orbit_groups_from_orbits(orbits)
    delta1, delta2 = epstein_nesbet_denominators(
        production_h, densities, spherical_orbit_groups=groups
    )
    representatives = {
        int(group): int(np.flatnonzero(groups == group)[0]) for group in np.unique(groups)
    }

    one_records = []
    for hole in sorted(modelspace.core):
        for particle in sorted(set(modelspace.valence) | set(modelspace.qspace)):
            raw = float(delta1[representatives[particle], representatives[hole]])
            production = float(
                _safe_denominator(np.asarray(raw), WHITE_DENOMINATOR_CUTOFF)
            )
            oracle = float(generator.Get1bDenominator(particle, hole))
            one_records.append((particle, hole, raw, production, oracle))

    two_records = []
    for channel_index in range(modelspace.GetNumberTwoBodyChannels()):
        channel = modelspace.GetTwoBodyChannel(channel_index)
        holes = channel.GetKetIndex_cc()
        particles = sorted(
            set(channel.GetKetIndex_qq())
            | set(channel.GetKetIndex_vv())
            | set(channel.GetKetIndex_qv())
        )
        for ket_index in holes:
            ket = channel.GetKet(ket_index)
            for bra_index in particles:
                bra = channel.GetKet(bra_index)
                raw = float(
                    delta2[
                        representatives[bra.p],
                        representatives[bra.q],
                        representatives[ket.p],
                        representatives[ket.q],
                    ]
                )
                production = float(
                    _safe_denominator(np.asarray(raw), WHITE_DENOMINATOR_CUTOFF)
                )
                oracle = float(
                    generator.Get2bDenominator(channel_index, bra_index, ket_index)
                )
                two_records.append(
                    (
                        channel_index,
                        bra.p,
                        bra.q,
                        ket.p,
                        ket.q,
                        raw,
                        production,
                        oracle,
                    )
                )

    def summarize(records: list[tuple[Any, ...]], value_offset: int) -> dict[str, Any]:
        differences = np.asarray(
            [record[value_offset + 1] - record[value_offset + 2] for record in records]
        )
        worst = int(np.argmax(np.abs(differences)))
        record = records[worst]
        raw_values = np.asarray([record[value_offset] for record in records])
        return {
            "count": len(records),
            "cutoff_trigger_count": int(
                np.count_nonzero(np.abs(raw_values) < WHITE_DENOMINATOR_CUTOFF)
            ),
            "max_abs": float(np.max(np.abs(differences))),
            "worst_labels": tuple(int(x) for x in record[:value_offset]),
            "raw_at_worst": float(record[value_offset]),
            "production_at_worst": float(record[value_offset + 1]),
            "imsrgpp_at_worst": float(record[value_offset + 2]),
        }

    return {
        "cutoff_mev": WHITE_DENOMINATOR_CUTOFF,
        "one_body": summarize(one_records, 2),
        "two_body": summarize(two_records, 5),
    }


def _load_production_flow(
    path: Path, nucleus: str, reference: Any
) -> tuple[MRHamiltonian, float, dict[str, Any]]:
    metadata_path = path / "metadata.json"
    if not metadata_path.is_file():
        raise ValueError(f"production flow metadata is unavailable: {metadata_path}")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("schema") != "mrimsrg_flow_v1":
        raise ValueError("production flow has an unsupported schema")
    if metadata.get("generator_implementation") != GENERATOR_IMPLEMENTATION:
        raise ValueError("production flow was not made by the current generator")
    denominator = metadata.get("generator_denominator", {})
    if denominator.get("cutoff_behavior") != (
        "absolute values below cutoff replaced by positive cutoff"
    ):
        raise ValueError("production flow does not record the current cutoff behavior")
    reference_metadata = metadata.get("reference_metadata", {})
    for key in ("A", "Z", "Nrefmax"):
        if reference_metadata.get(key) != reference.metadata.get(key):
            raise ValueError(f"production flow reference differs in {key}")
    if metadata.get("formula_basis", {}).get("equation_evaluation") != (
        "original basis already diagonal in gamma1"
    ):
        raise ValueError("strict SR flow comparison requires the identity natural basis")
    trajectory = metadata.get("trajectory")
    if not isinstance(trajectory, list) or not trajectory:
        raise ValueError("production flow has no trajectory")
    target_s = float(trajectory[-1]["s"])
    settings = metadata.get("flow_settings", {})
    if not np.isclose(target_s, float(settings.get("smax", np.nan)), rtol=0, atol=1e-10):
        raise ValueError("production flow stopped before the fixed comparison s")
    if metadata.get("flow_converged"):
        raise ValueError("production flow used an early residual stop")
    one_body = np.load(path / "final_mr_one_body.npy", allow_pickle=False)
    two_body = np.load(path / "final_mr_two_body.npy", allow_pickle=False)
    hamiltonian = MRHamiltonian(
        float(metadata["final_mr_zero_body"]), one_body, two_body
    )
    expected_shape = reference.one_body.shape
    if one_body.shape != expected_shape or two_body.shape != reference.two_body.shape:
        raise ValueError("production flow tensor shapes differ from the reference")
    return hamiltonian, target_s, metadata


def _full_flow_report(
    pyimsrg: Any,
    modelspace: Any,
    initial_h_j: Any,
    production_flow: Path,
    production_h_initial: MRHamiltonian,
    densities: Any,
    quanta: np.ndarray,
    groups: np.ndarray,
    reference: Any,
    nucleus: str,
    ode_tolerance: float,
    initial_step: float,
) -> dict[str, Any]:
    production_final, target_s, metadata = _load_production_flow(
        production_flow, nucleus, reference
    )
    settings = metadata["flow_settings"]
    if not np.isclose(float(settings["relative_tolerance"]), ode_tolerance):
        raise ValueError("production rtol differs from the requested C++ tolerance")
    if not np.isclose(float(settings["absolute_tolerance"]), ode_tolerance):
        raise ValueError("production atol differs from the requested C++ tolerance")

    solver = pyimsrg.IMSRGSolver(initial_h_j)
    solver.SetMethod("flow")
    solver.SetGenerator("white")
    solver.SetDenominatorPartitioning("Epstein_Nesbet")
    solver.SetDenominatorCutoff(WHITE_DENOMINATOR_CUTOFF)
    solver.SetEtaCriterion(0.0)
    solver.SetODETolerance(ode_tolerance)
    solver.SetDs(initial_step)
    solver.SetSmax(target_s)
    solver.Solve()
    imsrg_final_j = solver.GetH_s()

    production_eta = white_generator(
        production_final,
        densities,
        quanta,
        spherical_orbit_groups=groups,
    )
    production_rhs = commutator(production_eta, production_final, densities)
    imsrg_eta_j = pyimsrg.Operator(modelspace, 0, 0, 0, 2)
    imsrg_eta_j.SetAntiHermitian()
    generator = pyimsrg.Generator()
    generator.SetType("white")
    generator.SetDenominatorPartitioning("Epstein_Nesbet")
    generator.SetDenominatorCutoff(WHITE_DENOMINATOR_CUTOFF)
    generator.Update(imsrg_final_j, imsrg_eta_j)
    imsrg_rhs_j = pyimsrg.Commutator.Commutator(imsrg_eta_j, imsrg_final_j)
    imsrg_vacuum_j = imsrg_final_j.UndoNormalOrdering()
    imsrg_final = operator_to_mscheme(imsrg_final_j, reference.orbits)
    imsrg_eta = operator_to_mscheme(imsrg_eta_j, reference.orbits)
    imsrg_rhs = operator_to_mscheme(imsrg_rhs_j, reference.orbits)
    imsrg_vacuum = operator_to_mscheme(imsrg_vacuum_j, reference.orbits)
    production_vacuum = to_vacuum(production_final, densities)

    return {
        "production_flow": str(production_flow.resolve()),
        "target_s": target_s,
        "production_ode_method": metadata["ode_method"],
        "production_ode_tolerance": ode_tolerance,
        "imsrgpp_ode_method": "boost::odeint runge_kutta_dopri5",
        "imsrgpp_initial_step": initial_step,
        "h": _operator_comparison(imsrg_final, production_final),
        "eta": _operator_comparison(imsrg_eta, production_eta),
        "rhs": _operator_comparison(imsrg_rhs, production_rhs),
        "vacuum_h": _operator_comparison(imsrg_vacuum, production_vacuum),
        "denominators": _denominator_report(
            generator, modelspace, production_final, densities, reference.orbits
        ),
        "eta_norm": {
            "production_mscheme": float(
                np.sqrt(
                    np.vdot(production_eta.one_body, production_eta.one_body).real
                    + 0.25
                    * np.vdot(production_eta.two_body, production_eta.two_body).real
                )
            ),
            "imsrgpp_operator_norm": float(imsrg_eta_j.Norm()),
            "imsrgpp_one_body_norm": float(imsrg_eta_j.OneBodyNorm()),
            "imsrgpp_two_body_norm": float(imsrg_eta_j.TwoBodyNorm()),
        },
        "production_initial_zero_body": production_h_initial.zero_body,
        "production_final_zero_body": production_final.zero_body,
        "imsrgpp_final_zero_body": imsrg_final.zero_body,
    }


def run_check(args: argparse.Namespace) -> dict[str, Any]:
    repository = Path(__file__).resolve().parents[2]
    source_hashes = verify_oracle_sources(repository, args.imsrg_source)
    pyimsrg = import_pyimsrg(args.pyimsrg_dir)
    payload = read_jcoupled64(args.jcoupled64)
    reference = load_reference(args.reference, interaction_path=args.interaction)
    expected_a, expected_z = NUCLEI[args.nucleus]
    if (reference.metadata["A"], reference.metadata["Z"]) != (
        expected_a,
        expected_z,
    ):
        raise ValueError("reference A/Z does not match --nucleus")
    if reference.metadata["Nrefmax"] != 0:
        raise ValueError("strict SR check requires Nrefmax=0")

    densities = compute_densities(reference.determinants, reference.coefficients)
    validate_densities(densities, int(reference.metadata["A"]))
    lambda2_max = float(np.max(np.abs(densities.lambda2)))
    if lambda2_max > args.algebra_tolerance:
        raise ValueError(f"reference is not a Slater determinant: max|lambda2|={lambda2_max}")

    vacuum = VacuumHamiltonian(
        float(reference.metadata["zero_body"]),
        reference.one_body,
        reference.two_body,
    )
    production_h = normal_order(vacuum, densities)
    quanta = oscillator_quanta_from_orbits(reference.orbits)
    groups = spherical_orbit_groups_from_orbits(reference.orbits)
    production_eta = white_generator(
        production_h, densities, quanta, spherical_orbit_groups=groups
    )
    production_rhs_terms = commutator_contributions(
        production_eta, production_h, densities
    )
    production_rhs = commutator(production_eta, production_h, densities)

    modelspace, imsrg_vacuum = make_imsrg_operator(pyimsrg, payload, args.nucleus)
    imsrg_h_j = imsrg_vacuum.DoNormalOrdering()
    imsrg_eta_j = pyimsrg.Operator(modelspace, 0, 0, 0, 2)
    imsrg_eta_j.SetAntiHermitian()
    generator = pyimsrg.Generator()
    generator.SetType("white")
    generator.SetDenominatorPartitioning("Epstein_Nesbet")
    generator.Update(imsrg_h_j, imsrg_eta_j)
    denominator_report = _denominator_report(
        generator, modelspace, production_h, densities, reference.orbits
    )
    imsrg_rhs_j = pyimsrg.Commutator.Commutator(imsrg_eta_j, imsrg_h_j)
    contribution_report = _commutator_contribution_report(
        pyimsrg,
        modelspace,
        imsrg_eta_j,
        imsrg_h_j,
        production_rhs_terms,
        reference.orbits,
    )

    production_h_euler = MRHamiltonian(
        production_h.zero_body + args.euler_step * production_rhs.zero_body,
        production_h.one_body + args.euler_step * production_rhs.one_body,
        production_h.two_body + args.euler_step * production_rhs.two_body,
    )
    production_eta_euler = white_generator(
        production_h_euler, densities, quanta, spherical_orbit_groups=groups
    )
    production_rhs_euler = commutator(
        production_eta_euler, production_h_euler, densities
    )
    imsrg_h_euler_j = imsrg_h_j + args.euler_step * imsrg_rhs_j
    imsrg_eta_euler_j = pyimsrg.Operator(modelspace, 0, 0, 0, 2)
    imsrg_eta_euler_j.SetAntiHermitian()
    generator.Update(imsrg_h_euler_j, imsrg_eta_euler_j)
    imsrg_rhs_euler_j = pyimsrg.Commutator.Commutator(
        imsrg_eta_euler_j, imsrg_h_euler_j
    )

    def production_rhs_at(hamiltonian: MRHamiltonian) -> MRHamiltonian:
        eta = white_generator(
            hamiltonian, densities, quanta, spherical_orbit_groups=groups
        )
        return commutator(eta, hamiltonian, densities)

    def imsrg_rhs_at(hamiltonian: Any) -> Any:
        eta = pyimsrg.Operator(modelspace, 0, 0, 0, 2)
        eta.SetAntiHermitian()
        generator.Update(hamiltonian, eta)
        return pyimsrg.Commutator.Commutator(eta, hamiltonian)

    production_rk = production_h
    imsrg_rk_j = imsrg_h_j
    rk_pairs = []
    for step_index in range(1, args.rk4_steps + 1):
        production_rk = _mr_rk4_step(
            production_rk, args.rk4_step, production_rhs_at
        )
        imsrg_rk_j = _imsrg_rk4_step(imsrg_rk_j, args.rk4_step, imsrg_rhs_at)
        rk_pairs.append((step_index * args.rk4_step, production_rk, imsrg_rk_j))

    imsrg_vacuum_m = operator_to_mscheme(imsrg_vacuum, reference.orbits)
    imsrg_h_m = operator_to_mscheme(imsrg_h_j, reference.orbits)
    imsrg_eta_m = operator_to_mscheme(imsrg_eta_j, reference.orbits)
    imsrg_rhs_m = operator_to_mscheme(imsrg_rhs_j, reference.orbits)
    imsrg_h_euler_m = operator_to_mscheme(imsrg_h_euler_j, reference.orbits)
    imsrg_eta_euler_m = operator_to_mscheme(imsrg_eta_euler_j, reference.orbits)
    imsrg_rhs_euler_m = operator_to_mscheme(imsrg_rhs_euler_j, reference.orbits)
    rk4_checkpoints = [
        {
            "s": s,
            "h": _operator_comparison(
                operator_to_mscheme(imsrg_h_checkpoint, reference.orbits),
                production_h_checkpoint,
            ),
        }
        for s, production_h_checkpoint, imsrg_h_checkpoint in rk_pairs
    ]
    production_eta_rk = white_generator(
        production_rk, densities, quanta, spherical_orbit_groups=groups
    )
    production_rhs_rk = commutator(production_eta_rk, production_rk, densities)
    imsrg_eta_rk_j = pyimsrg.Operator(modelspace, 0, 0, 0, 2)
    imsrg_eta_rk_j.SetAntiHermitian()
    generator.Update(imsrg_rk_j, imsrg_eta_rk_j)
    imsrg_rhs_rk_j = pyimsrg.Commutator.Commutator(imsrg_eta_rk_j, imsrg_rk_j)

    result = {
        "schema": "mrimsrg_sr_imsrgpp_check_v1",
        "nucleus": args.nucleus,
        "reference": str(args.reference.resolve()),
        "jcoupled64": str(args.jcoupled64.resolve()),
        "interaction": str(args.interaction.resolve()),
        "pyimsrg_module": str(Path(pyimsrg.__file__).resolve()),
        "imsrg_source": str(args.imsrg_source.resolve()),
        "oracle_source_sha256": source_hashes,
        "lambda2_max_abs": lambda2_max,
        "mask": _mask_report(densities, quanta),
        "denominators": denominator_report,
        "vacuum_input": _operator_comparison(imsrg_vacuum_m, vacuum),
        "normal_ordered_h": _operator_comparison(imsrg_h_m, production_h),
        "eta_s0": _operator_comparison(imsrg_eta_m, production_eta),
        "rhs_s0": _operator_comparison(imsrg_rhs_m, production_rhs),
        "commutator_terms_s0": contribution_report,
        "eta_norm_s0": {
            "production_mscheme": float(
                np.sqrt(
                    np.vdot(production_eta.one_body, production_eta.one_body).real
                    + 0.25
                    * np.vdot(production_eta.two_body, production_eta.two_body).real
                )
            ),
            "imsrgpp_operator_norm": float(imsrg_eta_j.Norm()),
            "imsrgpp_one_body_norm": float(imsrg_eta_j.OneBodyNorm()),
            "imsrgpp_two_body_norm": float(imsrg_eta_j.TwoBodyNorm()),
        },
        "euler_step": args.euler_step,
        "h_after_euler": _operator_comparison(
            imsrg_h_euler_m, production_h_euler
        ),
        "eta_after_euler": _operator_comparison(
            imsrg_eta_euler_m, production_eta_euler
        ),
        "rhs_after_euler": _operator_comparison(
            imsrg_rhs_euler_m, production_rhs_euler
        ),
        "rk4_step": args.rk4_step,
        "rk4_checkpoints": rk4_checkpoints,
        "eta_after_rk4": _operator_comparison(
            operator_to_mscheme(imsrg_eta_rk_j, reference.orbits),
            production_eta_rk,
        ),
        "rhs_after_rk4": _operator_comparison(
            operator_to_mscheme(imsrg_rhs_rk_j, reference.orbits),
            production_rhs_rk,
        ),
    }
    if args.production_flow is not None:
        result["full_flow"] = _full_flow_report(
            pyimsrg,
            modelspace,
            imsrg_h_j,
            args.production_flow,
            production_h,
            densities,
            quanta,
            groups,
            reference,
            args.nucleus,
            args.full_flow_ode_tolerance,
            args.full_flow_initial_step,
        )
    return result


def _passes(
    report: dict[str, Any], algebra_tolerance: float, full_flow_tolerance: float
) -> bool:
    mask = report["mask"]
    if mask["sr_one_body_removed_by_delta_e"] or mask["sr_two_body_removed_by_delta_e"]:
        return False
    for rank in ("one_body", "two_body"):
        if report["denominators"][rank]["max_abs"] > algebra_tolerance:
            return False
    for name in (
        "vacuum_input",
        "normal_ordered_h",
        "eta_s0",
        "rhs_s0",
        "h_after_euler",
        "eta_after_euler",
        "rhs_after_euler",
        "eta_after_rk4",
        "rhs_after_rk4",
    ):
        comparison = report[name]
        for rank in ("zero_body", "one_body", "two_body"):
            if comparison[rank]["max_abs"] > algebra_tolerance:
                return False
    for checkpoint in report["rk4_checkpoints"]:
        for rank in ("zero_body", "one_body", "two_body"):
            if checkpoint["h"][rank]["max_abs"] > algebra_tolerance:
                return False
    for comparison in report["commutator_terms_s0"].values():
        for rank in ("zero_body", "one_body", "two_body"):
            if comparison[rank]["max_abs"] > algebra_tolerance:
                return False
    if "full_flow" in report:
        for name in ("h", "eta", "rhs", "vacuum_h"):
            comparison = report["full_flow"][name]
            for rank in ("zero_body", "one_body", "two_body"):
                if comparison[rank]["max_abs"] > full_flow_tolerance:
                    return False
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nucleus", choices=("He4", "O16"), required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--jcoupled64", type=Path, required=True)
    parser.add_argument(
        "--interaction",
        type=Path,
        default=Path(
            "/home/mengziyan/Forces/N2LO_opt/"
            "TwBME_N2LO_opt_hw20_emax2_e2max4.minipack"
        ),
    )
    parser.add_argument(
        "--pyimsrg-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "build-src" / "src",
    )
    parser.add_argument(
        "--imsrg-source", type=Path, default=Path(__file__).resolve().parents[2]
    )
    parser.add_argument("--algebra-tolerance", type=float, default=1e-10)
    parser.add_argument("--euler-step", type=float, default=1e-4)
    parser.add_argument("--rk4-step", type=float, default=1e-3)
    parser.add_argument("--rk4-steps", type=int, default=3)
    parser.add_argument("--production-flow", type=Path)
    parser.add_argument("--full-flow-ode-tolerance", type=float, default=1e-8)
    parser.add_argument("--full-flow-initial-step", type=float, default=1e-2)
    parser.add_argument("--full-flow-tolerance", type=float, default=1e-5)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    report = run_check(args)
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.json is not None:
        if args.json.exists():
            raise FileExistsError(f"refusing to overwrite existing report: {args.json}")
        args.json.write_text(rendered + "\n", encoding="utf-8")
    if not _passes(report, args.algebra_tolerance, args.full_flow_tolerance):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
