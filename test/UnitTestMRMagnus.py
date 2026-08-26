#!/usr/bin/env python3

"""Matrix-element oracle for the production scalar MR-Magnus BCH path."""

import sys
from pathlib import Path

import numpy as np
import pyIMSRG


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))

from prototype.mrimsrg.densities import compute_densities  # noqa: E402
from prototype.mrimsrg.jcoupling import (  # noqa: E402
    couple_scalar_two_body,
    extract_j_orbits,
)
from prototype.mrimsrg.magnus import (  # noqa: E402
    bch_product,
    bch_transform,
    magnus_derivative,
)
from prototype.mrimsrg.sr_imsrgpp_check import operator_to_mscheme  # noqa: E402


def random_scalar_operator(modelspace, layout, rng, adjoint_sign):
    operator = pyIMSRG.Operator(modelspace)
    if adjoint_sign > 0:
        operator.SetHermitian()
    else:
        operator.SetAntiHermitian()

    groups = {}
    for index in range(modelspace.GetNumberOrbits()):
        orbit = modelspace.GetOrbit(index)
        groups.setdefault((orbit.l, orbit.j2, orbit.tz2), []).append(index)
    for indices in groups.values():
        raw = rng.normal(size=(len(indices), len(indices)))
        matrix = 0.5 * (raw + adjoint_sign * raw.T)
        for row, p in enumerate(indices):
            for column, q in enumerate(indices):
                operator.SetOneBody(p, q, float(matrix[row, column]))

    for block in layout:
        channel = modelspace.GetTwoBodyChannelIndex(block.J, block.parity, block.Tz)
        two_body_channel = modelspace.GetTwoBodyChannel(channel)
        raw = rng.normal(size=block.matrix.shape)
        matrix = 0.5 * (raw + adjoint_sign * raw.T)
        for ibra, pair in enumerate(block.pairs):
            local_bra = two_body_channel.GetLocalIndex(*pair)
            first_ket = ibra if adjoint_sign > 0 else ibra + 1
            for iket in range(first_ket, len(block.pairs)):
                local_ket = two_body_channel.GetLocalIndex(*block.pairs[iket])
                operator.TwoBody.SetTBME_chij(
                    channel, channel, local_bra, local_ket, float(matrix[ibra, iket])
                )
    return operator


def errors(actual, expected):
    return (
        abs(actual.zero_body - expected.zero_body),
        float(np.max(np.abs(actual.one_body - expected.one_body))),
        float(np.max(np.abs(actual.two_body - expected.two_body))),
    )


def j_operator_error(left, right):
    return max(
        abs(left.ZeroBody - right.ZeroBody),
        (left.OneBody - right.OneBody).Norm(),
        (left.TwoBody - right.TwoBody).Norm(),
    )


def make_correlated_reference(modelspace):
    rows = []
    for index in range(modelspace.GetNumberOrbits()):
        orbit = modelspace.GetOrbit(index)
        for m2 in range(-orbit.j2, orbit.j2 + 1, 2):
            rows.append([index, orbit.n, orbit.l, orbit.j2, m2, orbit.tz2])
    orbits = np.asarray(rows, dtype=np.int32)

    # Both determinants fill a complete proton and neutron j=1/2 subshell,
    # hence each is J=0+, while their coherent mixture has nonzero lambda2.
    determinants = np.zeros((2, len(orbits)), dtype=np.uint8)
    for m_index, (_, n, l, j2, _, _) in enumerate(orbits):
        if (n, l, j2) == (0, 0, 1):
            determinants[0, m_index] = 1
        if (n, l, j2) == (0, 1, 1):
            determinants[1, m_index] = 1
    assert np.all(np.sum(determinants, axis=1) == 4)
    densities = compute_densities(
        determinants, np.array([np.sqrt(0.4), -np.sqrt(0.6)])
    )

    j_orbits = extract_j_orbits(orbits)
    occupations = [
        float(densities.gamma1[orbit.substates[0], orbit.substates[0]])
        for orbit in j_orbits
    ]
    modelspace.SetReferenceOcc(
        {index: occupation for index, occupation in enumerate(occupations)}
    )
    reference = pyIMSRG.MRReference(modelspace, 4, 2, 2, occupations)
    lambda_blocks = couple_scalar_two_body(densities.lambda2, orbits)
    for block in lambda_blocks:
        channel = modelspace.GetTwoBodyChannelIndex(block.J, block.parity, block.Tz)
        two_body_channel = modelspace.GetTwoBodyChannel(channel)
        for ibra, pair in enumerate(block.pairs):
            local_bra = two_body_channel.GetLocalIndex(*pair)
            for iket in range(ibra, len(block.pairs)):
                local_ket = two_body_channel.GetLocalIndex(*block.pairs[iket])
                reference.Lambda2.SetTBME_chij(
                    channel,
                    channel,
                    local_bra,
                    local_ket,
                    float(block.matrix[ibra, iket]),
                )
    reference.Validate(2e-10)
    assert reference.Lambda2.Norm() > 1e-6
    return reference, densities, orbits, lambda_blocks


def main():
    modelspace = pyIMSRG.ModelSpace(1, "He4", "He4")
    modelspace.SetHbarOmega(20.0)
    reference, densities, orbits, layout = make_correlated_reference(modelspace)

    rng = np.random.default_rng(20260826)
    hamiltonian = random_scalar_operator(modelspace, layout, rng, +1)
    omega = 0.015 * random_scalar_operator(modelspace, layout, rng, -1)
    d_omega = 0.001 * random_scalar_operator(modelspace, layout, rng, -1)
    oracle_h = operator_to_mscheme(hamiltonian, orbits)
    oracle_omega = operator_to_mscheme(omega, orbits)
    oracle_domega = operator_to_mscheme(d_omega, orbits)

    solver = pyIMSRG.IMSRGSolver(hamiltonian)
    solver.SetMRReference(reference)
    pyIMSRG.BCH.Set_BCH_Transform_Threshold(0.0)
    pyIMSRG.BCH.Set_BCH_Product_Threshold(0.0)
    try:
        actual_bch = operator_to_mscheme(
            solver.EvaluateBCHTransform(hamiltonian, omega), orbits
        )
        expected_bch = bch_transform(
            oracle_h, oracle_omega, densities, max_order=40
        ).operator
        bch_errors = errors(actual_bch, expected_bch)
        print(
            "MR-BCH random J->m: "
            f"zero={bch_errors[0]:.3e} one={bch_errors[1]:.3e} "
            f"two={bch_errors[2]:.3e}"
        )
        assert max(bch_errors) < 1e-10

        actual_product = operator_to_mscheme(
            solver.EvaluateBCHProduct(d_omega, omega), orbits
        )
        expected_product = bch_product(
            oracle_domega, oracle_omega, densities, threshold=0.0
        ).omega
        product_errors = errors(actual_product, expected_product)
        print(
            "MR-BCH product random J->m: "
            f"zero={product_errors[0]:.3e} one={product_errors[1]:.3e} "
            f"two={product_errors[2]:.3e}"
        )
        assert max(product_errors) < 1e-10

        actual_derivative = operator_to_mscheme(
            solver.EvaluateMagnusDerivative(omega, d_omega), orbits
        )
        expected_derivative = magnus_derivative(
            oracle_omega, oracle_domega, densities
        )
        derivative_errors = errors(actual_derivative, expected_derivative)
        print(
            "MR-Magnus derivative random J->m: "
            f"zero={derivative_errors[0]:.3e} "
            f"one={derivative_errors[1]:.3e} "
            f"two={derivative_errors[2]:.3e}"
        )
        assert max(derivative_errors) < 1e-10
    finally:
        pyIMSRG.BCH.Set_BCH_Transform_Threshold(1e-9)
        pyIMSRG.BCH.Set_BCH_Product_Threshold(1e-4)

    # The explicit MR dispatcher must reduce to the unmodified production SR
    # BCH path before any floating-point lambda2 additions are attempted.
    sr_modelspace = pyIMSRG.ModelSpace(1, "He4", "He4")
    sr_rows = []
    for index in range(sr_modelspace.GetNumberOrbits()):
        orbit = sr_modelspace.GetOrbit(index)
        for m2 in range(-orbit.j2, orbit.j2 + 1, 2):
            sr_rows.append([index, orbit.n, orbit.l, orbit.j2, m2, orbit.tz2])
    sr_orbits = np.asarray(sr_rows, dtype=np.int32)
    sr_layout = couple_scalar_two_body(
        np.zeros((len(sr_orbits),) * 4), sr_orbits
    )
    sr_rng = np.random.default_rng(8675309)
    sr_h = random_scalar_operator(sr_modelspace, sr_layout, sr_rng, +1)
    sr_omega = 0.015 * random_scalar_operator(
        sr_modelspace, sr_layout, sr_rng, -1
    )
    sr_domega = 0.001 * random_scalar_operator(
        sr_modelspace, sr_layout, sr_rng, -1
    )
    occupations = [
        sr_modelspace.GetOrbit(index).occ
        for index in range(sr_modelspace.GetNumberOrbits())
    ]
    slater_reference = pyIMSRG.MRReference(
        sr_modelspace, 4, 2, 0, occupations
    )
    slater_reference.Validate(1e-12)
    sr_solver = pyIMSRG.IMSRGSolver(sr_h)
    mr_solver = pyIMSRG.IMSRGSolver(sr_h)
    mr_solver.SetMRReference(slater_reference)
    pyIMSRG.BCH.Set_BCH_Transform_Threshold(0.0)
    pyIMSRG.BCH.Set_BCH_Product_Threshold(0.0)
    try:
        sr_bch = sr_solver.EvaluateBCHTransform(sr_h, sr_omega)
        mr_bch = mr_solver.EvaluateBCHTransform(sr_h, sr_omega)
        sr_product = sr_solver.EvaluateBCHProduct(sr_domega, sr_omega)
        mr_product = mr_solver.EvaluateBCHProduct(sr_domega, sr_omega)
    finally:
        pyIMSRG.BCH.Set_BCH_Transform_Threshold(1e-9)
        pyIMSRG.BCH.Set_BCH_Product_Threshold(1e-4)
    bch_sr_error = j_operator_error(sr_bch, mr_bch)
    product_sr_error = j_operator_error(sr_product, mr_product)
    print(
        f"MR->SR BCH degeneration: transform={bch_sr_error:.3e} "
        f"product={product_sr_error:.3e}"
    )
    # The two solvers evaluate the same OpenMP SR reductions independently;
    # their roundoff need not be bitwise identical even though lambda2 takes
    # the exact early-return branch in MRCommutator.
    assert bch_sr_error < 1e-14
    assert product_sr_error < 1e-14

    # Direct RK4 and the production first-order Magnus update represent the
    # same continuous MR flow.  Their one-step difference must therefore
    # vanish quadratically as the common interval is reduced.
    flow_hamiltonian = 0.1 * random_scalar_operator(
        modelspace, layout, np.random.default_rng(424242), +1
    )
    for p in range(modelspace.GetNumberOrbits()):
        orbit_p = modelspace.GetOrbit(p)
        for q in range(modelspace.GetNumberOrbits()):
            flow_hamiltonian.SetOneBody(p, q, 0.0)
        flow_hamiltonian.SetOneBody(
            p,
            p,
            5.0 * (2 * orbit_p.n + orbit_p.l)
            + 0.05 * p,
        )

    convergence_errors = []
    for step in (1e-2, 5e-3, 2.5e-3):
        direct_solver = pyIMSRG.IMSRGSolver(flow_hamiltonian)
        direct_solver.SetMRReference(reference)
        direct_solver.SetGenerator("white-ncsm")
        direct_solver.SetMethod("flow_RK4")
        direct_solver.SetEtaCriterion(0.0)
        direct_solver.SetDs(step)
        direct_solver.SetDsmax(step)
        direct_solver.SetSmax(step)
        direct_solver.Solve()

        magnus_solver = pyIMSRG.IMSRGSolver(flow_hamiltonian)
        magnus_solver.SetMRReference(reference)
        magnus_solver.SetGenerator("white-ncsm")
        magnus_solver.SetMethod("magnus")
        magnus_solver.SetMagnusAdaptive(False)
        magnus_solver.SetEtaCriterion(0.0)
        magnus_solver.SetDs(step)
        magnus_solver.SetSmax(step)
        pyIMSRG.BCH.Set_BCH_Transform_Threshold(0.0)
        pyIMSRG.BCH.Set_BCH_Product_Threshold(0.0)
        try:
            magnus_solver.Solve()
        finally:
            pyIMSRG.BCH.Set_BCH_Transform_Threshold(1e-9)
            pyIMSRG.BCH.Set_BCH_Product_Threshold(1e-4)
        convergence_errors.append(
            j_operator_error(direct_solver.GetH_s(), magnus_solver.GetH_s())
        )

    ratios = [
        convergence_errors[index + 1] / convergence_errors[index]
        for index in range(2)
    ]
    print(
        "MR direct/Magnus common limit: errors="
        + ",".join(f"{error:.3e}" for error in convergence_errors)
        + " ratios="
        + ",".join(f"{ratio:.3f}" for ratio in ratios)
    )
    assert all(0.15 < ratio < 0.35 for ratio in ratios)


if __name__ == "__main__":
    main()
