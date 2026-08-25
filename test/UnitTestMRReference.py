#!/usr/bin/env python3

import math
from pathlib import Path
import sys

import numpy as np
import pyIMSRG

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from prototype.mrimsrg.jcoupling import (  # noqa: E402
    CoupledBlock,
    couple_scalar_two_body,
    mr_lambda2_one_body_coupled,
)


def max_operator_difference(left, right):
    return max(
        abs(left.ZeroBody - right.ZeroBody),
        (left.OneBody - right.OneBody).Norm(),
        (left.TwoBody - right.TwoBody).Norm(),
    )


ms = pyIMSRG.ModelSpace(1, "He4", "He4")
occupations = [ms.GetOrbit(i).occ for i in range(ms.GetNumberOrbits())]
reference = pyIMSRG.MRReference(ms, 4, 2, 0, occupations)

# A Slater reference has lambda2=0 and must satisfy the cumulant contraction.
reference.Validate(1e-12)
assert reference.MaximumHermiticityViolation() == 0.0
assert reference.MaximumContractionViolation() == 0.0

# The proton trace follows the imsrg++ tz2=-1 convention.  An asymmetric
# nucleus prevents a proton/neutron swap from passing unnoticed.
ms_he8 = pyIMSRG.ModelSpace(1, "He8", "He8")
occupations_he8 = [ms_he8.GetOrbit(i).occ for i in range(ms_he8.GetNumberOrbits())]
pyIMSRG.MRReference(ms_he8, 8, 2, 0, occupations_he8).Validate(1e-12)

unit_test = pyIMSRG.UnitTest(ms)
unit_test.SetRandomSeed(8675309)
vacuum = unit_test.RandomOp(ms, 0, 0, 0, 2, +1)

# Exercise a nonzero normalized-pair lambda2. It is deliberately not required
# to be a physical density for this algebraic normal-ordering round trip.
lambda_source = unit_test.RandomOp(ms, 0, 0, 0, 2, +1)
reference.Lambda2 = lambda_source.TwoBody
mr_operator = reference.NormalOrder(vacuum)
recovered = reference.UndoNormalOrder(mr_operator)
error = max_operator_difference(recovered, vacuum)
assert error < 1e-10, f"MR normal-ordering round trip failed: {error}"

# The extra zero-body contribution is exactly the independently exposed
# C2-lambda2 contraction, while f and Gamma reuse the existing implementation.
sr_operator = vacuum.DoNormalOrdering()
expected_extra = reference.ContractLambda2(vacuum.TwoBody)
assert math.isclose(
    mr_operator.ZeroBody - sr_operator.ZeroBody,
    expected_extra,
    rel_tol=0.0,
    abs_tol=1e-11,
)
assert (mr_operator.OneBody - sr_operator.OneBody).Norm() < 1e-12
assert (mr_operator.TwoBody - sr_operator.TwoBody).Norm() < 1e-12

# Compare the slow C++ spherical-orbit IV/V/VI implementation to the independent
# Python J-scheme oracle.  The latter has its own direct m-scheme topology test.
rows = []
for index in range(ms.GetNumberOrbits()):
    orbit = ms.GetOrbit(index)
    for m2 in range(-orbit.j2, orbit.j2 + 1, 2):
        rows.append([index, orbit.n, orbit.l, orbit.j2, m2, orbit.tz2])
orbits = np.asarray(rows, dtype=np.int32)
layout = couple_scalar_two_body(np.zeros((len(rows),) * 4), orbits)
rng = np.random.default_rng(314159)


def random_blocks(hermiticity):
    blocks = []
    for block in layout:
        raw = rng.normal(size=block.matrix.shape)
        matrix = 0.5 * (raw + hermiticity * raw.T)
        blocks.append(CoupledBlock(block.J, block.parity, block.Tz, block.pairs, matrix))
    return tuple(blocks)


def fill_two_body(two_body, blocks, hermiticity):
    for block in blocks:
        channel = ms.GetTwoBodyChannelIndex(block.J, block.parity, block.Tz)
        two_body_channel = ms.GetTwoBodyChannel(channel)
        for ibra, pair in enumerate(block.pairs):
            assert two_body_channel.GetLocalIndex(*pair) == ibra
            first_ket = ibra + 1 if hermiticity < 0 else ibra
            for iket in range(first_ket, len(block.pairs)):
                two_body.SetTBME_chij(
                    channel, channel, ibra, iket, float(block.matrix[ibra, iket])
                )


x_blocks = random_blocks(-1)
y_blocks = random_blocks(+1)
lambda_blocks = random_blocks(+1)
x_operator = pyIMSRG.Operator(ms)
x_operator.SetAntiHermitian()
y_operator = pyIMSRG.Operator(ms)
y_operator.SetHermitian()
algebraic_reference = pyIMSRG.MRReference(ms, 4, 2, 0, occupations)
fill_two_body(x_operator.TwoBody, x_blocks, -1)
fill_two_body(y_operator.TwoBody, y_blocks, +1)
fill_two_body(algebraic_reference.Lambda2, lambda_blocks, +1)

cpp_parts = pyIMSRG.MR_comm221_lambda2_reference(
    x_operator, y_operator, algebraic_reference
)
python_parts = mr_lambda2_one_body_coupled(
    x_blocks, y_blocks, lambda_blocks, orbits, -1, +1
)
for cpp_name, python_name in (("IV", "iv"), ("V", "v"), ("VI", "vi")):
    cpp_matrix = getattr(cpp_parts, cpp_name)
    cpp_values = np.asarray(
        [
            [cpp_matrix(i, j) for j in range(ms.GetNumberOrbits())]
            for i in range(ms.GetNumberOrbits())
        ]
    )
    np.testing.assert_allclose(
        cpp_values, getattr(python_parts, python_name), atol=2e-11, rtol=0
    )

print("MRReference tests passed")
