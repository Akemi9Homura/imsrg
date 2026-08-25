#!/usr/bin/env python3

import math

import pyIMSRG


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

print("MRReference tests passed")
