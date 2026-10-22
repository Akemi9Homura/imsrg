#!/usr/bin/env python3

"""Configuration guards of the MR-IMSRG(2) path.

1. MR_Commutator refuses the global IMSRG(3) switch (it would return 3B parts
   without their MR terms) and still works once the switch is off.
2. An MR BCH transform that does not converge within the nested-commutator
   limit raises, while the SR path keeps its historical warning and a zero
   threshold (fixed-order expansion used by the oracle tests) is accepted.
3. MRReference.Validate refuses occupations in [1e-12, ModelSpace::OCC_CUT),
   which the SR contractions would round to zero; numerical zeros below
   1e-12 (He4 Nrefmax=2 has 0d occupations ~1e-25) are accepted.
"""

from pathlib import Path
import sys
import tempfile

import pyIMSRG

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))

from prototype.mrimsrg.export_jref import export_reference  # noqa: E402


def expect_failure(function, exception, needle):
    try:
        function()
    except exception as error:
        assert needle in str(error), str(error)
        return
    raise AssertionError(f"expected {exception.__name__} containing {needle!r}")


with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp) / "He4_Nrefmax2.jref"
    export_reference(REPOSITORY / "prototype/mrimsrg/data/He4_Nrefmax2", path)
    ms = pyIMSRG.ModelSpace(2, "He4", "He4")
    occupation_map = pyIMSRG.MRReference.ReadOccupationMap(ms, str(path))
    ms.SetReferenceOcc(occupation_map)
    reference = pyIMSRG.MRReference.ReadBinary(ms, str(path))

test = pyIMSRG.UnitTest(ms)
test.SetRandomSeed(20261022)
X = test.RandomOp(ms, 0, 0, 0, 2, -1)
Y = test.RandomOp(ms, 0, 0, 0, 2, +1)

# 1. IMSRG(3) switch.
baseline = pyIMSRG.MR_Commutator(X, Y, reference)
pyIMSRG.SetUseIMSRG3(True)
try:
    expect_failure(lambda: pyIMSRG.MR_Commutator(X, Y, reference),
                   ValueError, "MR-IMSRG(2) only")
finally:
    pyIMSRG.SetUseIMSRG3(False)
again = pyIMSRG.MR_Commutator(X, Y, reference)
# SR reductions over channels use dynamic OpenMP scheduling, so compare
# with a tolerance rather than bitwise.
assert abs(again.ZeroBody - baseline.ZeroBody) <= 1e-10 * max(1.0, abs(baseline.ZeroBody))
assert (again - baseline).Norm() <= 1e-10 * max(1.0, baseline.Norm())

# 2. Unconverged BCH: a large generator keeps the nested commutators growing.
omega = X * 50.0
mr_solver = pyIMSRG.IMSRGSolver(Y)
mr_solver.SetMRReference(reference)
expect_failure(lambda: mr_solver.EvaluateBCHTransform(Y, omega),
               RuntimeError, "did not converge")
sr_solver = pyIMSRG.IMSRGSolver(Y)
sr_solver.EvaluateBCHTransform(Y, omega)  # SR: warning only, no exception
pyIMSRG.BCH.Set_BCH_Transform_Threshold(0.0)
try:
    mr_solver.EvaluateBCHTransform(Y, omega)  # fixed 40 orders requested
finally:
    pyIMSRG.BCH.Set_BCH_Transform_Threshold(1e-9)

# 3. Occupations in [1e-12, OCC_CUT).
tiny = dict(occupation_map)
orbit = min((o for o, n in tiny.items() if 0.0 < n < 1.0), key=lambda o: tiny[o])
tiny[orbit] = 5e-7
ms_tiny = pyIMSRG.ModelSpace(2, "He4", "He4")
ms_tiny.SetReferenceOcc(tiny)
occupations = [tiny.get(i, 0.0) for i in range(ms_tiny.GetNumberOrbits())]
candidate = pyIMSRG.MRReference(ms_tiny, 4, 2, 2, occupations)
expect_failure(lambda: candidate.Validate(1e-10), RuntimeError, "OCC_CUT")

print("MR guards: IMSRG(3) switch, unconverged MR BCH and [1e-12, OCC_CUT) occupations refused")
