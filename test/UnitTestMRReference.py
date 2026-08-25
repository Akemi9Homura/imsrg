#!/usr/bin/env python3

import math
from pathlib import Path
import sys
import tempfile

import numpy as np
import pyIMSRG

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from prototype.mrimsrg.jcoupling import (  # noqa: E402
    CoupledBlock,
    couple_scalar_two_body,
    extract_j_orbits,
    mr_lambda2_one_body_coupled,
)
from prototype.mrimsrg.basis import prepare_natural_basis  # noqa: E402
from prototype.mrimsrg.densities import compute_densities  # noqa: E402
from prototype.mrimsrg.export_jref import export_reference  # noqa: E402
from prototype.mrimsrg.reference_io import load_reference  # noqa: E402


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

# The block-matrix implementation must match the slow spherical-orbit C++
# reference topology by topology.
optimized_parts = pyIMSRG.MR_comm221_lambda2(
    x_operator, y_operator, algebraic_reference
)
for name in ("IV", "V", "VI"):
    optimized = getattr(optimized_parts, name)
    reference_part = getattr(cpp_parts, name)
    maximum = max(
        abs(optimized(i, j) - reference_part(i, j))
        for i in range(ms.GetNumberOrbits())
        for j in range(ms.GetNumberOrbits())
    )
    assert maximum < 2e-11, f"optimized MR topology {name} differs by {maximum}"

# Exact lambda2=0 uses the existing C++ SR dispatcher and is bitwise unchanged.
slater_reference = pyIMSRG.MRReference(ms, 4, 2, 0, occupations)
sr_commutator = pyIMSRG.Commutator.Commutator(x_operator, y_operator)
mr_slater_commutator = pyIMSRG.MR_Commutator(
    x_operator, y_operator, slater_reference
)
assert max_operator_difference(sr_commutator, mr_slater_commutator) == 0.0

# For nonzero lambda2, the full entry adds the verified 1B term and contracts
# lambda2 with the already completed SR two-body commutator for the 0B term.
mr_commutator = pyIMSRG.MR_Commutator(
    x_operator, y_operator, algebraic_reference
)
expected_one_body = sr_commutator.OneBody + optimized_parts.Total()
assert (mr_commutator.OneBody - expected_one_body).Norm() < 2e-11
assert (mr_commutator.TwoBody - sr_commutator.TwoBody).Norm() == 0.0
assert math.isclose(
    mr_commutator.ZeroBody - sr_commutator.ZeroBody,
    algebraic_reference.ContractLambda2(sr_commutator.TwoBody),
    rel_tol=0.0,
    abs_tol=2e-11,
)

# Exercise the complete saved-wavefunction -> natural density -> compact
# J-scheme bridge -> C++ reader path with the physical correlated Be8 Nrefmax=0
# reference.  This also freezes the provenance digests from the acceptance log.
repository_root = Path(__file__).resolve().parents[1]
be8_path = repository_root / "prototype/mrimsrg/data/Be8_Nrefmax0_final"
be8_data = load_reference(be8_path)
be8_natural = prepare_natural_basis(
    compute_densities(be8_data.determinants, be8_data.coefficients)
)
be8_orbits = extract_j_orbits(be8_data.orbits)
be8_blocks = couple_scalar_two_body(
    be8_natural.densities.lambda2, be8_data.orbits
)
be8_ms = pyIMSRG.ModelSpace(2, "He4", "He4")
with tempfile.TemporaryDirectory() as temporary_directory:
    bridge_path = Path(temporary_directory) / "Be8_Nrefmax0.jref"
    summary = export_reference(be8_path, bridge_path)
    be8_ms.SetReferenceOcc(
        pyIMSRG.MRReference.ReadOccupationMap(be8_ms, str(bridge_path))
    )
    be8_reference = pyIMSRG.MRReference.ReadBinary(be8_ms, str(bridge_path))
assert summary["rdm_sha256"] == "9e8da10770cf143dcc87382c96aabf64961263fac86a02ac13ec6a468121c538"
assert summary["wavefunction_sha256"] == "20f36fd7da0461e2e7b7f92b98662d0241be69f4482df6a1655aaeae02be469b"
assert be8_reference.rdm_sha256 == summary["rdm_sha256"]
assert be8_reference.wavefunction_sha256 == summary["wavefunction_sha256"]
assert be8_reference.MaximumContractionViolation() < 1e-12
assert be8_reference.MaximumHermiticityViolation() < 1e-12
for i in range(be8_ms.GetNumberOrbits()):
    for j in range(be8_ms.GetNumberOrbits()):
        expected = 1.0 if i == j else 0.0
        assert abs(be8_reference.NaturalOrbitTransformation(i, j) - expected) < 1e-14
for block in be8_blocks:
    channel = be8_ms.GetTwoBodyChannelIndex(block.J, block.parity, block.Tz)
    two_body_channel = be8_ms.GetTwoBodyChannel(channel)
    for ibra, (a_file, b_file) in enumerate(block.pairs):
        oa, ob = be8_orbits[a_file], be8_orbits[b_file]
        a = be8_ms.GetOrbitIndex(oa.n, oa.l, oa.j2, oa.tz2)
        b = be8_ms.GetOrbitIndex(ob.n, ob.l, ob.j2, ob.tz2)
        local_bra = two_body_channel.GetLocalIndex(a, b)
        for iket, (c_file, d_file) in enumerate(block.pairs):
            oc, od = be8_orbits[c_file], be8_orbits[d_file]
            c = be8_ms.GetOrbitIndex(oc.n, oc.l, oc.j2, oc.tz2)
            d = be8_ms.GetOrbitIndex(od.n, od.l, od.j2, od.tz2)
            local_ket = two_body_channel.GetLocalIndex(c, d)
            actual = be8_reference.Lambda2.GetTBME_norm_chij(
                channel, channel, local_bra, local_ket
            )
            assert abs(actual - block.matrix[ibra, iket]) < 1e-13

# A correlated closed-shell Nrefmax=2 reference exercises the non-identity
# temporary natural-orbit rotation, including radial mixing in fixed (l,j,tz)
# blocks.  It must survive the same bridge without being mistaken for HO basis.
he4_nref2_path = repository_root / "prototype/mrimsrg/data/He4_Nrefmax2"
he4_nref2_ms = pyIMSRG.ModelSpace(2, "He4", "He4")
with tempfile.TemporaryDirectory() as temporary_directory:
    bridge_path = Path(temporary_directory) / "He4_Nrefmax2.jref"
    summary = export_reference(he4_nref2_path, bridge_path)
    he4_nref2_ms.SetReferenceOcc(
        pyIMSRG.MRReference.ReadOccupationMap(he4_nref2_ms, str(bridge_path))
    )
    he4_nref2_reference = pyIMSRG.MRReference.ReadBinary(
        he4_nref2_ms, str(bridge_path)
    )
assert not summary["natural_basis_is_identity"]
assert he4_nref2_reference.rdm_sha256 == "b160ba51f981af596d960a03ace8a53ba777093c77afe325ac740e6337b90df2"
assert he4_nref2_reference.MaximumContractionViolation() < 1e-11
assert he4_nref2_reference.Lambda2.Norm() > 1e-3
maximum_off_diagonal_transformation = max(
    abs(he4_nref2_reference.NaturalOrbitTransformation(i, j))
    for i in range(he4_nref2_ms.GetNumberOrbits())
    for j in range(he4_nref2_ms.GetNumberOrbits())
    if i != j
)
assert maximum_off_diagonal_transformation > 1e-3

print("MRReference tests passed")
