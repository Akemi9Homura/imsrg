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
    reconstruct_scalar_two_body,
)
from prototype.mrimsrg.basis import (  # noqa: E402
    prepare_natural_basis,
    transform_hamiltonian,
)
from prototype.mrimsrg.densities import compute_densities  # noqa: E402
from prototype.mrimsrg.export_jref import export_reference  # noqa: E402
from prototype.mrimsrg.generator import (  # noqa: E402
    oscillator_quanta_from_orbits,
    spherical_orbit_groups_from_orbits,
    white_generator,
)
from prototype.mrimsrg.normal_order import MRHamiltonian  # noqa: E402
from prototype.mrimsrg.reference_io import load_reference  # noqa: E402
from prototype.mrimsrg.sr_imsrgpp_check import operator_to_mscheme  # noqa: E402


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

# The new explicit entry must reduce to the existing White generator for a
# Slater reference whenever the relaxed Delta-e mask selects the same ph and
# pp-hh channels (He4 in this emax=1 test space).
eta_sr = pyIMSRG.Operator(ms)
eta_sr.SetAntiHermitian()
eta_ncsm_sr = pyIMSRG.Operator(ms)
eta_ncsm_sr.SetAntiHermitian()
generator_sr = pyIMSRG.Generator()
generator_sr.SetType("white")
generator_ncsm_sr = pyIMSRG.Generator()
generator_ncsm_sr.SetType("white-ncsm")
generator_sr.Update(vacuum, eta_sr)
generator_ncsm_sr.Update(vacuum, eta_ncsm_sr)
assert max_operator_difference(eta_sr, eta_ncsm_sr) < 1e-12

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
assert max_operator_difference(sr_commutator, mr_slater_commutator) < 1e-11

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

# The solver owns no hidden density state: without a reference it dispatches
# to the original SR commutator, and after explicit attachment it dispatches
# to the same verified MR entry.
dispatcher_solver = pyIMSRG.IMSRGSolver(y_operator)
dispatcher_sr = dispatcher_solver.EvaluateCommutator(x_operator, y_operator)
assert max_operator_difference(dispatcher_sr, sr_commutator) < 1e-11

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

# A fixed Nmax-truncated reference can be embedded in a larger Hamiltonian
# space without manufacturing correlations on the added orbits.  The C++
# helper must preserve every low-space value, make the added density blocks
# exactly zero/identity, and survive its own binary writer/reader.
target_digest = "1" * 64
be8_ms_emax4 = pyIMSRG.ModelSpace(4, "He4", "He4")
be8_ms_emax4.SetHbarOmega(20.0)
embedded_reference = be8_reference.EmbedInModelSpace(
    be8_ms_emax4, target_digest
)
embedded_reference.Validate(1e-10)
assert embedded_reference.emax == 4
assert embedded_reference.e2max == 8
assert embedded_reference.interaction_sha256 == target_digest
assert embedded_reference.rdm_sha256 == be8_reference.rdm_sha256
assert embedded_reference.wavefunction_sha256 == be8_reference.wavefunction_sha256
assert embedded_reference.MaximumContractionViolation() < 1e-12

source_to_target = {}
for p in range(be8_ms.GetNumberOrbits()):
    orbit = be8_ms.GetOrbit(p)
    target = be8_ms_emax4.GetOrbitIndex(orbit.n, orbit.l, orbit.j2, orbit.tz2)
    source_to_target[p] = target
    assert abs(
        embedded_reference.occupations[target] - be8_reference.occupations[p]
    ) < 1e-15
for p in range(be8_ms.GetNumberOrbits()):
    for q in range(be8_ms.GetNumberOrbits()):
        assert abs(
            embedded_reference.NaturalOrbitTransformation(
                source_to_target[p], source_to_target[q]
            )
            - be8_reference.NaturalOrbitTransformation(p, q)
        ) < 1e-15
for p in range(be8_ms_emax4.GetNumberOrbits()):
    if p in source_to_target.values():
        continue
    assert embedded_reference.occupations[p] == 0.0
    for q in range(be8_ms_emax4.GetNumberOrbits()):
        expected = 1.0 if p == q else 0.0
        assert abs(
            embedded_reference.NaturalOrbitTransformation(p, q) - expected
        ) < 1e-15
        assert abs(
            embedded_reference.NaturalOrbitTransformation(q, p) - expected
        ) < 1e-15
for ch in range(be8_ms.GetNumberTwoBodyChannels()):
    source_channel = be8_ms.GetTwoBodyChannel(ch)
    target_channel_index = be8_ms_emax4.GetTwoBodyChannelIndex(
        source_channel.J, source_channel.parity, source_channel.Tz
    )
    target_channel = be8_ms_emax4.GetTwoBodyChannel(target_channel_index)
    for ibra in range(source_channel.GetNumberKets()):
        bra = source_channel.GetKet(ibra)
        target_bra = target_channel.GetLocalIndex(
            source_to_target[bra.p], source_to_target[bra.q]
        )
        for iket in range(source_channel.GetNumberKets()):
            ket = source_channel.GetKet(iket)
            target_ket = target_channel.GetLocalIndex(
                source_to_target[ket.p], source_to_target[ket.q]
            )
            source_value = be8_reference.Lambda2.GetTBME_norm_chij(
                ch, ch, ibra, iket
            )
            target_value = embedded_reference.Lambda2.GetTBME_norm_chij(
                target_channel_index, target_channel_index, target_bra, target_ket
            )
            assert abs(source_value - target_value) < 1e-15

with tempfile.TemporaryDirectory() as temporary_directory:
    embedded_path = Path(temporary_directory) / "Be8_Nrefmax0_emax4.jref"
    embedded_reference.WriteBinary(str(embedded_path))
    try:
        embedded_reference.WriteBinary(str(embedded_path))
        raise AssertionError("MRReference.WriteBinary overwrote an existing file")
    except RuntimeError as error:
        assert "refusing to overwrite" in str(error)
    roundtrip_ms = pyIMSRG.ModelSpace(4, "He4", "He4")
    roundtrip_ms.SetHbarOmega(20.0)
    roundtrip_ms.SetReferenceOcc(
        pyIMSRG.MRReference.ReadOccupationMap(roundtrip_ms, str(embedded_path))
    )
    roundtrip_reference = pyIMSRG.MRReference.ReadBinary(
        roundtrip_ms, str(embedded_path)
    )
assert roundtrip_reference.interaction_sha256 == target_digest
assert np.max(
    np.abs(
        np.asarray(roundtrip_reference.occupations)
        - np.asarray(embedded_reference.occupations)
    )
) == 0.0
assert (
    roundtrip_reference.NaturalOrbitTransformation
    - embedded_reference.NaturalOrbitTransformation
).Norm() == 0.0
assert (
    roundtrip_reference.Lambda2 - embedded_reference.Lambda2
).Norm() == 0.0

# Compare the complete fractional-occupation White-NCSM generator against the
# Python m-scheme implementation.  The C++ path uses only spherical monopoles
# and J blocks; conversion here is solely the independent test oracle.
be8_unit_test = pyIMSRG.UnitTest(be8_ms)
be8_unit_test.SetRandomSeed(20260825)
be8_hamiltonian = be8_unit_test.RandomOp(be8_ms, 0, 0, 0, 2, +1)
file_to_model = {
    orbit.index: be8_ms.GetOrbitIndex(orbit.n, orbit.l, orbit.j2, orbit.tz2)
    for orbit in be8_orbits
}
be8_one_body_m = np.zeros((be8_data.norb, be8_data.norb))
for p, row_p in enumerate(be8_data.orbits):
    a = int(row_p[0])
    for q, row_q in enumerate(be8_data.orbits):
        b = int(row_q[0])
        if (
            row_p[2] == row_q[2]
            and row_p[3] == row_q[3]
            and row_p[4] == row_q[4]
            and row_p[5] == row_q[5]
        ):
            be8_one_body_m[p, q] = be8_hamiltonian.OneBody(
                file_to_model[a], file_to_model[b]
            )
be8_hamiltonian_blocks = []
for layout_block in be8_blocks:
    channel = be8_ms.GetTwoBodyChannelIndex(
        layout_block.J, layout_block.parity, layout_block.Tz
    )
    two_body_channel = be8_ms.GetTwoBodyChannel(channel)
    matrix = np.zeros_like(layout_block.matrix)
    for ibra, (a_file, b_file) in enumerate(layout_block.pairs):
        local_bra = two_body_channel.GetLocalIndex(
            file_to_model[a_file], file_to_model[b_file]
        )
        for iket, (c_file, d_file) in enumerate(layout_block.pairs):
            local_ket = two_body_channel.GetLocalIndex(
                file_to_model[c_file], file_to_model[d_file]
            )
            matrix[ibra, iket] = be8_hamiltonian.TwoBody.GetTBME_norm_chij(
                channel, channel, local_bra, local_ket
            )
    be8_hamiltonian_blocks.append(
        CoupledBlock(
            layout_block.J,
            layout_block.parity,
            layout_block.Tz,
            layout_block.pairs,
            matrix,
        )
    )
be8_two_body_m = reconstruct_scalar_two_body(
    tuple(be8_hamiltonian_blocks), be8_data.orbits
)
be8_python_eta = white_generator(
    MRHamiltonian(
        be8_hamiltonian.ZeroBody, be8_one_body_m, be8_two_body_m
    ),
    be8_natural.densities,
    oscillator_quanta_from_orbits(be8_data.orbits),
    spherical_orbit_groups=spherical_orbit_groups_from_orbits(be8_data.orbits),
)
be8_cpp_eta = pyIMSRG.Operator(be8_ms)
be8_cpp_eta.SetAntiHermitian()
be8_generator = pyIMSRG.Generator()
be8_generator.SetType("white-ncsm")
be8_generator.Update(be8_hamiltonian, be8_cpp_eta)
be8_dispatcher = pyIMSRG.IMSRGSolver(be8_hamiltonian)
be8_dispatcher.SetMRReference(be8_reference)
be8_dispatched_rhs = be8_dispatcher.EvaluateCommutator(
    be8_cpp_eta, be8_hamiltonian
)
be8_direct_rhs = pyIMSRG.MR_Commutator(
    be8_cpp_eta, be8_hamiltonian, be8_reference
)
assert max_operator_difference(be8_dispatched_rhs, be8_direct_rhs) < 1e-11
for a in be8_orbits:
    for b in be8_orbits:
        if (a.l, a.j2, a.tz2) != (b.l, b.j2, b.tz2):
            continue
        expected = be8_python_eta.one_body[a.substates[0], b.substates[0]]
        actual = be8_cpp_eta.OneBody(file_to_model[a.index], file_to_model[b.index])
        assert abs(actual - expected) < 2e-11
be8_python_eta_blocks = couple_scalar_two_body(
    be8_python_eta.two_body, be8_data.orbits
)
for block in be8_python_eta_blocks:
    channel = be8_ms.GetTwoBodyChannelIndex(block.J, block.parity, block.Tz)
    two_body_channel = be8_ms.GetTwoBodyChannel(channel)
    for ibra, (a_file, b_file) in enumerate(block.pairs):
        local_bra = two_body_channel.GetLocalIndex(
            file_to_model[a_file], file_to_model[b_file]
        )
        for iket, (c_file, d_file) in enumerate(block.pairs):
            local_ket = two_body_channel.GetLocalIndex(
                file_to_model[c_file], file_to_model[d_file]
            )
            actual = be8_cpp_eta.TwoBody.GetTBME_norm_chij(
                channel, channel, local_bra, local_ket
            )
            assert abs(actual - block.matrix[ibra, iket]) < 2e-11

# One complete fixed-step RK4 update must use the MR commutator at all four
# stages.  Reconstruct the same step explicitly and compare every operator
# rank, rather than accepting agreement of the scalar energy alone.
rk4_step = 1e-4


def be8_eta_at(hamiltonian):
    eta = pyIMSRG.Operator(be8_ms)
    eta.SetAntiHermitian()
    generator = pyIMSRG.Generator()
    generator.SetType("white-ncsm")
    generator.Update(hamiltonian, eta)
    return eta


rk4_initial = pyIMSRG.Operator(be8_hamiltonian)
rk4_k1 = pyIMSRG.MR_Commutator(
    be8_eta_at(rk4_initial), rk4_initial, be8_reference
)
rk4_stage2 = rk4_initial + 0.5 * rk4_step * rk4_k1
rk4_k2 = pyIMSRG.MR_Commutator(
    be8_eta_at(rk4_stage2), rk4_stage2, be8_reference
)
rk4_stage3 = rk4_initial + 0.5 * rk4_step * rk4_k2
rk4_k3 = pyIMSRG.MR_Commutator(
    be8_eta_at(rk4_stage3), rk4_stage3, be8_reference
)
rk4_stage4 = rk4_initial + rk4_step * rk4_k3
rk4_k4 = pyIMSRG.MR_Commutator(
    be8_eta_at(rk4_stage4), rk4_stage4, be8_reference
)
rk4_expected = rk4_initial + rk4_step / 6.0 * (
    rk4_k1 + 2.0 * rk4_k2 + 2.0 * rk4_k3 + rk4_k4
)
rk4_solver = pyIMSRG.IMSRGSolver(rk4_initial)
rk4_solver.SetMRReference(be8_reference)
rk4_solver.SetGenerator("white-ncsm")
rk4_solver.SetMethod("flow_RK4")
rk4_solver.SetEtaCriterion(0.0)
rk4_solver.SetDsmax(rk4_step)
rk4_solver.SetSmax(rk4_step)
rk4_solver.Solve()
assert max_operator_difference(rk4_solver.GetH_s(), rk4_expected) < 3e-11

# A correlated closed-shell Nrefmax=2 reference exercises the non-identity
# temporary natural-orbit rotation, including radial mixing in fixed (l,j,tz)
# blocks.  It must survive the same bridge without being mistaken for HO basis.
he4_nref2_path = repository_root / "prototype/mrimsrg/data/He4_Nrefmax2"
he4_nref2_data = load_reference(he4_nref2_path)
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

# The reusable J-scheme 0/1/2B basis transformation must agree tensor by
# tensor with the independent m-scheme covariance helper.  Use the physical
# He4 Nrefmax=2 natural orbitals so this exercises actual radial mixing rather
# than an artificial permutation.  Production stores only the J-scheme result;
# the dense tensors below exist solely as a small-space test oracle.
he4_transform_test = pyIMSRG.UnitTest(he4_nref2_ms)
he4_transform_test.SetRandomSeed(424242)
he4_ho_operator = he4_transform_test.RandomOp(
    he4_nref2_ms, 0, 0, 0, 2, +1
)
he4_U_j = np.asarray(
    [
        [
            he4_nref2_reference.NaturalOrbitTransformation(i, j)
            for j in range(he4_nref2_ms.GetNumberOrbits())
        ]
        for i in range(he4_nref2_ms.GetNumberOrbits())
    ]
)
he4_nat_operator = he4_ho_operator.TransformOneAndTwoBody(
    he4_nref2_reference.NaturalOrbitTransformation
)
he4_round_trip = he4_nat_operator.TransformOneAndTwoBody(
    he4_nref2_reference.NaturalOrbitTransformation.t()
)
assert max_operator_difference(he4_round_trip, he4_ho_operator) < 2e-11

he4_U_m = np.zeros((he4_nref2_data.norb, he4_nref2_data.norb))
for p, old_orbit in enumerate(he4_nref2_data.orbits):
    old_parent, _, old_l, old_j2, old_m2, old_tz2 = (
        int(value) for value in old_orbit
    )
    for q, new_orbit in enumerate(he4_nref2_data.orbits):
        new_parent, _, new_l, new_j2, new_m2, new_tz2 = (
            int(value) for value in new_orbit
        )
        if (old_l, old_j2, old_m2, old_tz2) == (
            new_l,
            new_j2,
            new_m2,
            new_tz2,
        ):
            he4_U_m[p, q] = he4_U_j[old_parent, new_parent]
he4_ho_m = operator_to_mscheme(he4_ho_operator, he4_nref2_data.orbits)
he4_expected_nat_m = transform_hamiltonian(
    he4_ho_m, he4_U_m, to_natural=True
)
he4_actual_nat_m = operator_to_mscheme(
    he4_nat_operator, he4_nref2_data.orbits
)
np.testing.assert_allclose(
    he4_actual_nat_m.one_body,
    he4_expected_nat_m.one_body,
    atol=2e-11,
    rtol=0,
)
np.testing.assert_allclose(
    he4_actual_nat_m.two_body,
    he4_expected_nat_m.two_body,
    atol=3e-11,
    rtol=0,
)

# The lossless J-coupled bridge used by the production driver must preserve
# every stored rank exactly enough for the 1e-10 algebra gates.  no2bpack is
# intentionally not used here because its established payload is float32.
with tempfile.TemporaryDirectory() as temporary_directory:
    j64_path = Path(temporary_directory) / "round_trip.jcoupled64"
    j64_writer = pyIMSRG.ReadWrite()
    j64_writer.Write_jcoupled64(str(j64_path), he4_ho_operator)
    j64_recovered = pyIMSRG.Operator(he4_nref2_ms)
    j64_reader = pyIMSRG.ReadWrite()
    j64_reader.Read_jcoupled64(str(j64_path), j64_recovered)
assert max_operator_difference(j64_recovered, he4_ho_operator) < 1e-13

print("MRReference tests passed")
