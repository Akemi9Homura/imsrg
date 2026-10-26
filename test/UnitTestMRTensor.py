#!/usr/bin/env python3

"""MR-IMSRG(2) commutator and MR normal ordering for tensor operators.

With the correlated He4 Nrefmax=2 reference (lambda2 != 0, fractional
natural occupations), the C++ [Omega, O] for a scalar anti-Hermitian Omega
and a Hermitian rank-L tensor O (L = 1 odd and even parity, L = 2) is
expanded to an m-scheme component and compared with the dense m-scheme
oracle (Hergert Eqs. (49)-(51), lambda3 = 0).  This covers the SR scalar-
tensor terms with fractional occupations and the tensor lambda2 one-body
terms, which vanish in the single-reference limit.  MRReference.NormalOrder
and UndoNormalOrder of a tensor are checked the same way.

Reduced matrix elements use the imsrg++ (Edmonds) Wigner-Eckart convention
<J M|T^L_mu|J' M'> = (-1)^(J-M) (J L J'; -M mu M') <J||T||J'>.
"""

from functools import lru_cache
from pathlib import Path
import sys
import tempfile

import numpy as np
import pyIMSRG
from sympy import Rational
from sympy.physics.wigner import clebsch_gordan, wigner_3j

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))

from prototype.mrimsrg.commutator import commutator_contributions  # noqa: E402
from prototype.mrimsrg.densities import Densities  # noqa: E402
from prototype.mrimsrg.export_jref import export_reference  # noqa: E402
from prototype.mrimsrg.normal_order import MRHamiltonian, VacuumHamiltonian, normal_order  # noqa: E402


@lru_cache(maxsize=None)
def cg(j1, m1, j2, m2, J, M):  # doubled arguments
    if m1 + m2 != M or abs(M) > J:
        return 0.0
    return float(clebsch_gordan(*(Rational(x, 2) for x in (j1, j2, J, m1, m2, M))))


@lru_cache(maxsize=None)
def wigner_eckart(J, M, L, mu, Jp, Mp):  # doubled; L = 0 means ordinary ME
    if L == 0:
        return 1.0 if (J == Jp and M == Mp) else 0.0
    if -M + mu + Mp != 0:
        return 0.0
    phase = -1.0 if ((J - M) // 2) % 2 else 1.0
    return phase * float(wigner_3j(*(Rational(x, 2) for x in (J, L, Jp, -M, mu, Mp))))


def m_states(modelspace):
    states = []
    for index in range(modelspace.GetNumberOrbits()):
        orbit = modelspace.GetOrbit(index)
        for m2 in range(-orbit.j2, orbit.j2 + 1, 2):
            states.append((index, m2))
    return states


def expand_one_body(op, modelspace, states, L2, mu2):
    j2 = [modelspace.GetOrbit(i).j2 for i in range(modelspace.GetNumberOrbits())]
    size = len(states)
    out = np.zeros((size, size))
    for p, (a, ma) in enumerate(states):
        for q, (b, mb) in enumerate(states):
            weight = wigner_eckart(j2[a], ma, L2, mu2, j2[b], mb)
            if weight:
                out[p, q] = weight * op.GetOneBody(a, b)
    return out


def expand_two_body(two_body, modelspace, states, L2, mu2):
    """Unnormalized coupled TwoBodyME::GetTBME_J -> antisymmetric m-scheme."""
    size = len(states)
    out = np.zeros((size, size, size, size))
    norbits = modelspace.GetNumberOrbits()
    j2 = [modelspace.GetOrbit(i).j2 for i in range(norbits)]
    by_orbit = {}
    for p, (a, m) in enumerate(states):
        by_orbit.setdefault(a, []).append((p, m))
    for a in range(norbits):
        for b in range(norbits):
            for c in range(norbits):
                for d in range(norbits):
                    for Jb in range(abs(j2[a] - j2[b]) // 2, (j2[a] + j2[b]) // 2 + 1):
                        for Jk in range(abs(j2[c] - j2[d]) // 2, (j2[c] + j2[d]) // 2 + 1):
                            if L2 == 0 and Jb != Jk:
                                continue
                            if not abs(2 * Jb - 2 * Jk) <= L2 <= 2 * (Jb + Jk):
                                continue
                            value = two_body.GetTBME_J(Jb, Jk, a, b, c, d)
                            if value == 0.0:
                                continue
                            for p, ma in by_orbit[a]:
                                for q, mb in by_orbit[b]:
                                    cb = cg(j2[a], ma, j2[b], mb, 2 * Jb, ma + mb)
                                    if cb == 0.0:
                                        continue
                                    for r, mc in by_orbit[c]:
                                        for s, md in by_orbit[d]:
                                            if mc + md != ma + mb - mu2:
                                                continue
                                            ck = cg(j2[c], mc, j2[d], md, 2 * Jk, mc + md)
                                            if ck == 0.0:
                                                continue
                                            w = wigner_eckart(2 * Jb, ma + mb, L2, mu2, 2 * Jk, mc + md)
                                            out[p, q, r, s] += cb * ck * w * value
    return out


def expand(op, modelspace, states, mu2):
    L2 = 2 * op.GetJRank()
    return (expand_one_body(op, modelspace, states, L2, mu2),
            expand_two_body(op.TwoBody, modelspace, states, L2, mu2))


def main():
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "He4_Nrefmax2.jref"
        export_reference(REPOSITORY / "prototype/mrimsrg/data/He4_Nrefmax2", path)
        modelspace = pyIMSRG.ModelSpace(2, "He4", "He4")
        modelspace.SetReferenceOcc(pyIMSRG.MRReference.ReadOccupationMap(modelspace, str(path)))
        reference = pyIMSRG.MRReference.ReadBinary(modelspace, str(path))
    assert reference.Lambda2.Norm() > 1e-3

    states = m_states(modelspace)
    occupations = np.array([modelspace.GetOrbit(a).occ for a, _ in states])
    lambda2 = expand_two_body(reference.Lambda2, modelspace, states, 0, 0)
    densities = Densities(np.diag(occupations), np.zeros_like(lambda2), lambda2)

    test = pyIMSRG.UnitTest(modelspace)
    test.SetRandomSeed(20261026)
    omega = test.RandomOp(modelspace, 0, 0, 0, 2, -1)
    x1, x2 = expand(omega, modelspace, states, 0)

    worst = 0.0
    for L, parity, mu2 in ((1, 1, 0), (1, 0, 2), (2, 0, 2)):
        op = test.RandomOp(modelspace, L, 0, parity, 2, +1)
        y1, y2 = expand(op, modelspace, states, mu2)
        scale = max(np.abs(y1).max(), np.abs(y2).max())

        # Commutator.
        result = pyIMSRG.MR_Commutator(omega, op, reference)
        z1, z2 = expand(result, modelspace, states, mu2)
        terms = commutator_contributions(
            MRHamiltonian(0.0, x1, x2), MRHamiltonian(0.0, y1, y2), densities
        )
        e1 = sum(t.one_body for t in terms.values())
        e2 = sum(t.two_body for t in terms.values())
        lam = terms["mr_lambda2_one_body"].one_body
        one_error = float(np.abs(z1 - e1).max())
        two_error = float(np.abs(z2 - e2).max())
        print(
            f"L={L} parity={parity}: [Omega,O] one-body {one_error:.2e} two-body {two_error:.2e}"
            f"  (|lambda2 one-body| {np.abs(lam).max():.2e}, |O| {scale:.2e})"
        )
        assert np.abs(lam).max() > 1e-3
        worst = max(worst, one_error, two_error)

        # MR normal ordering and its inverse.
        normal = reference.NormalOrder(op)
        n1, n2 = expand(normal, modelspace, states, mu2)
        expected = normal_order(VacuumHamiltonian(0.0, y1, y2), densities)
        no_error = max(float(np.abs(n1 - expected.one_body).max()),
                       float(np.abs(n2 - expected.two_body).max()))
        back = reference.UndoNormalOrder(normal)
        b1, b2 = expand(back, modelspace, states, mu2)
        undo_error = max(float(np.abs(b1 - y1).max()), float(np.abs(b2 - y2).max()))
        print(f"           NormalOrder {no_error:.2e}  UndoNormalOrder round trip {undo_error:.2e}")
        worst = max(worst, no_error, undo_error)

    assert worst < 1e-10, worst
    print("MR tensor commutator and normal ordering match the m-scheme oracle")


if __name__ == "__main__":
    main()
