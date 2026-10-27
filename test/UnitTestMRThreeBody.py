#!/usr/bin/env python3

"""MR-NO2B normal ordering of 3N interactions and the lambda3 terms.

Uses the correlated He4 Nrefmax=2 reference (lambda2, lambda3 != 0) and checks
everything against quantities computed directly from its m-scheme
wavefunction, so no part of the J-coupled machinery checks itself:

1. jref version 2: lambda3 and the reference energy survive the C++
   reader/writer; the version-1 export is byte-identical to before.
2. MR-NO2B (MR_NormalOrderThreeBody) with an exactly antisymmetric random
   3N operator: E3 = <Psi|V3|Psi> including the lambda3 term, f3 and Gamma3
   against m-scheme contractions with the wavefunction densities.  With
   the N3LO lnl 3N file (emax6, hw16; only used as a source of realistic
   matrix elements) the no2b and full files agree to their single-precision
   storage, and only the full file includes the lambda3 term.
3. The zero-body lambda3 term of the MR commutator: MR_Commutator(X, Y)
   zero body = <Psi|[X, Y]|Psi> exactly; without the term it is off by it.
4. imsrg++ with 3bme (no2b and full): log lines, metadata, and the initial
   zero body equal to the library MR-NO2B plus the NN part.
"""

import argparse
from collections import defaultdict
from functools import lru_cache
import hashlib
from itertools import combinations, product
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import pyIMSRG

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))
sys.path.insert(0, str(REPOSITORY / "test"))

from prototype.mrimsrg.export_jref import export_reference  # noqa: E402
from prototype.mrimsrg.lambda3 import annihilation_matrix, cg2  # noqa: E402
from prototype.mrimsrg.reference_io import load_reference  # noqa: E402
from UnitTestMRDriver import vacuum_operator  # noqa: E402
from UnitTestMRTensor import expand_one_body, expand_two_body  # noqa: E402

FIXTURE = REPOSITORY / "prototype/mrimsrg/data/He4_Nrefmax2"
THREE_BODY = REPOSITORY / "prototype/mrimsrg/data/N3LO_lnl_3N"
FULL_FILE = THREE_BODY / "ThBME_full_srg2.0_N3LO_EM500_lnl_N3max40_hw16_emax6_e2max12_e3max6.me3j.gz"
NO2B_FILE = THREE_BODY / "ThBME_NO2B_srg2.0_N3LO_EM500_lnl_N3max40_hw16_emax6_e2max12_e3max6.me3j.gz"
V1_SHA256 = "1c9934aa92b50cc55122f15434b0093a860d99ae50a745d89848f76f1f608898"


def load(path, e3max=6):
    modelspace = pyIMSRG.ModelSpace(2, "He4", "He4")
    modelspace.SetE3max(e3max)
    modelspace.SetReferenceOcc(pyIMSRG.MRReference.ReadOccupationMap(modelspace, str(path)))
    return modelspace, pyIMSRG.MRReference.ReadBinary(modelspace, str(path), 1e-10)


class Wavefunction:
    """The fixture determinants in the oscillator m-scheme, ModelSpace orbit labels."""

    def __init__(self, modelspace):
        orbits = np.load(FIXTURE / "orbits.npy")  # jindex, n, l, j2, m2, tz2
        self.dets = np.load(FIXTURE / "determinants.npy").astype(np.uint8)
        self.coef = np.load(FIXTURE / "coefficients.npy").astype(float)
        label = {}
        for i in range(modelspace.GetNumberOrbits()):
            o = modelspace.GetOrbit(i)
            label[(o.n, o.l, o.j2, o.tz2)] = i
        self.mstates = [(label[(int(r[1]), int(r[2]), int(r[3]), int(r[5]))], int(r[4])) for r in orbits]
        self.size = len(self.mstates)
        self.index = {state: k for k, state in enumerate(self.mstates)}
        self.j2 = [modelspace.GetOrbit(i).j2 for i in range(modelspace.GetNumberOrbits())]

    def density(self, rank):
        columns = list(combinations(range(self.size), rank))
        A = annihilation_matrix(self.dets, self.coef, columns)
        return columns, (A.T @ A)


class MschemeThreeBody:
    """Antisymmetric m-scheme <abc|V|def> from ThreeBodyME::GetME_pn (UnitTest convention)."""

    def __init__(self, operator, wavefunction):
        self.three_body = operator.ThreeBody
        self.wf = wavefunction
        self.cache = {}

    def coupled(self, *key):
        if key not in self.cache:
            self.cache[key] = self.three_body.GetME_pn(*key)
        return self.cache[key]

    def __call__(self, x, y, z, u, v, w):
        st, j2 = self.wf.mstates, self.wf.j2
        (a, ma), (b, mb), (c, mc) = st[x], st[y], st[z]
        (d, md), (e, me), (f, mf) = st[u], st[v], st[w]
        M = ma + mb + mc
        if M != md + me + mf:
            return 0.0
        total = 0.0
        for J12 in range(abs(j2[a] - j2[b]) // 2, (j2[a] + j2[b]) // 2 + 1):
            if a == b and J12 % 2:
                continue
            c1 = cg2(j2[a], ma, j2[b], mb, 2 * J12, ma + mb)
            if c1 == 0.0:
                continue
            for J45 in range(abs(j2[d] - j2[e]) // 2, (j2[d] + j2[e]) // 2 + 1):
                if d == e and J45 % 2:
                    continue
                c2 = cg2(j2[d], md, j2[e], me, 2 * J45, md + me)
                if c2 == 0.0:
                    continue
                low = max(abs(2 * J12 - j2[c]), abs(2 * J45 - j2[f]), abs(M))
                for twoJ in range(low, min(2 * J12 + j2[c], 2 * J45 + j2[f]) + 1, 2):
                    c3 = cg2(2 * J12, ma + mb, j2[c], mc, twoJ, M)
                    c4 = cg2(2 * J45, md + me, j2[f], mf, twoJ, M)
                    if c3 and c4:
                        total += c1 * c2 * c3 * c4 * self.coupled(J12, J45, twoJ, a, b, c, d, e, f)
        return total


def matrix(arma, n):
    return np.array([[arma(i, j) for j in range(n)] for i in range(n)])


def check_reference_files(root):
    v1 = root / "v1.jref"
    v2 = root / "v2.jref"
    export_reference(FIXTURE, v1)
    summary = export_reference(FIXTURE, v2, lambda3=True)
    assert hashlib.sha256(v1.read_bytes()).hexdigest() == V1_SHA256, "version-1 export changed"
    _, ref1 = load(v1)
    _, ref2 = load(v2)
    assert not ref1.has_lambda3 and not ref1.has_reference_energy
    assert ref2.has_lambda3 and ref2.Lambda3Norm() > 1e-2, ref2.Lambda3Norm()
    reference_energy = float(load_reference(FIXTURE).metadata["reference_energy"])
    assert ref2.has_reference_energy and abs(ref2.reference_energy - reference_energy) < 1e-12
    copy = root / "v2_copy.jref"
    ref2.WriteBinary(str(copy))
    _, ref3 = load(copy)
    blocks2, blocks3 = ref2.Lambda3Blocks(), ref3.Lambda3Blocks()
    assert len(blocks2) == len(blocks3) == summary["lambda3 channels"]
    for b2, b3 in zip(blocks2, blocks3):
        assert b2[:4] == b3[:4] and np.array_equal(np.array(b2[4]), np.array(b3[4]))
    print(f"jref v2: {summary['lambda3 channels']} lambda3 channels, |lambda3| = {ref2.Lambda3Norm():.4e}, "
          f"round trip exact; version 1 unchanged")
    return v1, v2


def check_mr_no2b(v2):
    modelspace, reference = load(v2)
    norbits = modelspace.GetNumberOrbits()
    C = matrix(reference.NaturalOrbitTransformation, norbits)
    wf = Wavefunction(modelspace)

    unit = pyIMSRG.UnitTest(modelspace)
    unit.SetRandomSeed(7)
    random_3b = unit.RandomOp(modelspace, 0, 0, 0, 3, +1)
    N3, report = pyIMSRG.MR_NormalOrderThreeBody(random_3b, reference, 1e-10)
    assert report.storage == "full" and report.lambda3_term_included
    V = MschemeThreeBody(random_3b, wf)

    triples, gamma3 = wf.density(3)
    gamma3 = gamma3.tocoo()
    exact = sum(g * V(*triples[i], *triples[j]) for i, j, g in zip(gamma3.row, gamma3.col, gamma3.data))
    e3_error = abs(N3.ZeroBody - exact)
    print(f"MR-NO2B random 3N: E3 {N3.ZeroBody:.12f} vs <Psi|V3|Psi> {exact:.12f} (diff {e3_error:.2e}); "
          f"lambda3 term {report.e3_lambda3:.6f}")
    assert abs(report.e3_lambda3) > 1e-3 and e3_error < 1e-10

    # f3 = 1/4 sum v gamma2, oscillator basis
    pairs, gamma2 = wf.density(2)
    gamma2 = gamma2.tocoo()
    f_oscillator = C @ matrix(N3.OneBody, norbits) @ C.T
    f_error = 0.0
    for p, q in product(range(norbits), repeat=2):
        op, oq = modelspace.GetOrbit(p), modelspace.GetOrbit(q)
        if p > q or (op.l, op.j2, op.tz2) != (oq.l, oq.j2, oq.tz2):
            continue
        x, y = wf.index[(p, op.j2)], wf.index[(q, oq.j2)]
        value = sum(g * V(x, *pairs[i], y, *pairs[j])
                    for i, j, g in zip(gamma2.row, gamma2.col, gamma2.data)
                    if x not in pairs[i] and y not in pairs[j])
        f_error = max(f_error, abs(value - f_oscillator[p, q]))

    # Gamma3 = sum v gamma1, coupled to J, sampled elements
    gamma1 = wf.density(1)[1].toarray()
    N3_oscillator = N3.TransformOneAndTwoBody(reference.NaturalOrbitTransformation.t(), 1e-10)
    rng = np.random.default_rng(3)
    g_error = 0.0
    j2 = wf.j2
    for _ in range(30):
        channel = modelspace.GetTwoBodyChannel(int(rng.integers(modelspace.GetNumberTwoBodyChannels())))
        nkets = channel.GetNumberKets()
        if nkets == 0:
            continue
        bra, ket = channel.GetKet(int(rng.integers(nkets))), channel.GetKet(int(rng.integers(nkets)))
        J = channel.J
        p, q, r, s = bra.p, bra.q, ket.p, ket.q
        value = 0.0
        for mp in range(-j2[p], j2[p] + 1, 2):
            mq = -mp
            if abs(mq) > j2[q] or (p == q and mp == mq):
                continue
            cb = cg2(j2[p], mp, j2[q], mq, 2 * J, 0)
            for mr in range(-j2[r], j2[r] + 1, 2):
                ms_ = -mr
                if abs(ms_) > j2[s] or (r == s and mr == ms_) or cb == 0.0:
                    continue
                ck = cg2(j2[r], mr, j2[s], ms_, 2 * J, 0)
                if ck == 0.0:
                    continue
                X1, X2 = wf.index[(p, mp)], wf.index[(q, mq)]
                Y1, Y2 = wf.index[(r, mr)], wf.index[(s, ms_)]
                for xx, yy in zip(*np.nonzero(gamma1)):
                    if xx in (X1, X2) or yy in (Y1, Y2):
                        continue
                    value += cb * ck * gamma1[xx, yy] * V(X1, X2, int(xx), Y1, Y2, int(yy))
        g_error = max(g_error, abs(value - N3_oscillator.TwoBody.GetTBME_J(J, J, p, q, r, s)))
    print(f"MR-NO2B random 3N: f3 max diff {f_error:.2e}, Gamma3 max diff {g_error:.2e}")
    assert f_error < 1e-10 and g_error < 1e-10

    # Real 3N file in its two formats.
    full = pyIMSRG.Operator(modelspace, 0, 0, 0, 3)
    pyIMSRG.ReadWrite().Read_Darmstadt_3body(str(FULL_FILE), full, 6, 12, 6)
    no2b = pyIMSRG.Operator(modelspace, 0, 0, 0, 3)
    no2b.ThreeBody.SetMode("no2b")
    no2b.ThreeBody.ReadFile([str(NO2B_FILE)], [6, 12, 6, 6])
    N_full, r_full = pyIMSRG.MR_NormalOrderThreeBody(full, reference, 1e-10)
    N_no2b, r_no2b = pyIMSRG.MR_NormalOrderThreeBody(no2b, reference, 1e-10)
    assert r_full.storage == "full" and r_full.lambda3_term_included
    assert r_no2b.storage == "no2b" and not r_no2b.lambda3_term_included and "no2b" in r_no2b.lambda3_reason
    one_diff = np.abs(matrix(N_full.OneBody, norbits) - matrix(N_no2b.OneBody, norbits)).max()
    two_diff = (N_full - N_no2b).TwoBodyNorm() / N_full.TwoBodyNorm()
    zero_diff = abs(N_full.ZeroBody - r_full.e3_lambda3 - N_no2b.ZeroBody)
    V_real = MschemeThreeBody(full, wf)
    exact_real = sum(g * V_real(*triples[i], *triples[j])
                     for i, j, g in zip(gamma3.row, gamma3.col, gamma3.data))
    print(f"N3LO lnl 3N: full vs no2b 1B {one_diff:.1e}, 2B {two_diff:.1e} (relative), 0B {zero_diff:.1e}; "
          f"full E3 vs <Psi|V3|Psi> {abs(N_full.ZeroBody - exact_real):.1e} (lambda3 term {r_full.e3_lambda3:.4f})")
    # Single-precision storage; the Darmstadt (ab)J12 c states are overcomplete,
    # so orderings that should agree differ at the float level.
    assert one_diff < 1e-6 and two_diff < 1e-6 and zero_diff < 1e-6
    assert abs(N_full.ZeroBody - exact_real) < 1e-7
    return modelspace, reference


def apply(operator, state, size):
    _, one, two = operator
    out = defaultdict(float)
    nonzero = [(p, q, one[p, q]) for p in range(size) for q in range(size) if one[p, q] != 0.0]
    for det, c in state.items():
        occupied = [k for k in range(size) if det >> k & 1]
        for p, q, v in nonzero:
            if not det >> q & 1:
                continue
            d1 = det ^ (1 << q)
            if d1 >> p & 1:
                continue
            sign = bin(det & ((1 << q) - 1)).count("1") + bin(d1 & ((1 << p) - 1)).count("1")
            out[d1 | (1 << p)] += (-1) ** sign * v * c
        for r, s in combinations(occupied, 2):
            d1 = det ^ (1 << r)
            d2 = d1 ^ (1 << s)
            base = bin(det & ((1 << r) - 1)).count("1") + bin(d1 & ((1 << s) - 1)).count("1")
            free = [k for k in range(size) if not d2 >> k & 1]
            block = two[:, :, r, s]
            for p, q in combinations(free, 2):
                v = block[p, q]
                if v == 0.0:
                    continue
                d3 = d2 | (1 << q)
                sign = base + bin(d2 & ((1 << q) - 1)).count("1") + bin(d3 & ((1 << p) - 1)).count("1")
                out[d3 | (1 << p)] += (-1) ** sign * v * c   # a^dag_p a^dag_q a_s a_r
    return out


def check_commutator(modelspace, reference):
    unit = pyIMSRG.UnitTest(modelspace)
    unit.SetRandomSeed(11)
    X = unit.RandomOp(modelspace, 0, 0, 0, 2, -1)
    Y = unit.RandomOp(modelspace, 0, 0, 0, 2, +1)
    Z = pyIMSRG.MR_Commutator(X, Y, reference)
    lambda3_term = pyIMSRG.MR_comm220_lambda3(X, Y, reference)

    wf = Wavefunction(modelspace)
    back = reference.NaturalOrbitTransformation.t()

    def mscheme(op):
        vacuum = reference.UndoNormalOrder(op).TransformOneAndTwoBody(back, 1e-10)
        return (vacuum.ZeroBody, expand_one_body(vacuum, modelspace, wf.mstates, 0, 0),
                expand_two_body(vacuum.TwoBody, modelspace, wf.mstates, 0, 0))

    weights = 1 << np.arange(wf.size, dtype=np.uint64)
    psi = {int(row @ weights): c for row, c in zip(wf.dets.astype(np.uint64), wf.coef)}
    x_psi = apply(mscheme(X), psi, wf.size)
    y_psi = apply(mscheme(Y), psi, wf.size)
    exact = -2.0 * sum(x_psi[k] * y_psi.get(k, 0.0) for k in x_psi)  # X anti-Hermitian, Y Hermitian
    error = abs(Z.ZeroBody - exact)
    print(f"MR commutator zero body {Z.ZeroBody:.12f} vs <Psi|[X,Y]|Psi> {exact:.12f} (diff {error:.2e}); "
          f"lambda3 term {lambda3_term:.6f}")
    assert abs(lambda3_term) > 1e-2 and error < 1e-9 * max(1.0, abs(exact))


def check_driver(executable, root, v1, v2):
    data = load_reference(FIXTURE)
    modelspace_nn = pyIMSRG.ModelSpace(2, "He4", "He4")
    modelspace_nn.SetHbarOmega(20.0)
    nn_file = root / "He4_bare.jcoupled64"
    pyIMSRG.ReadWrite().Write_jcoupled64(str(nn_file), vacuum_operator(modelspace_nn, data))

    for label, reference_file, three_body, three_body_type in (
        ("no2b", v2, NO2B_FILE, "no2b"),
        ("full", v2, FULL_FILE, "full"),
        ("full_v1", v1, FULL_FILE, "full"),
    ):
        output = root / f"H_{label}.jcoupled64"
        command = [
            str(executable), f"2bme={nn_file}", "fmt2=jcoupled64",
            f"3bme={three_body}", f"3bme_type={three_body_type}",
            "e3max=6", "file3e1max=6", "file3e2max=12", "file3e3max=6",
            "reference=He4", "valence_space=He4", "A=4", "hw=20", "emax=2",
            "basis=oscillator", "method=magnus_adaptive", "core_generator=white-ncsm",
            "denominator_partitioning=Epstein_Nesbet", "BetaCM=0", "nucleon_mass_correction=false",
            f"mr_reference_file={reference_file}", "smax=0.05", "eta_criterion=1e-12",
            f"scratch={root}", f"flowfile={root / ('flow_' + label + '.dat')}",
            f"intfile={root / ('int_' + label)}", f"write_H_jcoupled64={output}",
        ]
        environment = dict(os.environ, OMP_NUM_THREADS="2")
        completed = subprocess.run(command, text=True, capture_output=True, env=environment, check=False)
        if completed.returncode != 0:
            print(completed.stdout[-4000:])
            print(completed.stderr[-4000:], file=sys.stderr)
            raise RuntimeError(f"MR driver with 3N ({label}) failed")
        log = completed.stdout
        meta = json.loads(Path(str(output) + ".mrmeta.json").read_text())

        # Library reference for the initial zero body: NN part + MR-NO2B 3N.
        modelspace, reference = load(reference_file)
        hbare = pyIMSRG.Operator(modelspace, 0, 0, 0, 3)
        rw = pyIMSRG.ReadWrite()
        rw.Read_jcoupled64(str(nn_file), hbare)
        if three_body_type == "full":
            rw.Read_Darmstadt_3body(str(three_body), hbare, 6, 12, 6)
        else:
            hbare.ThreeBody.SetMode("no2b")
            hbare.ThreeBody.ReadFile([str(three_body)], [6, 12, 6, 6])
        N3, report = pyIMSRG.MR_NormalOrderThreeBody(hbare, reference, 1e-10)
        nn = pyIMSRG.Operator(modelspace, 0, 0, 0, 2)
        rw.Read_jcoupled64(str(nn_file), nn)
        expected = reference.NormalOrder(nn.TransformOneAndTwoBody(reference.NaturalOrbitTransformation, 1e-10))
        expected_zero = expected.ZeroBody + N3.ZeroBody
        line = next(x for x in log.splitlines() if x.startswith("MR initial natural-basis Hamiltonian: E="))
        driver_zero = float(line.split("E=")[1].split()[0])
        assert abs(driver_zero - expected_zero) < 1e-9, (label, driver_zero, expected_zero)

        lambda3_present = reference_file == v2
        assert f"MR reference lambda3: {'present' if lambda3_present else 'absent'}" in log
        assert meta["three_body"]["input"] == three_body_type
        assert meta["lambda3"]["in_reference"] == lambda3_present
        assert meta["lambda3"]["commutator_zero_body_term"] == lambda3_present
        included = three_body_type == "full" and lambda3_present
        assert meta["three_body"]["lambda3_zero_body_term"].startswith("included" if included else "omitted")
        assert ("lambda3 x V3 zero-body term included" in log) == included
        if lambda3_present:
            assert "initial E - E_ref" in log
        print(f"imsrg++ 3bme_type={three_body_type} reference {'v2' if lambda3_present else 'v1'}: "
              f"initial E {driver_zero:.10f} (library {expected_zero:.10f}), "
              f"lambda3xV3 {'included' if included else 'omitted'}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("imsrg_executable", type=Path)
    args = parser.parse_args()
    for path in (FULL_FILE, NO2B_FILE):
        if not path.is_file():
            raise SystemExit(f"missing 3N test file {path} (see doc/guide/测试.md)")
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        v1, v2 = check_reference_files(root)
        modelspace, reference = check_mr_no2b(v2)
        check_commutator(modelspace, reference)
        check_driver(args.imsrg_executable, root, v1, v2)
    print("MR-NO2B and lambda3 checks passed")


if __name__ == "__main__":
    main()
