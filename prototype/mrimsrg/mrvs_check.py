#!/usr/bin/env python3
"""MR-VS check: VS-IMSRG(2) with the ensemble reference vs with a correlated (MR) reference.

Be8, emax2, NNLOopt hw20, NN only, BetaCM = 0, valence space p shell above a He4 core, default VS
generators (atan core step, shell-model-atan valence step) and method=magnus in every run.  Nmax 0 of
emax2 Be8 is exactly the He4 core times the p shell, so BIGSTICK in Nmax 0 both diagonalizes the valence
space of an exported vacuum H(s) and supplies the valence-space reference (bigstick_jref HF-basis chain).

  sr        ensemble reference, basis=HF: vacuum H(s) and Rp2(s), .snt, C.mat, occupations
  gate      MR-VS with the ensemble occupations, lambda2 = 0 and C = C_HF: must reproduce sr
            (run with mr_validation_tolerance=0.3: lambda2 = 0 is not N-representable)
  mr<k>     MR-VS with the reference from the Nmax 0 ground state of the previous H(s) (k = 1: sr)

Each H(s) is diagonalized in Nmax 0 (E, 2+, <Rp2>); the .snt route (core energy + KSHELL) must give the
same energies, which checks the MR -> core re-normal-ordering.  Exact Be8 values: BIGSTICK Nmax 12
(full emax2), doc/records/工作记录.md, 2026-10-27 Be8 gate.

  mrvs_check.py --workdir DIR --fixture-dir DATA [--iterations 3]
"""

import argparse
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
from types import SimpleNamespace

import numpy as np

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "test"), str(REPOSITORY / "build/src")]

from prototype.mrimsrg.bigstick_jref import basis_transformation_for, convert, read_no2bpack_orbits  # noqa: E402
from prototype.mrimsrg.diagnose_collapse import ground_state  # noqa: E402
from prototype.mrimsrg.export_jref import write_jref  # noqa: E402
from prototype.mrimsrg.jref_format import MAGIC  # noqa: E402
from prototype.mrimsrg.postdiag_bigstick import expectation_values, parse_res  # noqa: E402
from prototype.mrimsrg.reference_io import load_reference  # noqa: E402

EXACT = {"E": -26.325, "Rp2": 4.209, "Ex2+": 3.340}
COMMON = ["fmt2=jcoupled64", "3bme=none", "reference=Be8", "valence_space=p-shell", "A=8", "hw=20", "emax=2",
          "BetaCM=0", "nucleon_mass_correction=false", "method=magnus", "Operators=Rp2"]


def read_jref(path):
    """Header, orbits, occupations, transformation and lambda2 blocks of a jref (the tail is not needed)."""
    data = Path(path).read_bytes()
    if not data.startswith(MAGIC):
        raise ValueError(f"{path}: not a jref")
    pos = len(MAGIC) + 8
    *header, norbits, nchannels = struct.unpack_from("<7idQQ", data, pos)
    pos += struct.calcsize("<7idQQ")
    digests = tuple(data[pos + 64 * k:pos + 64 * (k + 1)].decode() for k in range(3))
    pos += 192
    orbits, occupations = [], []
    for _ in range(norbits):
        *orbit, occ = struct.unpack_from("<5id", data, pos)
        pos += struct.calcsize("<5id")
        orbits.append(tuple(orbit))
        occupations.append(occ)
    C = np.frombuffer(data, "<f8", norbits * norbits, pos).reshape(norbits, norbits).copy()
    pos += 8 * norbits * norbits
    blocks = []
    for _ in range(nchannels):
        J, parity, Tz, _, npairs = struct.unpack_from("<3iIQ", data, pos)
        pos += struct.calcsize("<3iIQ")
        pairs = [struct.unpack_from("<II", data, pos + 8 * k) for k in range(npairs)]
        pos += 8 * npairs
        matrix = np.frombuffer(data, "<f8", npairs * npairs, pos).reshape(npairs, npairs).copy()
        pos += 8 * npairs * npairs
        blocks.append((J, parity, Tz, pairs, matrix))
    return {"header": tuple(header), "digests": digests, "orbits": orbits,
            "occupations": np.array(occupations), "C": C, "blocks": blocks}


def write_arma(path, matrix):
    """Armadillo arma_binary double matrix, the format of imsrg++ <intfile>C.mat."""
    rows, cols = matrix.shape
    Path(path).write_bytes(f"ARMA_MAT_BIN_FN008\n{rows} {cols}\n".encode()
                           + np.asarray(matrix, "<f8").flatten(order="F").tobytes())


def run_imsrg(imsrg, bare, run, extra, env):
    run.mkdir(parents=True, exist_ok=True)
    command = [str(imsrg), f"2bme={bare}", *COMMON, f"flowfile={run}/flow.dat", f"intfile={run}/run",
               f"scratch={run}/", f"write_H_no2bpack={run}/H.no2bpack", *extra]
    with open(run / "imsrg.log", "w") as log:
        status = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env).returncode
    if status != 0:
        raise SystemExit(f"imsrg++ failed in {run}; see {run / 'imsrg.log'}")


def nmax0(bigstick, run):
    """BIGSTICK Nmax 0 (= valence space) of the exported H(s): E, 2+ excitation energy, <Rp2>."""
    work = run / "nmax0"
    energy, _ = ground_state(bigstick, run / "H.no2bpack", work, 0, (4, 4), 0)
    states = [{"E": energy}]
    expectation_values(SimpleNamespace(operator=[("Rp2", run / "H_Rp2.no2bpack")], hw=20, bigstick=bigstick),
                       work, "gs", states)
    ex2 = next(s["Ex"] for s in parse_res(work / "gs.res") if abs(s["J"] - 2) < 0.1)
    return {"E": energy, "Ex2+": ex2, "Rp2": states[0]["expectation"]["Rp2"]}


def snt_values(snt):
    """Zero-body term and all numbers of the .snt body (comments dropped)."""
    zero, numbers = None, []
    for line in Path(snt).read_text().splitlines():
        if "zero body" in line.lower():
            zero = float(line.split(":")[-1])
        if line.lstrip().startswith(("!", "#")):
            continue
        numbers += [float(x) for x in line.split()]
    return zero, np.array(numbers)


def kshell(snt, work):
    """Two lowest positive-parity KSHELL energies of the .snt plus its zero-body term (total energies)."""
    work.mkdir(parents=True, exist_ok=True)
    (work / "H.snt").write_text(Path(snt).read_text())
    answers = "N\nH.snt\nBe8\n\n+2\n0\nbeta_cm = 0\n\nN\n\n"
    ui = Path.home() / "kshell/bin/kshell_ui.py"
    with open(work / "kshell_ui.log", "w") as log:
        subprocess.run(["python3", str(ui)], input=answers, text=True, cwd=work, stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    with open(work / "kshell.log", "w") as log:
        subprocess.run(["bash", "Be8_H.sh"], cwd=work, stdout=log, stderr=subprocess.STDOUT, check=True)
    # final eigenvalues from the log (5 decimals; the summary needs `python` and keeps 3)
    energies = [float(line.split()[2]) for line in (work / "log_Be8_H_m0p.txt").read_text().splitlines()
                if "<H>:" in line and "<JJ>:" in line][-2:]
    zero, _ = snt_values(snt)
    return [zero + e for e in energies]


def gate_reference(template, occ_file, c_mat, labels, output):
    """lambda2 = 0 jref with the ensemble occupations of the SR run and C = C_HF (template: pair tables)."""
    occ = {}
    for line in Path(occ_file).read_text().splitlines():
        n, l, j2, tz2, value = line.split()
        occ[(int(n), int(l), int(j2), int(tz2))] = float(value)
    orbits = template["orbits"]
    labels4 = [o[1:] for o in orbits]
    write_jref(output, header=template["header"], orbits=orbits,
               occupations=np.array([occ.get(o, 0.0) for o in labels4]),
               transformation=basis_transformation_for(labels4, c_mat, labels),
               lambda_blocks=[(J, p, Tz, pairs, np.zeros_like(m)) for J, p, Tz, pairs, m in template["blocks"]],
               digests=template["digests"])


def next_reference(bigstick, run, jref_in, output):
    """Nmax 0 ground state of the exported (natural-basis) H(s), written as a jref in that basis."""
    work = run / "reference"
    ground_state(bigstick, run / "H.no2bpack", work, 0, (4, 4), 0)
    labels = read_no2bpack_orbits(run / "H.no2bpack")
    if jref_in is None:  # SR run: imsrg++ wrote HO -> HF
        c_mat = run / "runC.mat"
    else:  # MR run: H(s) is in the natural basis of jref_in; write its HO -> natural map in no2bpack order
        ref = read_jref(jref_in)
        index = {o[1:]: k for k, o in enumerate(ref["orbits"])}
        order = [index[label] for label in labels]
        c_mat = run / "C_total.mat"
        write_arma(c_mat, ref["C"][np.ix_(order, order)])
    convert(work / "gs.bas", work / "gs.wfn", run / "H.no2bpack", output, nrefmax=0, hw=20, e2max=4,
            lambda3=False, basis_transformation=c_mat, basis_orbits=run / "H.no2bpack")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--fixture-dir", type=Path, default=REPOSITORY / "prototype/mrimsrg/data")
    parser.add_argument("--iterations", type=int, default=3)
    parser.add_argument("--imsrg", type=Path, default=REPOSITORY / "build/src/imsrg++")
    parser.add_argument("--bigstick", type=Path, default=Path.home() / "BigstickPublick/bin/bigstick-openmp.x")
    args = parser.parse_args()
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="2")
    bigstick = str(args.bigstick)

    bare = work / "Be8_bare.jcoupled64"
    if not bare.exists():
        import pyIMSRG
        from UnitTestMRDriver import vacuum_operator
        modelspace = pyIMSRG.ModelSpace(2, "Be8", "Be8")
        modelspace.SetHbarOmega(20.0)
        fixture = load_reference(args.fixture_dir / "Be8_Nrefmax0_final")
        pyIMSRG.ReadWrite().Write_jcoupled64(str(bare), vacuum_operator(modelspace, fixture))

    results = {}

    def record(label, run, jref=None):
        r = nmax0(bigstick, run)
        r["kshell"] = kshell(run / "run.snt", run / "kshell")
        if jref is not None:
            ref = read_jref(jref)
            r["ref occ"] = {f"{o[1]}{o[2]}{o[3]}{o[4]:+d}": round(float(x), 6)
                            for o, x in zip(ref["orbits"], ref["occupations"]) if 1e-8 < x < 1 - 1e-8}
            r["|lambda2|"] = float(np.sqrt(sum(np.sum(b[4] ** 2) for b in ref["blocks"])))
            # <Phi|H(s)|Phi> (MR runs print it); at self-consistency it equals the Nmax 0 ground state
            r["<Phi|H(s)|Phi>"] = next(float(line.split("=")[1].split()[0]) for line in open(run / "imsrg.log")
                                       if line.startswith("IM-FCIQMC export check: H(s).ZeroBody"))
        results[label] = r
        print(f"{label:6s} E={r['E']:9.4f} ({r['E'] - EXACT['E']:+.4f})  2+={r['Ex2+']:.4f}  "
              f"Rp2={r['Rp2']:.4f} ({r['Rp2'] - EXACT['Rp2']:+.4f})  KSHELL-BIGSTICK={r['kshell'][0] - r['E']:+.1e}"
              + (f"  <Phi|H(s)|Phi>={r['<Phi|H(s)|Phi>']:9.4f}  |lambda2|={r['|lambda2|']:.3f}" if jref else ""),
              flush=True)

    print(f"exact  E={EXACT['E']:9.4f}  2+={EXACT['Ex2+']:.4f}  Rp2={EXACT['Rp2']:.4f}")
    sr = work / "sr"
    run_imsrg(args.imsrg, bare, sr, ["basis=HF", "write_omega=true"], env)
    record("sr", sr)

    jref = work / "ref1.jref"
    if not jref.exists():
        next_reference(bigstick, sr, None, jref)
    gate_jref = work / "gate.jref"
    if not gate_jref.exists():
        gate_reference(read_jref(jref), sr / "run_occ.dat", sr / "runC.mat", sr / "H.no2bpack", gate_jref)
    mr_extra = lambda j: ["basis=oscillator", f"mr_reference_file={j}", "mr_export_basis=natural"]
    gate = work / "gate"
    # lambda2 = 0 with fractional occupations violates the fixed-N contraction identity by n(1 - n) = 0.25;
    # that is the ensemble approximation itself, so only this wiring gate loosens the reference check.
    run_imsrg(args.imsrg, bare, gate, mr_extra(gate_jref) + ["mr_validation_tolerance=0.3"], env)
    record("gate", gate, gate_jref)
    zero_sr, snt_sr = snt_values(sr / "run.snt")
    zero_gate, snt_gate = snt_values(gate / "run.snt")
    results["gate .snt max|diff|"] = float(np.max(np.abs(snt_sr - snt_gate))) if snt_sr.shape == snt_gate.shape \
        else None
    results["gate zero-body diff"] = zero_gate - zero_sr
    print(f"gate vs sr: .snt max|diff| = {results['gate .snt max|diff|']}, zero body {zero_gate - zero_sr:+.2e}")

    for k in range(1, args.iterations + 1):
        run = work / f"mr{k}"
        run_imsrg(args.imsrg, bare, run, mr_extra(jref), env)
        record(f"mr{k}", run, jref)
        new = work / f"ref{k + 1}.jref"
        if not new.exists():
            next_reference(bigstick, run, jref, new)
        jref = new
    (work / "results.json").write_text(json.dumps({"exact": EXACT, **results}, indent=1))


if __name__ == "__main__":
    main()
