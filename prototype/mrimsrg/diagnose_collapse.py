#!/usr/bin/env python3
"""Collapse diagnostic: N-sector weights of IM-NCSM eigenstates for an Nref=0 and an Nref=2 reference.

He4, emax2, NNLOopt (fixtures in prototype/mrimsrg/data).  For each s the MR White-NCSM flow is run,
H(s) is exported in the HO basis (where the Nref space is defined) and BIGSTICK diagonalizes it in
the Nmax=4 space and inside the reference space (Nmax=2).  The ground state is split into sectors
N(D) = sum over occupied orbits of 2n+l minus the minimum.  With the operator-level Delta-e mask an
Nref>0 run can lose the N=0 <-> N=2 mixing inside its own reference space (collapse); the exact
full-space (Nmax=8) ground state, projected onto N<=2, gives the weight an ideal block decoupling keeps.

  diagnose_collapse.py --workdir DIR [--s 0 1 10 100 1000] [--imsrg build/src/imsrg++] [--bigstick ...]
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "test")]

from prototype.mrimsrg.bigstick_jref import BigstickBasis, configurations, read_bigstick_wfn  # noqa: E402
from prototype.mrimsrg.export_jref import export_reference  # noqa: E402
from prototype.mrimsrg.reference_io import load_reference  # noqa: E402

FIXTURES = {"Nref0": "He4_Nrefmax0_final", "Nref2": "He4_Nrefmax2"}


def weights(occupation, amplitudes, e_of_mstate):
    n = occupation.astype(np.int64) @ e_of_mstate
    n -= n.min()
    w = {int(N): float(np.sum(amplitudes[n == N] ** 2)) for N in np.unique(n)}
    total = sum(w.values())
    return {N: v / total for N, v in w.items()}


def ground_state(bigstick, h_path, workdir, nmax):
    """BIGSTICK ground state of a no2bpack in the He4 emax2 Nmax space: energy and N weights."""
    workdir.mkdir(parents=True, exist_ok=True)
    link = workdir / "h.bin"
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(Path(h_path).resolve())
    space = f"auto\n2\n2 2\n0\n+\ny\n{nmax}\n"
    solver = "ex\n2\n" if nmax <= 2 else "ld\n3 300\n"
    with open(workdir / "n.log", "w") as log:
        subprocess.run([bigstick], input=f"n\ngs\n{space}no2\nh.bin\n20 0\nend\n{solver}", text=True,
                       cwd=workdir, stdout=log, stderr=subprocess.STDOUT, check=True)
    if not (workdir / "gs.bas").exists():
        with open(workdir / "ba.log", "w") as log:
            subprocess.run([bigstick], input=f"ba\n{space}gs\n", text=True, cwd=workdir, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
    basis = BigstickBasis(workdir / "gs.bas")
    energy, _, amplitudes, _ = read_bigstick_wfn(workdir / "gs.wfn", basis, 1)
    mstates, occupation, _ = configurations(basis)
    e_of = np.array([2 * basis.orbits[o][0] + basis.orbits[o][1] for o, _ in mstates], dtype=np.int64)
    return energy, weights(occupation, amplitudes, e_of)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--s", type=float, nargs="+", default=[0, 0.3, 1, 3, 10, 30, 100, 1000])
    parser.add_argument("--imsrg", type=Path, default=REPOSITORY / "build/src/imsrg++")
    parser.add_argument("--bigstick", type=Path, default=Path.home() / "BigstickPublick/bin/bigstick-openmp.x")
    args = parser.parse_args()
    sys.path.insert(0, str(args.imsrg.parent))
    import pyIMSRG
    from UnitTestMRDriver import vacuum_operator

    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1")
    modelspace = pyIMSRG.ModelSpace(2, "He4", "He4")
    modelspace.SetHbarOmega(20.0)
    bare = work / "He4_bare.jcoupled64"
    data0 = load_reference(REPOSITORY / "prototype/mrimsrg/data" / FIXTURES["Nref0"])
    pyIMSRG.ReadWrite().Write_jcoupled64(str(bare), vacuum_operator(modelspace, data0))
    summary = {}
    for label, fixture in FIXTURES.items():
        data = load_reference(REPOSITORY / "prototype/mrimsrg/data" / fixture)
        orbits = np.asarray(data.orbits)
        ref_w = weights(np.asarray(data.determinants), np.asarray(data.coefficients),
                        (2 * orbits[:, 1] + orbits[:, 2]).astype(np.int64))
        jref = work / f"{label}.jref"
        if not jref.exists():
            export_reference(REPOSITORY / "prototype/mrimsrg/data" / fixture, jref)
        rows = []
        for s in args.s:
            tag = f"{label}_s{s:g}"
            omega = work / f"omega_{tag}"
            omega.mkdir(exist_ok=True)
            h = work / f"H_{tag}.no2bpack"
            run = subprocess.run(
                [str(args.imsrg), f"2bme={bare}", "fmt2=jcoupled64", "3bme=none", "reference=He4",
                 "valence_space=He4", "A=4", "hw=20", "emax=2", "basis=oscillator", "method=magnus_adaptive",
                 "nsteps=1", "core_generator=white-ncsm", "denominator_partitioning=Epstein_Nesbet", "BetaCM=0",
                 "nucleon_mass_correction=false", f"mr_reference_file={jref}", f"smax={s:g}",
                 "eta_criterion=1e-12", f"scratch={omega}", f"flowfile={work / ('flow_' + tag + '.dat')}",
                 f"intfile={work / tag}", f"write_H_no2bpack={h}"],
                capture_output=True, text=True, env=env)
            if run.returncode != 0:
                raise SystemExit(f"imsrg++ failed for {tag}:\n{run.stdout[-2000:]}")
            flow = [line.split() for line in open(work / f"flow_{tag}.dat") if line.strip()]
            e4, w4 = ground_state(str(args.bigstick), h, work / f"N4_{tag}", 4)
            e2, w2 = ground_state(str(args.bigstick), h, work / f"N2_{tag}", 2)
            rows.append({"s": s, "E_s": float(flow[-1][2]), "Nmax4": {"E": e4, "weights": w4},
                         "Nmax2": {"E": e2, "weights": w2}})
            print(f"{label} s={s:<6g} E(s)={rows[-1]['E_s']:9.4f} | Nmax4 E={e4:9.4f} N0 {w4.get(0, 0):.4f} "
                  f"N2 {w4.get(2, 0):.4f} N4 {w4.get(4, 0):.4f} | Nmax2 E={e2:9.4f} N2 {w2.get(2, 0):.4f}", flush=True)
        summary[label] = {"reference_weights": ref_w, "rows": rows}
    e8, w8 = ground_state(str(args.bigstick), work / "H_Nref0_s0.no2bpack", work / "exact_Nmax8", 8)
    p = w8[0] + w8[2]
    summary["exact"] = {"E": e8, "weights": w8, "projected_N2": w8[2] / p}
    print(f"exact (Nmax 8) E={e8:.5f}; N<=2 projection keeps N2 {w8[2] / p:.4f}; "
          f"Nref2 reference N2 {summary['Nref2']['reference_weights'].get(2, 0):.4f}")
    (work / "collapse.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
