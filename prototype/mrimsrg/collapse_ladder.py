#!/usr/bin/env python3
"""Collapse ladder: IM-NCSM with Nref = 0, 2, ... up to the full space (emax2, NNLOopt; He4 or Be8).

Each reference is the BIGSTICK Nmax = Nref ground state of the bare H (HO basis, jref with lambda3).  It is
flowed with the MR White-NCSM generator, H(s) is exported in the HO basis and diagonalized in P
(Nmax = Nref) and in the full space (He4 Nmax 8, Be8 Nmax 12); the ground state is split into N sectors,
its intrinsic point-proton <Rp2> is read with BIGSTICK `x` from the evolved Rp2 (it connects Delta N = 2
directly), and the lowest 2+ excitation energy is taken from the same diagonalization.  An ideal P/Q block
decoupling would make the P results approach the exact ones as Nref grows and leave Nref = full untouched
(Q empty, generator zero).  The Delta-e mask instead rotates the in-space N mixing away for every Nref
(collapse): He4 Nref >= 4 end on nearly the same H(s) and Nref = 8 moves away from the exact energy.

The bare H, the exact values and each reference are computed once and reused, so runs with different
--nref or --s can share a workdir (the bare files can also be copied in, then pyIMSRG is not needed).
Each (Nref, s) point writes Nref<n>_s<s>/result.json.

  collapse_ladder.py --workdir DIR [--nucleus He4|Be8] [--nref 0 2 4 6 8] [--s 0 1 10 100]
                     [--imsrg build/src/imsrg++] [--bigstick ...]
"""

import argparse
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "test")]

from prototype.mrimsrg.bigstick_jref import convert  # noqa: E402
from prototype.mrimsrg.diagnose_collapse import ground_state, run_flow  # noqa: E402
from prototype.mrimsrg.export_jref import export_reference  # noqa: E402
from prototype.mrimsrg.postdiag_bigstick import expectation_values, parse_res  # noqa: E402
from prototype.mrimsrg.reference_io import load_reference  # noqa: E402

# fixture with the bare interaction, (Z, N), full-space Nmax in emax2, largest Nmax diagonalized exactly
NUCLEI = {"He4": ("He4_Nrefmax0_final", (2, 2), 8, 2), "Be8": ("Be8_Nrefmax0_final", (4, 4), 12, 0)}


def rp2(bigstick, workdir, operator, energy):
    """<Rp2> of the ground state that ground_state() left in workdir (BIGSTICK `x`)."""
    states = [{"E": energy}]
    expectation_values(SimpleNamespace(operator=[("Rp2", operator)], hw=20, bigstick=bigstick), workdir, "gs", states)
    return states[0]["expectation"]["Rp2"]


def excitation_2plus(workdir):
    """Excitation energy of the lowest J = 2 state in workdir/gs.res (None if not among the computed states)."""
    states = parse_res(workdir / "gs.res")
    return next((state["Ex"] for state in states if abs(state["J"] - 2) < 0.1), None)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--nucleus", choices=sorted(NUCLEI), default="He4")
    parser.add_argument("--nref", type=int, nargs="+", default=[0, 2, 4, 6, 8])
    parser.add_argument("--s", type=float, nargs="+", default=[0, 1, 10, 100])
    parser.add_argument("--imsrg", type=Path, default=REPOSITORY / "build/src/imsrg++")
    parser.add_argument("--bigstick", type=Path, default=Path.home() / "BigstickPublick/bin/bigstick-openmp.x")
    args = parser.parse_args()
    fixture_name, nucleons, full_nmax, exact_nmax = NUCLEI[args.nucleus]
    A = sum(nucleons)

    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1")
    bigstick = str(args.bigstick)
    fixture = REPOSITORY / "prototype/mrimsrg/data" / fixture_name
    bare = work / f"{args.nucleus}_bare.jcoupled64"
    if not bare.exists():
        sys.path.insert(0, str(args.imsrg.parent))
        import pyIMSRG
        from UnitTestMRDriver import vacuum_operator
        modelspace = pyIMSRG.ModelSpace(2, args.nucleus, args.nucleus)
        modelspace.SetHbarOmega(20.0)
        pyIMSRG.ReadWrite().Write_jcoupled64(str(bare), vacuum_operator(modelspace, load_reference(fixture)))
    hamiltonian = work / "H_bare.no2bpack"
    if not hamiltonian.exists():
        export_reference(fixture, work / "fixture_Nref0.jref")
        run_flow(args.imsrg, bare, work / "fixture_Nref0.jref", 0, hamiltonian, work / "bare",
                 work / "bare/flow.dat", work / "bare/", env, operators=("Rp2",), nucleus=args.nucleus, A=A)

    def diagonalize(h, workdir, nmax, rp2_operator):
        energy, w = ground_state(bigstick, h, workdir, nmax, nucleons, exact_nmax)
        return {"E": energy, "Rp2": rp2(bigstick, workdir, rp2_operator, energy),
                "Ex2+": excitation_2plus(workdir), "weights": w}

    cache = work / "exact.json"
    if not cache.exists():
        cache.write_text(json.dumps(diagonalize(hamiltonian, work / "exact", full_nmax, work / "H_bare_Rp2.no2bpack")))
    exact = json.loads(cache.read_text())
    print(f"exact (Nmax {full_nmax}) E = {exact['E']:.5f}, <Rp2> = {exact['Rp2']:.5f} fm^2, "
          f"Ex(2+) = {exact['Ex2+']}", flush=True)
    for nref in args.nref:
        reference = work / f"ref{nref}"
        jref = work / f"Nref{nref}.jref"
        if not jref.exists():
            _, w_ref = ground_state(bigstick, hamiltonian, reference, nref, nucleons, exact_nmax)
            (reference / "weights.json").write_text(json.dumps(w_ref))
            convert(reference / "gs.bas", reference / "gs.wfn", hamiltonian, jref, nrefmax=nref, hw=20, e2max=4)
        w_ref = json.loads((reference / "weights.json").read_text())
        print(f"Nref{nref} reference " + " ".join(f"N{k}:{v:.4f}" for k, v in sorted(w_ref.items(), key=lambda x: int(x[0]))),
              flush=True)
        for s in args.s:
            run = work / f"Nref{nref}_s{s:g}"
            h = run / "H.no2bpack"
            e_s = run_flow(args.imsrg, bare, jref, s, h, run, run / "flow.dat", f"{run}/", env,
                           operators=("Rp2",), nucleus=args.nucleus, A=A)
            p = diagonalize(h, run / "P", nref, run / "H_Rp2.no2bpack")
            full = diagonalize(h, run / "full", full_nmax, run / "H_Rp2.no2bpack")
            (run / "result.json").write_text(json.dumps({"nucleus": args.nucleus, "Nref": nref, "s": s, "E_s": e_s,
                                                         "reference_weights": w_ref, "exact": exact,
                                                         "P": p, "full": full}, indent=1))
            ex = lambda d: f"{d['Ex2+']:.3f}" if d["Ex2+"] is not None else "-"
            print(f"  s={s:<6g} E(s)={e_s:9.4f} | P: E={p['E']:9.4f} ({p['E'] - exact['E']:+.4f}) Rp2 {p['Rp2']:.4f} "
                  f"({p['Rp2'] - exact['Rp2']:+.4f}) 2+ {ex(p)} N0 {p['weights'].get(0, 0):.4f} | full: "
                  f"E={full['E']:9.4f} ({full['E'] - exact['E']:+.4f}) Rp2 {full['Rp2']:.4f} 2+ {ex(full)} "
                  f"N0 {full['weights'].get(0, 0):.4f}", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
