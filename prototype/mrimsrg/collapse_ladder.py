#!/usr/bin/env python3
"""Collapse ladder: IM-NCSM with Nref = 0, 2, ... up to the full space (He4 emax2, NNLOopt, full space = Nmax 8).

Each reference is the BIGSTICK Nmax = Nref ground state of the bare H (HO basis, jref with lambda3).  It is
flowed with the MR White-NCSM generator, H(s) is exported in the HO basis and diagonalized in P
(Nmax = Nref) and in the full space; the ground state is split into N sectors.  An ideal P/Q block
decoupling would make the P energy approach the exact one as Nref grows and leave Nref = 8 untouched
(Q empty, generator zero).  The Delta-e mask instead rotates the in-space N mixing away for every Nref
(collapse): Nref >= 4 end on nearly the same H(s) and Nref = 8 moves away from the exact energy.

  collapse_ladder.py --workdir DIR [--nref 0 2 4 6 8] [--s 0 1 10 100] [--imsrg build/src/imsrg++] [--bigstick ...]
"""

import argparse
import json
import os
from pathlib import Path
import sys

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "test")]

from prototype.mrimsrg.bigstick_jref import convert  # noqa: E402
from prototype.mrimsrg.diagnose_collapse import FIXTURES, ground_state, run_flow  # noqa: E402
from prototype.mrimsrg.export_jref import export_reference  # noqa: E402
from prototype.mrimsrg.reference_io import load_reference  # noqa: E402

FULL_NMAX = 8


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--nref", type=int, nargs="+", default=[0, 2, 4, 6, 8])
    parser.add_argument("--s", type=float, nargs="+", default=[0, 1, 10, 100])
    parser.add_argument("--imsrg", type=Path, default=REPOSITORY / "build/src/imsrg++")
    parser.add_argument("--bigstick", type=Path, default=Path.home() / "BigstickPublick/bin/bigstick-openmp.x")
    args = parser.parse_args()
    sys.path.insert(0, str(args.imsrg.parent))
    import pyIMSRG
    from UnitTestMRDriver import vacuum_operator

    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1")
    bigstick = str(args.bigstick)
    fixture = REPOSITORY / "prototype/mrimsrg/data" / FIXTURES["Nref0"]
    bare = work / "He4_bare.jcoupled64"
    if not bare.exists():
        modelspace = pyIMSRG.ModelSpace(2, "He4", "He4")
        modelspace.SetHbarOmega(20.0)
        pyIMSRG.ReadWrite().Write_jcoupled64(str(bare), vacuum_operator(modelspace, load_reference(fixture)))
    hamiltonian = work / "H_bare.no2bpack"
    if not hamiltonian.exists():
        export_reference(fixture, work / "fixture_Nref0.jref")
        run_flow(args.imsrg, bare, work / "fixture_Nref0.jref", 0, hamiltonian, work / "bare",
                 work / "bare/flow.dat", work / "bare/", env)
    exact, _ = ground_state(bigstick, hamiltonian, work / "exact", FULL_NMAX)
    summary = {"exact": exact}
    for nref in args.nref:
        reference = work / f"ref{nref}"
        _, w_ref = ground_state(bigstick, hamiltonian, reference, nref)
        jref = work / f"Nref{nref}.jref"
        if not jref.exists():
            convert(reference / "gs.bas", reference / "gs.wfn", hamiltonian, jref, nrefmax=nref, hw=20, e2max=4)
        print(f"Nref{nref} reference " + " ".join(f"N{k}:{v:.4f}" for k, v in sorted(w_ref.items())), flush=True)
        rows = []
        for s in args.s:
            run = work / f"Nref{nref}_s{s:g}"
            h = run / "H.no2bpack"
            e_s = run_flow(args.imsrg, bare, jref, s, h, run, run / "flow.dat", f"{run}/", env)
            e_p, w_p = ground_state(bigstick, h, run / "P", nref)
            e_full, w_full = ground_state(bigstick, h, run / "full", FULL_NMAX)
            rows.append({"s": s, "E_s": e_s, "P": {"E": e_p, "weights": w_p},
                         "full": {"E": e_full, "weights": w_full}})
            print(f"  s={s:<6g} E(s)={e_s:9.4f} | P: E={e_p:9.4f} ({e_p - exact:+.4f}) N0 {w_p.get(0, 0):.4f} "
                  f"N2 {w_p.get(2, 0):.4f} | full: E={e_full:9.4f} ({e_full - exact:+.4f}) N0 {w_full.get(0, 0):.4f}",
                  flush=True)
        summary[nref] = {"reference_weights": w_ref, "rows": rows}
        (work / "ladder.json").write_text(json.dumps(summary, indent=1))
    print(f"exact (Nmax {FULL_NMAX}) E = {exact:.5f}")


if __name__ == "__main__":
    main()
