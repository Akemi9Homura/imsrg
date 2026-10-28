#!/usr/bin/env python3
"""Workflow check of the valence three-body output (valence_3b_file) in the p shell above a He4 core.

N2LO_opt bare hw20, emax2 (me2j file truncated), NN only, BetaCM = 0, basis=HF, method=magnus.
  Be8  perturbative W with one gathered Omega (default) | per Magnus segment | optional IMSRG(3) flow
  C12  perturbative W, gathered Omega: the ensemble reference is the Slater determinant p3/2^4 x p3/2^4, on which
       <Phi0| :W:_ref |Phi0> must vanish exactly (any error in the lower-body terms or the file convention shows).
Checks: the standard .snt reproduces the earlier KSHELL/BIGSTICK numbers; the lower-body terms written by
imsrg++ (difference of the two .snt) equal -sum n W, +1/2 sum n n W, -1/6 sum n n n W built here in the
m-scheme from the three-body file; the C12 determinant check.  Physics (E, 2+, <Rp2>) is printed next to the
exact Be8 values (BIGSTICK Nmax 12 = full emax2; mr-imsrg-dev doc/records/工作记录.md 2026-10-27), but at
emax2 the three-body is expected to be tiny (Johnson, Flow & Tell 2026: Be8 emax2 net 0.08 MeV).

  check_be8.py --workdir DIR [--imsrg build/src/imsrg++]
"""

import argparse
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pspace  # noqa: E402

REPOSITORY = Path(__file__).resolve().parents[2]
NN = "/home/mengziyan/Forces/N2LO_opt/TwBME_N2LO_opt_bare_hw20_emax12_e2max24.me2j.gz"
EXACT = {"E": -26.325, "Ex2+": 3.340, "Rp2": 4.209}
P32 = {(0, 1, 3, -1), (0, 1, 3, 1)}


def run(imsrg, work, tag, nucleus, A, extra, scratch=True):
    """One imsrg++ run; with scratch the Omega segments are flushed to disk and W reads them back."""
    d = work / tag
    d.mkdir(parents=True, exist_ok=True)
    if scratch:
        (d / "scratch").mkdir(exist_ok=True)
        extra = extra + [f"scratch={d}/scratch"]
    cmd = [str(imsrg), f"2bme={NN}", "fmt2=me2j", "file2e1max=12", "file2e2max=24", "3bme=none", "hw=20",
           "emax=2", "e3max=6", "BetaCM=0", "nucleon_mass_correction=false", "basis=HF", "method=magnus",
           "valence_space=p-shell", f"reference={nucleus}", f"A={A}", "Operators=Rp2",
           f"flowfile={d}/flow.dat", f"intfile={d}/run", f"valence_3b_file={d}/v3b", *extra]
    if not (d / "v3b_3b.int").exists():
        with open(d / "imsrg.log", "w") as log:
            if subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                              env=dict(os.environ, OMP_NUM_THREADS="2")).returncode != 0:
                raise SystemExit(f"imsrg++ failed: {d / 'imsrg.log'}")
    return d


def lowest(states):
    gs = states[0]
    ex2 = next((s["E"] - gs["E"] for s in states[1:] if abs(s["J"] - 2) < 0.1), None)
    return gs, ex2


def base_snt(d):
    """The two-body Hamiltonian that v3b.snt adds the lower-body terms of W to."""
    return d / "v3b_2b.snt" if (d / "v3b_2b.snt").exists() else d / "run.snt"


def check_lower_body(d, occ_by_orbit):
    """imsrg++ lower-body terms (v3b.snt minus its base) against the m-scheme formula from the three-body file."""
    std, v3b = pspace.read_snt(base_snt(d)), pspace.read_snt(d / "v3b.snt")
    space = pspace.Space(std["orbits"])
    one_s, two_s, _ = pspace.hamiltonian(std, space)
    one_v, two_v, _ = pspace.hamiltonian(v3b, space)
    three = pspace.mscheme_from_coupled(space, pspace.three_body_blocks(
        pspace.read_3b(d / "v3b_3b.int", std["orbits"])), 3)
    occ = [occ_by_orbit.get(std["orbits"][o], 0.0) for o, m in space.mstates]
    zero, one, two = pspace.lower_body_from_three(three, occ)
    d1 = max(abs(one_v.get(k, 0) - one_s.get(k, 0) - one.get(k, 0)) for k in set(one_v) | set(one_s) | set(one))
    d2 = max(abs(two_v.get(k, 0) - two_s.get(k, 0) - two.get(k, 0)) for k in set(two_v) | set(two_s) | set(two))
    d0 = v3b["zero"] - std["zero"] - zero
    scale = max([abs(v) for v in three.values()] + [1e-300])
    return {"zero": d0, "one": d1, "two": d2, "max|W|": scale, "E0p": v3b["zero"] - std["zero"]}


def decomposition(d, Zv, Nv):
    """Expectation values in the ground state of the base Hamiltonian: three-body, lower-body one+two-body, zero-body."""
    std, v3b = pspace.read_snt(base_snt(d)), pspace.read_snt(d / "v3b.snt")
    states, (space, basis, X, H) = pspace.spectrum(base_snt(d), Zv, Nv)
    one_s, two_s, _ = pspace.hamiltonian(std, space)
    one_v, two_v, three = pspace.hamiltonian(v3b, space, pspace.read_3b(d / "v3b_3b.int", std["orbits"]))
    diff1 = {k: one_v.get(k, 0) - one_s.get(k, 0) for k in set(one_v) | set(one_s)}
    diff2 = {k: two_v.get(k, 0) - two_s.get(k, 0) for k in set(two_v) | set(two_s)}
    x = X[:, 0]
    m3 = pspace.many_body_matrix(basis, None, None, three)
    m12 = pspace.many_body_matrix(basis, diff1, diff2, None)
    return {"3B": float(x @ m3 @ x), "1+2B": float(x @ m12 @ x), "0B": v3b["zero"] - std["zero"]}


def compare_w(a, b):
    """|W_a - W_b| / |W_b| over the J-coupled elements of two three-body files."""
    va = {e[:9]: e[9] for e in pspace.read_3b(a)}
    vb = {e[:9]: e[9] for e in pspace.read_3b(b)}
    keys = set(va) | set(vb)
    num = sum((va.get(k, 0.0) - vb.get(k, 0.0)) ** 2 for k in keys)
    den = sum(v ** 2 for v in vb.values())
    return (num / den) ** 0.5 if den > 0 else float("nan"), den ** 0.5


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--imsrg", type=Path, default=REPOSITORY / "build/src/imsrg++")
    parser.add_argument("--with-imsrg3", action="store_true",
                        help="also run the IMSRG(3) flow benchmark (slow: hours on two threads)")
    args = parser.parse_args()
    work = args.workdir.resolve()
    runs = {
        "Be8_gathered": ("Be8", 8, ["valence_3b_omega=gathered"]),
        "Be8_segments": ("Be8", 8, ["valence_3b_omega=segments"]),
        "C12_gathered": ("C12", 12, ["valence_3b_omega=gathered"]),
        "Be8_at_end": ("Be8", 8, ["valence_3b_source=at_end"]),
        "C12_at_end": ("C12", 12, ["valence_3b_source=at_end"]),
    }
    if args.with_imsrg3:
        runs["Be8_imsrg3"] = ("Be8", 8, ["IMSRG3=true", "valence_3b_source=flow"])
    dirs = {tag: run(args.imsrg, work, tag, *spec) for tag, spec in runs.items()}
    memory = run(args.imsrg, work, "Be8_gathered_memory", "Be8", 8, ["valence_3b_omega=gathered"], scratch=False)
    same = all((memory / f).read_text().split("\n")[1:] == (dirs["Be8_gathered"] / f).read_text().split("\n")[1:]
               for f in ("v3b.snt", "v3b_3b.int"))
    print(f"Omega from scratch vs in memory: v3b.snt and v3b_3b.int identical apart from the header line: {same}")

    print(f"exact Be8: E {EXACT['E']:.4f}  2+ {EXACT['Ex2+']:.4f}  Rp2 {EXACT['Rp2']:.4f}")
    for tag in [t for t in ("Be8_gathered", "Be8_segments", "Be8_at_end", "Be8_imsrg3") if t in dirs]:
        d = dirs[tag]
        ops = [("Rp2", d / "runRp2.snt")]
        std, _ = pspace.spectrum(d / "run.snt", 2, 2, operator_paths=ops)
        w3, _ = pspace.spectrum(d / "v3b.snt", 2, 2, three_path=d / "v3b_3b.int", operator_paths=ops)
        (g0, x0), (g1, x1) = lowest(std), lowest(w3)
        print(f"{tag:13s} without W: E {g0['E']:.4f}  2+ {x0:.4f}  Rp2 {g0['Rp2']:.4f} | with W: E {g1['E']:.4f}  "
              f"2+ {x1:.4f}  Rp2 {g1['Rp2']:.4f} | shift E {g1['E'] - g0['E']:+.4f}  2+ {x1 - x0:+.4f}")
        if base_snt(d) != d / "run.snt":
            b2, _ = pspace.spectrum(base_snt(d), 2, 2, operator_paths=ops)
            gb, xb = lowest(b2)
            print(f"{'':13s} two-body feedback only (v3b_2b.snt): E {gb['E']:.4f}  2+ {xb:.4f}  Rp2 {gb['Rp2']:.4f}"
                  f" | shift E {gb['E'] - g0['E']:+.4f}  2+ {xb - x0:+.4f}")
        dec = decomposition(d, 2, 2)
        print(f"{'':13s} in the ground state without W: 3B {dec['3B']:+.5f}  1+2B {dec['1+2B']:+.5f}  "
              f"0B {dec['0B']:+.5f}  sum {sum(dec.values()):+.5f}")
    for tag, occ in (("Be8_gathered", {o: 0.5 for o in P32}), ("C12_gathered", {o: 1.0 for o in P32}),
                     ("Be8_at_end", {o: 0.5 for o in P32}), ("C12_at_end", {o: 1.0 for o in P32})):
        c = check_lower_body(dirs[tag], occ)
        print(f"{tag:13s} lower-body terms, imsrg++ minus m-scheme formula: zero {c['zero']:+.1e}  one {c['one']:.1e}  "
              f"two {c['two']:.1e}  (max |W| {c['max|W|']:.2e}, E0' {c['E0p']:+.6f})")
    # W relative to the at_end (and, if run, the IMSRG(3) flow) three-body
    for ref in [r for r in ("Be8_at_end", "Be8_imsrg3") if r in dirs]:
        for tag in [t for t in ("Be8_gathered", "Be8_segments", "Be8_at_end", "Be8_imsrg3") if t in dirs and t != ref]:
            rel, norm = compare_w(dirs[tag] / "v3b_3b.int", dirs[ref] / "v3b_3b.int")
            print(f"|W({tag}) - W({ref})| / |W({ref})| = {rel:.3f}   (|W({ref})| = {norm:.4f})")
    # C12: <Phi0|:W:|Phi0> = 0 for the reference determinant p3/2^4 x p3/2^4
    for tag in ("C12_gathered", "C12_at_end"):
        c12_determinant_check(dirs[tag])


def c12_determinant_check(d):
    std, v3b = pspace.read_snt(base_snt(d)), pspace.read_snt(d / "v3b.snt")
    space = pspace.Space(std["orbits"])
    phi0 = sum(1 << k for k, (o, m) in enumerate(space.mstates) if std["orbits"][o] in P32)
    basis = pspace.Basis(space, 4, 4, 0)
    e = np.zeros(len(basis.dets))
    e[basis.lookup[phi0]] = 1.0
    hs = pspace.many_body_matrix(basis, *pspace.hamiltonian(std, space)[:2], None, std["zero"])
    hv = pspace.many_body_matrix(basis, *pspace.hamiltonian(v3b, space, pspace.read_3b(d / "v3b_3b.int",
                                                                                         std["orbits"])), v3b["zero"])
    print(f"{d.name}: <Phi0|base + :W:|Phi0> - <Phi0|base|Phi0> = {e @ hv @ e - e @ hs @ e:+.2e} MeV (must vanish)")


if __name__ == "__main__":
    main()
