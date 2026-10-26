#!/usr/bin/env python3
"""IM-NCSM post-diagonalization of an exported Hamiltonian with BIGSTICK.

The MR driver writes H(s) as a vacuum 0/1/2-body no2bpack (minipack layout);
BIGSTICK's `no2` reader takes it directly (zero body included, one body
folded into the two-body part, any single-particle basis).  This script runs
one BIGSTICK Lanczos per Nmax, parses the spectrum and writes JSON.  The
oscillator centre-of-mass strength is fixed to 0: imsrg++ writes zero H_cm
components, so no Lawson term can be added.

Usage:
  postdiag_bigstick.py --no2bpack H.no2bpack --Z 8 --N 8 --emax 8 --hw 20 \
      --nmax 0 2 4 [--states 5] [--max-iter 300] [--parity +] [--jz2 0] \
      [--bigstick bigstick-openmp.x] --workdir DIR --output result.json
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess
import time

EXACT_MAX = 2000   # M-scheme dimension up to which BIGSTICK diagonalizes fully


def parse_res(path: Path):
    """States (E, Ex, J, T, parity) from a BIGSTICK .res file."""
    states = []
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.split()[:3] == ["State", "E", "Ex"]:
            for row in lines[i + 1:]:
                fields = row.split()
                if len(fields) < 6 or not fields[0].isdigit():
                    break
                states.append({"E": float(fields[1]), "Ex": float(fields[2]),
                               "J": float(fields[3]), "T": float(fields[4]), "parity": int(fields[5])})
            break
    return states


def run_nmax(args, nmax: int):
    work = args.workdir / f"Nmax{nmax}"
    work.mkdir(parents=True, exist_ok=True)
    link = work / "hamiltonian_s.bin"   # BIGSTICK takes names starting with "no2" for the format code
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(args.no2bpack.resolve())
    name = f"Z{args.Z}N{args.N}_Nmax{nmax}"
    start = time.time()
    # Lanczos ("ld"); small spaces break it (a closed-shell Nmax=0 determinant
    # fails outright, the 81-dimensional O20 Nmax=0 space gave a non-J-pure
    # state 26 MeV below the exact ground state), so spaces up to EXACT_MAX
    # use full diagonalization ("ex").
    for method, solver in (("lanczos", ["ld", f"{args.states} {args.max_iter}"]),
                           ("exact", ["ex", str(args.states)])):
        deck = "\n".join(["n", name, "auto", str(args.emax), f"{args.Z} {args.N}", str(args.jz2),
                          args.parity, "y", str(nmax), "no2", link.name, f"{args.hw} 0", "end",
                          *solver, ""])
        (work / "input.txt").write_text(deck)
        res = work / f"{name}.res"
        if res.exists():
            res.unlink()
        with (work / "input.txt").open() as stdin, (work / "bigstick.log").open("w") as log:
            status = subprocess.run([str(args.bigstick)], stdin=stdin, stdout=log,
                                    stderr=subprocess.STDOUT, cwd=work, check=False).returncode
        states = parse_res(res) if res.is_file() else []
        match = re.search(r"Total basis\s*=\s*(\d+)", (work / "bigstick.log").read_text(errors="replace"))
        dimension = int(match.group(1)) if match else None
        if status == 0 and states and (method == "exact" or (dimension or 0) > EXACT_MAX):
            break
    else:
        raise RuntimeError(f"BIGSTICK failed for Nmax={nmax}; see {work / 'bigstick.log'}")
    return {"Nmax": nmax, "dimension": dimension, "method": method,
            "seconds": round(time.time() - start, 1), "states": states}


def extrapolate(nmax_values, energies):
    """E(N) = E_inf + a exp(-b N) through the last three points (None if not monotone)."""
    if len(energies) < 3:
        return None
    (n1, n2, n3), (e1, e2, e3) = nmax_values[-3:], energies[-3:]
    if not (n2 - n1 == n3 - n2):
        return None
    d1, d2 = e2 - e1, e3 - e2
    if d1 == 0.0 or d2 / d1 <= 0.0 or d2 / d1 >= 1.0:
        return None
    ratio = d2 / d1
    return e3 + d2 * ratio / (1.0 - ratio)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--no2bpack", required=True, type=Path)
    parser.add_argument("--Z", required=True, type=int)
    parser.add_argument("--N", required=True, type=int)
    parser.add_argument("--emax", required=True, type=int)
    parser.add_argument("--hw", required=True, type=float)
    parser.add_argument("--nmax", required=True, type=int, nargs="+")
    parser.add_argument("--states", type=int, default=5)
    parser.add_argument("--max-iter", type=int, default=300)
    parser.add_argument("--parity", default="+")
    parser.add_argument("--jz2", type=int, default=0)
    parser.add_argument("--bigstick", type=Path, default=Path(os.environ.get(
        "BIGSTICK", "/home/mengziyan/BigstickPublick/bin/bigstick-openmp.x")))
    parser.add_argument("--workdir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    results = [run_nmax(args, nmax) for nmax in sorted(args.nmax)]
    ground = [r["states"][0]["E"] for r in results]
    summary = {
        "no2bpack": str(args.no2bpack.resolve()), "Z": args.Z, "N": args.N, "emax": args.emax, "hw": args.hw,
        "results": results,
        "ground_state": dict(zip([r["Nmax"] for r in results], ground)),
        "extrapolated_ground_state": extrapolate([r["Nmax"] for r in results], ground),
    }
    args.output.write_text(json.dumps(summary, indent=1))
    for r in results:
        s = r["states"]
        print(f"Nmax={r['Nmax']:2d} dim={r['dimension']} E0={s[0]['E']:.4f} "
              f"lowest: {[(round(x['E'], 3), round(x['J'], 2)) for x in s[:4]]}")
    if summary["extrapolated_ground_state"] is not None:
        print(f"exponential extrapolation (last three Nmax): {summary['extrapolated_ground_state']:.3f}")


if __name__ == "__main__":
    main()
