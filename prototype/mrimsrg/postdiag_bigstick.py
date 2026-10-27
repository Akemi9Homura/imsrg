#!/usr/bin/env python3
"""IM-NCSM post-diagonalization of an exported Hamiltonian with BIGSTICK.

The MR driver writes H(s) as a vacuum 0/1/2-body no2bpack (minipack layout);
BIGSTICK's `no2` reader takes it directly (zero body included, one body
folded into the two-body part, any single-particle basis).  This script runs
one BIGSTICK Lanczos per Nmax, parses the spectrum and writes JSON.  The
oscillator centre-of-mass strength is fixed to 0: imsrg++ writes zero H_cm
components (a Lawson term goes in with no2bpack_combine.py instead).

Scalar operators exported with H (e.g. H_HCM.no2bpack, H_Rp2.no2bpack from
`Operators=HCM,Rp2`) are evaluated in every eigenstate with BIGSTICK's `x`
option, which applies the operator to the saved wavefunctions: the full
zero-, one- and two-body expectation value <psi|O(s)|psi>.

Usage:
  postdiag_bigstick.py --no2bpack H.no2bpack --Z 8 --N 8 --emax 8 --hw 20 \
      --nmax 0 2 4 [--states 5] [--max-iter 300] [--parity +] [--jz2 0] \
      [--operator HCM=H_HCM.no2bpack ...] \
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


def link(path: Path, target: Path):
    if path.is_symlink() or path.exists():
        path.unlink()
    path.symlink_to(target.resolve())


def expectation_values(args, work: Path, name: str, states):
    """<psi_i|O|psi_i> for the states saved in <name>.wfn, one BIGSTICK `x` run per operator."""
    values = {}
    for op_name, op_file in args.operator:
        link(work / f"op_{op_name}.bin", op_file)
        deck = "\n".join(["x", name, f"x_{op_name}", "no2", f"op_{op_name}.bin", f"{args.hw} 0", "end", ""])
        with (work / f"x_{op_name}.log").open("w") as log:
            status = subprocess.run([str(args.bigstick)], input=deck, text=True, stdout=log,
                                    stderr=subprocess.STDOUT, cwd=work, check=False).returncode
        res = work / f"x_{op_name}.res"
        if status != 0 or not res.is_file():
            raise RuntimeError(f"BIGSTICK x failed for {op_name}; see {work / f'x_{op_name}.log'}")
        lines = res.read_text().splitlines()
        start = next(i for i, line in enumerate(lines) if "<H >" in line) + 1
        rows = []
        for line in lines[start:]:
            fields = line.split()
            if len(fields) != 6 or not fields[0].isdigit():
                break
            rows.append(fields)
        for fields in rows:
            index, energy, norm, value = int(fields[0]) - 1, float(fields[1]), float(fields[5]), float(fields[4])
            if index >= len(states):
                continue
            # same wavefunction file as the spectrum: energies agree to the printed precision
            if abs(energy - states[index]["E"]) > 2e-3 or abs(norm - 1.0) > 1e-4:
                raise RuntimeError(f"BIGSTICK x for {op_name}: state {index + 1} has E={energy}, norm={norm}; "
                                   f"spectrum E={states[index]['E']}")
            values.setdefault(index, {})[op_name] = value
    for index, state in enumerate(states):
        if index in values:
            state["expectation"] = values[index]


def run_nmax(args, nmax: int):
    work = args.workdir / f"Nmax{nmax}"
    work.mkdir(parents=True, exist_ok=True)
    hamiltonian = work / "hamiltonian_s.bin"   # BIGSTICK takes names starting with "no2" for the format code
    link(hamiltonian, args.no2bpack)
    name = f"Z{args.Z}N{args.N}_Nmax{nmax}"
    start = time.time()
    # Lanczos ("ld"); small spaces break it (a closed-shell Nmax=0 determinant
    # fails outright, the 81-dimensional O20 Nmax=0 space gave a non-J-pure
    # state 26 MeV below the exact ground state), so spaces up to EXACT_MAX
    # use full diagonalization ("ex").
    for method, solver in (("lanczos", ["ld", f"{args.states} {args.max_iter}"]),
                           ("exact", ["ex", str(args.states)])):
        deck = "\n".join(["n", name, "auto", str(args.emax), f"{args.Z} {args.N}", str(args.jz2),
                          args.parity, "y", str(nmax), "no2", hamiltonian.name, f"{args.hw} 0", "end",
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
    if args.operator:
        expectation_values(args, work, name, states)
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
    parser.add_argument("--operator", action="append", default=[], metavar="NAME=FILE",
                        help="scalar operator no2bpack exported with H (same basis); repeatable")
    parser.add_argument("--workdir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    operators = []
    for item in args.operator:
        op_name, sep, op_file = item.partition("=")
        if not sep or not op_name.isidentifier() or not Path(op_file).is_file():
            parser.error(f"--operator needs NAME=FILE with an existing file, got {item}")
        operators.append((op_name, Path(op_file)))
    args.operator = operators
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
        if "expectation" in s[0]:
            print("         ground state: " + ", ".join(f"<{k}>={v:.6f}" for k, v in s[0]["expectation"].items()))
    if summary["extrapolated_ground_state"] is not None:
        print(f"exponential extrapolation (last three Nmax): {summary['extrapolated_ground_state']:.3f}")


if __name__ == "__main__":
    main()
