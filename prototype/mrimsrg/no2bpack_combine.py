#!/usr/bin/env python3
"""Linear combination of no2bpack files with the same layout.

  no2bpack_combine.py --term 1 H.no2bpack --term 0.5 H_HCM.no2bpack [--constant -15] --output out.no2bpack

writes  sum_i c_i X_i + constant  (the constant goes into the zero body).  All
inputs must come from the same export (same orbits, same record sequence);
this is checked.  Use: a Lawson term beta (Hcm - 3/2 hw) for BIGSTICK, or
H +- beta O for a Hellmann-Feynman check of operator expectation values.

Layout (ReadWrite::Write_no2bpack): header double hw, int emax, norbits, nobme,
ntbme; norbits x 4 ints (n, l, 2j, 2tz); double zero body; nobme float32;
ntbme records of 5 ints (a, b, c, d, J) and a float32, followed by a float32
centre-of-mass placeholder when |f_ab - f_cd| is 0 or 2 (f = sum of 2n + l).
"""

from __future__ import annotations

import argparse
from pathlib import Path
import struct

import numpy as np


def layout(data: bytes):
    """Byte offsets of the zero body, the one-body floats and the two-body floats."""
    hw, emax, norbits, nobme, ntbme = struct.unpack_from("<d4i", data, 0)
    pos = 24
    orbits = np.frombuffer(data, dtype="<i4", count=4 * norbits, offset=pos).reshape(norbits, 4)
    pos += 16 * norbits
    zero = pos
    pos += 8
    one = pos
    pos += 4 * nobme
    f = 2 * orbits[:, 0] + orbits[:, 1]
    two = np.empty(ntbme, dtype=np.int64)
    for i in range(ntbme):
        a, b, c, d = struct.unpack_from("<4i", data, pos)
        two[i] = pos + 20
        pos += 24
        if abs(int(f[a] + f[b] - f[c] - f[d])) in (0, 2):
            pos += 4
    if pos != len(data):
        raise ValueError(f"no2bpack layout ends at byte {pos}, file has {len(data)}")
    return {"header": data[:24 + 16 * norbits], "zero": zero, "one": one, "nobme": nobme, "two": two}


def floats_at(data: bytes, offsets: np.ndarray) -> np.ndarray:
    raw = np.frombuffer(data, dtype=np.uint8)
    idx = offsets[:, None] + np.arange(4)
    return raw[idx].copy().view("<f4").ravel()


def combine(terms, constant: float, output: Path):
    blobs = [(float(c), Path(p).read_bytes()) for c, p in terms]
    ref = layout(blobs[0][1])
    # every byte that is not a matrix element (header, orbits, record indices,
    # centre-of-mass placeholders) must agree between the inputs
    structure = np.ones(len(blobs[0][1]), dtype=bool)
    structure[ref["zero"]:ref["zero"] + 8] = False
    structure[ref["one"]:ref["one"] + 4 * ref["nobme"]] = False
    structure[(ref["two"][:, None] + np.arange(4)).ravel()] = False
    first = np.frombuffer(blobs[0][1], dtype=np.uint8)
    out = bytearray(blobs[0][1])
    zero = constant
    one = np.zeros(ref["nobme"], dtype=np.float64)
    two = np.zeros(len(ref["two"]), dtype=np.float64)
    for coefficient, data in blobs:
        if len(data) != len(first) or not np.array_equal(np.frombuffer(data, dtype=np.uint8)[structure],
                                                          first[structure]):
            raise ValueError("inputs do not share the same no2bpack layout (orbits, records)")
        zero += coefficient * struct.unpack_from("<d", data, ref["zero"])[0]
        one += coefficient * np.frombuffer(data, dtype="<f4", count=ref["nobme"], offset=ref["one"])
        two += coefficient * floats_at(data, ref["two"])
    struct.pack_into("<d", out, ref["zero"], zero)
    out[ref["one"]:ref["one"] + 4 * ref["nobme"]] = one.astype("<f4").tobytes()
    raw = np.frombuffer(out, dtype=np.uint8).copy()
    raw[ref["two"][:, None] + np.arange(4)] = two.astype("<f4").view(np.uint8).reshape(-1, 4)
    output.write_bytes(raw.tobytes())


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--term", nargs=2, action="append", required=True, metavar=("COEFFICIENT", "FILE"))
    parser.add_argument("--constant", type=float, default=0.0, help="added to the zero body")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    combine(args.term, args.constant, args.output)


if __name__ == "__main__":
    main()
