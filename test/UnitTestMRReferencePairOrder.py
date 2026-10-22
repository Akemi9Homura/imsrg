#!/usr/bin/env python3

"""MRReference::ReadBinary must apply the exchange phase to reversed pairs.

The jref file lists each two-body pair as (a, b) in file orbit numbering.  When
that pair maps to b < a in the ModelSpace, |ba;J> = (-1)^{j_a+j_b+J+1} |ab;J>.
This test rewrites the He4 Nrefmax=2 reference with every a != b pair reversed
and the lambda2 rows/columns multiplied by that phase, i.e. the same physical
cumulant written in the other pair order, and requires both files to read back
to identical Lambda2.  Without the phase, entries between reversed and
unreversed pairs change sign while the symmetry check still passes.
"""

from pathlib import Path
import struct
import sys
import tempfile

import numpy as np
import pyIMSRG

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))

from prototype.mrimsrg.export_jref import export_reference  # noqa: E402

HEADER = struct.Struct("<16sII7idQQ")
ORBIT = struct.Struct("<5id")
CHANNEL = struct.Struct("<3iIQ")
PAIR = struct.Struct("<II")


def reverse_pairs(source: Path, target: Path) -> int:
    data = source.read_bytes()
    out = bytearray()
    pos = 0
    header = HEADER.unpack_from(data, pos)
    norbits, nchannels = header[-2], header[-1]
    pos += HEADER.size + 3 * 64
    j2 = []
    for _ in range(norbits):
        orbit = ORBIT.unpack_from(data, pos)
        j2.append(orbit[3])
        pos += ORBIT.size
    pos += 8 * norbits * norbits
    out += data[:pos]
    reversed_count = 0
    for _ in range(nchannels):
        J, parity, Tz, reserved, npairs = CHANNEL.unpack_from(data, pos)
        out += data[pos:pos + CHANNEL.size]
        pos += CHANNEL.size
        phase = np.ones(npairs)
        for i in range(npairs):
            a, b = PAIR.unpack_from(data, pos)
            pos += PAIR.size
            if a != b:
                out += PAIR.pack(b, a)
                jsum = (j2[a] + j2[b]) // 2
                phase[i] = 1.0 if (jsum + J + 1) % 2 == 0 else -1.0
                reversed_count += 1
            else:
                out += PAIR.pack(a, b)
        size = 8 * npairs * npairs
        matrix = np.frombuffer(data[pos:pos + size], dtype="<f8").reshape(npairs, npairs)
        pos += size
        out += np.ascontiguousarray(phase[:, None] * matrix * phase[None, :], dtype="<f8").tobytes()
    assert pos == len(data)
    target.write_bytes(bytes(out))
    return reversed_count


def read(path: Path):
    ms = pyIMSRG.ModelSpace(2, "He4", "He4")
    ms.SetReferenceOcc(pyIMSRG.MRReference.ReadOccupationMap(ms, str(path)))
    return ms, pyIMSRG.MRReference.ReadBinary(ms, str(path))


with tempfile.TemporaryDirectory() as tmp:
    original = Path(tmp) / "He4_Nrefmax2.jref"
    export_reference(REPOSITORY / "prototype/mrimsrg/data/He4_Nrefmax2", original)
    swapped = Path(tmp) / "He4_Nrefmax2_reversed_pairs.jref"
    count = reverse_pairs(original, swapped)
    assert count > 0
    ms_a, ref_a = read(original)
    ms_b, ref_b = read(swapped)
    worst = 0.0
    nonzero_sign_sensitive = 0
    for ch in range(ms_a.GetNumberTwoBodyChannels()):
        tbc = ms_a.GetTwoBodyChannel(ch)
        n = tbc.GetNumberKets()
        for i in range(n):
            for j in range(n):
                va = ref_a.Lambda2.GetTBME_norm_chij(ch, ch, i, j)
                vb = ref_b.Lambda2.GetTBME_norm_chij(ch, ch, i, j)
                worst = max(worst, abs(va - vb))
                if abs(va) > 1e-8:
                    nonzero_sign_sensitive += 1
    assert nonzero_sign_sensitive > 0
    assert worst < 1e-14, worst
    print(f"reversed {count} pairs; max |Lambda2 difference| = {worst:.3e}")
