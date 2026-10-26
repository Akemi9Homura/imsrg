#!/usr/bin/env python3

"""BIGSTICK wavefunction -> jref agrees with the simpleFCI -> jref route.

The He4 Nrefmax=2 ground state from BIGSTICK (test/data/bigstick_He4_Nrefmax2,
single-precision Lanczos vectors) is converted by bigstick_jref.py and read by
the C++ MRReference reader with the production validation tolerance.  Natural
occupations, the natural-orbit transformation and lambda2 must match the jref
built from the double-precision simpleFCI fixture to the BIGSTICK precision.
"""

from pathlib import Path
import struct
import sys
import tempfile

import numpy as np
import pyIMSRG

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))

from prototype.mrimsrg.bigstick_jref import convert  # noqa: E402
from prototype.mrimsrg.export_jref import export_reference  # noqa: E402

FIXTURE = REPOSITORY / "test/data/bigstick_He4_Nrefmax2"


def read(path):
    modelspace = pyIMSRG.ModelSpace(2, "He4", "He4")
    modelspace.SetReferenceOcc(pyIMSRG.MRReference.ReadOccupationMap(modelspace, str(path)))
    reference = pyIMSRG.MRReference.ReadBinary(modelspace, str(path), 1e-10)
    occupations = np.array([modelspace.GetOrbit(i).occ for i in range(modelspace.GetNumberOrbits())])
    return modelspace, reference, occupations


def transformation(path):
    """(orbit labels, natural-orbit matrix) straight from the jref file."""
    data = path.read_bytes()
    pos = 16 + 8 + 28
    _hw, norbits, _nchannels = struct.unpack_from("<dQQ", data, pos)
    pos += 24 + 3 * 64
    labels = []
    for _ in range(norbits):
        labels.append(struct.unpack_from("<5i", data, pos)[1:])
        pos += 28
    return labels, np.frombuffer(data, "<f8", norbits * norbits, pos).reshape(norbits, norbits)


with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    bigstick = root / "bigstick.jref"
    simplefci = root / "simplefci.jref"
    summary = convert(FIXTURE / "he4_nref2.bas", FIXTURE / "he4_nref2.wfn", FIXTURE / "he4_nref2.wfn",
                      bigstick, nrefmax=2, hw=20.0)
    export_reference(REPOSITORY / "prototype/mrimsrg/data/He4_Nrefmax2", simplefci)
    assert abs(summary["energy"] + 17.0134908) < 1e-6, summary
    assert abs(summary["pairs"] - 6.0) < 1e-12, summary

    ms_b, ref_b, occ_b = read(bigstick)
    ms_s, ref_s, occ_s = read(simplefci)
    occupation_error = float(np.max(np.abs(occ_b - occ_s)))

    labels_b, C_b = transformation(bigstick)
    labels_s, C_s = transformation(simplefci)
    permutation = [labels_b.index(label) for label in labels_s]
    transformation_error = float(np.max(np.abs(C_s - C_b[np.ix_(permutation, permutation)])))

    lambda_error = 0.0
    n = ms_b.GetNumberOrbits()
    for a in range(n):
        for b in range(n):
            ja, jb = ms_b.GetOrbit(a).j2, ms_b.GetOrbit(b).j2
            for c in range(n):
                for d in range(n):
                    for J in range(abs(ja - jb) // 2, (ja + jb) // 2 + 1):
                        lambda_error = max(lambda_error, abs(ref_b.Lambda2.GetTBME_J(J, J, a, b, c, d)
                                                             - ref_s.Lambda2.GetTBME_J(J, J, a, b, c, d)))

print(f"BIGSTICK vs simpleFCI jref: occupations {occupation_error:.1e}, "
      f"natural orbits {transformation_error:.1e}, lambda2 {lambda_error:.1e}")
assert occupation_error < 1e-7 and transformation_error < 1e-7 and lambda_error < 5e-7
