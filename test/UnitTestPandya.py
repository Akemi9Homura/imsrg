#!/usr/bin/env python3
"""Validate Pandya and zeroed intermediates against independent contractions."""
import itertools
import math
import pyIMSRG as imsrg

ms = imsrg.ModelSpace(1, "He6", "He6")
ms.PreCalculateSixJ()
ms.PreCalculateNineJ()
ut = imsrg.UnitTest(ms)
ut.SetRandomSeed(222)
largest = 0.0
count = 0


def check(name, x, y, initial, references):
    global largest, count
    imsrg.Commutator.SetSingleThread(True)
    expected = imsrg.Operator(initial)
    for reference in references:
        getattr(imsrg.ReferenceImplementations, reference)(x, y, expected)
    serial = None
    for single in (True, False):
        imsrg.Commutator.SetSingleThread(single)
        actual = imsrg.Operator(initial)
        getattr(imsrg.Commutator, name)(x, y, actual)
        error = (actual - expected).Norm()
        assert math.isfinite(error) and error < 1e-10, (name, count, single, error)
        largest = max(largest, error)
        if serial is not None:
            assert (actual - serial).Norm() < 1e-12
        serial = actual
    count += 1


try:
    for j, parity, trank, hx, hy in itertools.product((0, 1, 2), (0, 1), (0, 1, 2), (-1, 1), (-1, 1)):
        x = ut.RandomOp(ms, 0, 0, 0, 2, hx)
        y = ut.RandomOp(ms, j, trank, parity, 2, hy)
        initial = ut.RandomOp(ms, j, trank, parity, 2, -hx * hy)
        y.MakeReduced()
        initial.MakeReduced()
        check("comm222_phst", x, y, initial, ["comm222_phst"])
        check("comm222_pp_hh_221st", x, y, initial, ["comm222_pp_hhst", "comm221st"])
    for hx, hy in itertools.product((-1, 1), repeat=2):
        x = ut.RandomOp(ms, 0, 0, 0, 2, hx)
        y = ut.RandomOp(ms, 0, 0, 0, 2, hy)
        initial = ut.RandomOp(ms, 0, 0, 0, 2, -hx * hy)
        check("comm222_phss", x, y, initial, ["comm222_phss"])
        check("comm222_pp_hh_221ss", x, y, initial, ["comm222_pp_hhss", "comm221ss"])

    # Include distinct radial orbits with the same angular quantum numbers.
    radial_ms = imsrg.ModelSpace(2, "He6", "He6")
    radial_ms.PreCalculateSixJ()
    radial_ms.PreCalculateNineJ()
    radial_ut = imsrg.UnitTest(radial_ms)
    radial_ut.SetRandomSeed(2222)
    x = radial_ut.RandomOp(radial_ms, 0, 0, 0, 2, -1)
    y = radial_ut.RandomOp(radial_ms, 1, 1, 1, 2, 1)
    initial = radial_ut.RandomOp(radial_ms, 1, 1, 1, 2, 1)
    check("comm222_phst", x, y, initial, ["comm222_phst"])
    check("comm222_pp_hh_221st", x, y, initial, ["comm222_pp_hhst", "comm221st"])
finally:
    imsrg.Commutator.SetSingleThread(False)

print(f"Pandya/intermediates: {count} cases, serial/parallel; max norm error {largest:.3e}")
