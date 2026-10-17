#!/usr/bin/env python3
"""Check scalar-tensor 122 contractions against the independent reference."""
import itertools
import math

import pyIMSRG


ms = pyIMSRG.ModelSpace(2, "He6", "He6")
ms.PreCalculateSixJ()
ut = pyIMSRG.UnitTest(ms)
ut.SetRandomSeed(122)
largest_error = 0.0


def check(x, y, initial, case):
    reference = pyIMSRG.Operator(initial)
    pyIMSRG.Commutator.SetSingleThread(True)
    pyIMSRG.ReferenceImplementations.comm122st(x, y, reference)
    serial = None
    largest = 0.0
    for single_thread in (True, False):
        pyIMSRG.Commutator.SetSingleThread(single_thread)
        actual = pyIMSRG.Operator(initial)
        pyIMSRG.Commutator.comm122st(x, y, actual)
        error = (actual - reference).Norm()
        assert math.isfinite(error) and error < 1e-10, (
            case, single_thread, error
        )
        largest = max(largest, error)
        if serial is None:
            serial = actual
        else:
            assert (actual - serial).Norm() == 0.0
    return largest

try:
    for j, parity, trank, hx, hy in itertools.product(
        (1, 2), (0, 1), (0, 1), (-1, 1), (-1, 1)
    ):
        x = ut.RandomOp(ms, 0, 0, 0, 2, hx)
        y = ut.RandomOp(ms, j, trank, parity, 2, hy)
        x.MakeNotReduced()
        y.MakeReduced()
        # Nonzero output also checks that the contraction adds to the result.
        initial = ut.RandomOp(ms, j, trank, parity, 2, -hx * hy)
        largest_error = max(largest_error, check(x, y, initial, (j, parity, trank, hx, hy)))

    # Exercise a channel spanning several temporary ket blocks, including
    # its final partial block, with the symmetry used in Magnus evolution.
    large_ms = pyIMSRG.ModelSpace(6, "O18", "sd-shell")
    large_ms.PreCalculateSixJ()
    max_kets = max(large_ms.GetTwoBodyChannel(ch).GetNumberKets()
                   for ch in range(large_ms.GetNumberTwoBodyChannels()))
    assert max_kets > 256
    large_ut = pyIMSRG.UnitTest(large_ms)
    large_ut.SetRandomSeed(1226)
    x = large_ut.RandomOp(large_ms, 0, 0, 0, 2, -1)
    y = large_ut.RandomOp(large_ms, 1, 0, 0, 2, 1)
    x.MakeNotReduced()
    y.MakeReduced()
    initial = 0.0 * y
    largest_error = max(largest_error, check(x, y, initial, "large channel"))
finally:
    pyIMSRG.Commutator.SetSingleThread(False)

print(f"comm122st: 32 tensor/symmetry cases plus large-channel case, serial and parallel; max error {largest_error:.3e}")
