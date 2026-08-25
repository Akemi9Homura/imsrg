#!/usr/bin/env python3

"""Compare every selected White-NCSM EN denominator with the m-scheme oracle."""

import json
from pathlib import Path
import sys
import tempfile

import numpy as np
import pyIMSRG


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))
sys.path.insert(0, str(REPOSITORY / "test"))

from UnitTestMRCorrelatedDriver import projected_reference_densities  # noqa: E402
from UnitTestMRDriver import vacuum_operator  # noqa: E402
from prototype.mrimsrg.export_jref import export_reference  # noqa: E402
from prototype.mrimsrg.generator import epstein_nesbet_denominators  # noqa: E402
from prototype.mrimsrg.reference_io import load_reference  # noqa: E402
from prototype.mrimsrg.sr_imsrgpp_check import operator_to_mscheme  # noqa: E402


def cutoff(value, threshold=1e-6):
    return threshold if abs(value) < threshold else float(value)


def m_orbits_by_spherical_orbit(orbits):
    result = {}
    for index, row in enumerate(orbits):
        result.setdefault(int(row[0]), []).append(index)
    return result


def one_body_representatives(orbits, i, j):
    i_states = {int(orbits[p, 4]): p for p in range(len(orbits)) if int(orbits[p, 0]) == i}
    j_states = {int(orbits[p, 4]): p for p in range(len(orbits)) if int(orbits[p, 0]) == j}
    common_m = min(set(i_states) & set(j_states))
    return i_states[common_m], j_states[common_m]


def pair_representatives(groups, a, b):
    for p in groups[a]:
        for q in groups[b]:
            if p != q:
                return p, q
    raise AssertionError(f"no antisymmetric m-scheme representative for ({a},{b})")


def main():
    if len(sys.argv) not in (1, 2):
        raise SystemExit("usage: UnitTestMRDenominators.py [report.json]")
    report_path = Path(sys.argv[1]).resolve() if len(sys.argv) == 2 else None
    fixtures = (
        ("He4", "He4_Nrefmax2"),
        ("Be8", "Be8_Nrefmax0_final"),
        ("C12", "C12_Nrefmax0_final"),
        ("O16", "O16_Nrefmax2"),
    )
    report = {
        "schema": "mrimsrg_white_ncsm_denominators_v1",
        "comparison_input": "same float64 J-representable Hamiltonian and RDM",
        "cutoff_mev": 1e-6,
        "tolerance_max_abs_mev": 1e-10,
        "one_body_columns": [
            "source_orbit", "target_orbit", "source_m", "target_m",
            "actual_mev", "expected_mev", "abs_error_mev",
        ],
        "two_body_columns": [
            "J", "parity", "Tz", "channel_index", "source_pair_index",
            "target_pair_index", "a", "b", "c", "d", "p", "q", "r", "s",
            "actual_mev", "expected_mev", "abs_error_mev",
        ],
        "systems": {},
    }
    with tempfile.TemporaryDirectory() as temporary_directory:
        root = Path(temporary_directory)
        for nucleus, fixture in fixtures:
            data = load_reference(REPOSITORY / "prototype/mrimsrg/data" / fixture)
            reference_file = root / f"{nucleus}.jref"
            export_reference(REPOSITORY / "prototype/mrimsrg/data" / fixture, reference_file)
            modelspace = pyIMSRG.ModelSpace(2, nucleus, nucleus)
            modelspace.SetHbarOmega(20.0)
            modelspace.SetReferenceOcc(
                pyIMSRG.MRReference.ReadOccupationMap(modelspace, str(reference_file))
            )
            reference = pyIMSRG.MRReference.ReadBinary(modelspace, str(reference_file))
            vacuum = vacuum_operator(modelspace, data)
            initial = reference.NormalOrder(
                vacuum.TransformOneAndTwoBody(reference.NaturalOrbitTransformation)
            )
            oracle_h = operator_to_mscheme(initial, data.orbits)
            oracle_densities = projected_reference_densities(reference, modelspace, data)
            delta1, delta2 = epstein_nesbet_denominators(
                oracle_h,
                oracle_densities,
                spherical_orbit_groups=np.asarray(data.orbits[:, 0], dtype=np.int64),
            )
            eta = pyIMSRG.Operator(modelspace)
            eta.SetAntiHermitian()
            generator = pyIMSRG.Generator()
            generator.SetType("white-ncsm")
            generator.Update(initial, eta)

            one_body = []
            for i in range(modelspace.GetNumberOrbits()):
                oi = modelspace.GetOrbit(i)
                for j in range(i + 1, modelspace.GetNumberOrbits()):
                    oj = modelspace.GetOrbit(j)
                    if (oi.l, oi.j2, oi.tz2) != (oj.l, oj.j2, oj.tz2):
                        continue
                    if 2 * oi.n + oi.l == 2 * oj.n + oj.l:
                        continue
                    p, q = one_body_representatives(data.orbits, i, j)
                    for source, target, mp, mq in ((i, j, p, q), (j, i, q, p)):
                        actual = generator.Get1bDenominatorWhiteNCSM(source, target)
                        expected = cutoff(delta1[mp, mq])
                        one_body.append([
                            source, target, mp, mq, actual, expected,
                            abs(actual - expected),
                        ])

            groups = m_orbits_by_spherical_orbit(data.orbits)
            two_body = []
            for channel_index in range(modelspace.GetNumberTwoBodyChannels()):
                channel = modelspace.GetTwoBodyChannel(channel_index)
                for ibra in range(channel.GetNumberKets()):
                    bra = channel.GetKet(ibra)
                    oi = modelspace.GetOrbit(bra.p)
                    oj = modelspace.GetOrbit(bra.q)
                    bra_quanta = 2 * oi.n + oi.l + 2 * oj.n + oj.l
                    for iket in range(ibra + 1, channel.GetNumberKets()):
                        ket = channel.GetKet(iket)
                        ok = modelspace.GetOrbit(ket.p)
                        ol = modelspace.GetOrbit(ket.q)
                        ket_quanta = 2 * ok.n + ok.l + 2 * ol.n + ol.l
                        if bra_quanta == ket_quanta:
                            continue
                        p, q = pair_representatives(groups, bra.p, bra.q)
                        r, s = pair_representatives(groups, ket.p, ket.q)
                        for source, target, indices in (
                            (ibra, iket, (p, q, r, s)),
                            (iket, ibra, (r, s, p, q)),
                        ):
                            actual = generator.Get2bDenominatorWhiteNCSM(
                                channel_index, source, target
                            )
                            expected = cutoff(delta2[indices])
                            two_body.append([
                                channel.J, channel.parity, channel.Tz,
                                channel_index, source, target,
                                bra.p, bra.q, ket.p, ket.q,
                                *indices, actual, expected, abs(actual - expected),
                            ])

            maximum = max(
                [entry[-1] for entry in one_body + two_body],
                default=0.0,
            )
            print(
                f"{nucleus}: 1B={len(one_body)} 2B={len(two_body)} "
                f"max={maximum:.3e} MeV"
            )
            assert maximum < 1e-10
            report["systems"][nucleus] = {
                "fixture": fixture,
                "maximum_abs_error_mev": maximum,
                "one_body": one_body,
                "two_body": two_body,
            }
    if report_path is not None:
        report_path.write_text(
            json.dumps(report, separators=(",", ":"), sort_keys=True) + "\n",
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
