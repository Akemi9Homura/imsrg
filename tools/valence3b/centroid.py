#!/usr/bin/env python3
"""Configuration centroids of a valence Hamiltonian with and without the valence three-body (level-1 diagnostic).

For a configuration (occupation n_a of every valence orbit) the centroid is the average of <D|H|D> over all
determinants D of the configuration.  For a k-body operator it equals sum over orbit multisets {a,b,..} of
(number of k-particle subsets with those orbits) x (average diagonal element on the antisymmetric k-particle
states of those orbits), so no many-body solver is needed.  The difference "with W - without W" is what the
valence three-body (with its lower-body terms relative to the core) does to the cost of each configuration.

  centroid.py <run.snt> <v3b.snt> <v3b_3b.int> <config> [<config> ...]
  config: name=orbit:n,orbit:n,...  with orbit = p|n + n l 2j, e.g.  40Ca_0p0h=p1s1:2,p0d3:4,n1s1:2,n0d3:4
"""

from itertools import combinations
from math import comb
import sys

import numpy as np

import pspace


def averages(snt, three_entries=None):
    """Average diagonal elements per orbit, orbit pair and orbit triple (keys: sorted .snt orbit indices)."""
    space = pspace.Space(snt["orbits"])
    one, two, three = pspace.hamiltonian(snt, space, three_entries)
    orb = [o for o, m in space.mstates]
    avg1, avg2, avg3 = {}, {}, {}
    for o in snt["orbits"]:
        ks = [k for k, (oo, m) in enumerate(space.mstates) if oo == o]
        avg1[o] = float(np.mean([one.get((k, k), 0.0) for k in ks]))
    sums2, sums3 = {}, {}
    for (bra, ket), v in two.items():
        if bra == ket:
            key = tuple(sorted(orb[k] for k in bra))
            sums2[key] = sums2.get(key, 0.0) + v
    if three:
        for (bra, ket), v in three.items():
            if bra == ket:
                key = tuple(sorted(orb[k] for k in bra))
                sums3[key] = sums3.get(key, 0.0) + v
    deg = {o: snt["orbits"][o][2] + 1 for o in snt["orbits"]}

    def dimension(key):
        d = 1
        for o in set(key):
            d *= comb(deg[o], key.count(o))
        return d
    for key, s in sums2.items():
        avg2[key] = s / dimension(key)
    for key, s in sums3.items():
        avg3[key] = s / dimension(key)
    return avg1, avg2, avg3


def centroid(snt, avg, occ):
    """occ: {orbit index: n}."""
    avg1, avg2, avg3 = avg
    e = snt["zero"] + sum(n * avg1[o] for o, n in occ.items())
    for key, v in avg2.items():
        e += v * _subsets(key, occ)
    for key, v in avg3.items():
        e += v * _subsets(key, occ)
    return e


def _subsets(key, occ):
    count = 1
    for o in set(key):
        count *= comb(int(round(occ.get(o, 0))), key.count(o))
    return count


def parse_config(text, orbits):
    name, spec = text.split("=")
    label = {}
    for k, (n, l, j2, tz2) in orbits.items():
        label[f"{'p' if tz2 < 0 else 'n'}{n}{'spdfgh'[l]}{j2}"] = k
    occ = {}
    for item in spec.split(","):
        o, n = item.split(":")
        occ[label[o]] = int(n)
    return name, occ


def main():
    std_path, v3b_path, three_path, *configs = sys.argv[1:]
    std, v3b = pspace.read_snt(std_path), pspace.read_snt(v3b_path)
    a_std = averages(std)
    a_v3b = averages(v3b, pspace.read_3b(three_path, std["orbits"]))
    print(f"{'configuration':18s} {'without W':>12s} {'with W':>12s} {'W shift':>10s}")
    for text in configs:
        name, occ = parse_config(text, std["orbits"])
        e0, e1 = centroid(std, a_std, occ), centroid(v3b, a_v3b, occ)
        print(f"{name:18s} {e0:12.4f} {e1:12.4f} {e1 - e0:+10.4f}")


if __name__ == "__main__":
    main()
