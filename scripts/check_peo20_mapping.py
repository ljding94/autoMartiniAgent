"""Sanity-check the PEO-20 AA→CG mapping against a single AA frame.

Projects equil.gro through reference/polyply_PEO20/PEO20_mapping.json,
asserts each bead is the mass-weighted COM of its atoms, prints the
20 bead positions, and reports the 1-2 bond lengths (should sit near the
Polyply equilibrium 0.36 nm) and 1-2-3 bond angles (Polyply equilibrium
123 deg). A single frame won't match perfectly, but anything wildly
off would expose a wrong atom-index grouping.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MAPPING = REPO / "reference/polyply_PEO20/PEO20_mapping.json"
GRO = REPO / "reference/gromacs/equil.gro"
ATOMIC_MASS = {"O": 15.9994, "C": 12.011, "H": 1.008, "N": 14.007, "S": 32.06}


def parse_gro(gro: Path):
    text = gro.read_text().splitlines()
    n = int(text[1].strip())
    coords = {}
    for line in text[2 : 2 + n]:
        resnr = int(line[0:5])
        resname = line[5:10].strip()
        aname = line[10:15].strip()
        idx = int(line[15:20])
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        coords[idx] = {
            "resnr": resnr,
            "resname": resname,
            "name": aname,
            "xyz": (x, y, z),
        }
    return coords


def element_of(atom_name: str) -> str:
    return atom_name[0].upper()


def project(mapping, coords):
    beads = []
    for b in mapping["beads"]:
        ms, mx, my, mz = 0.0, 0.0, 0.0, 0.0
        for idx, name in zip(b["atom_indices"], b["atom_names"]):
            x, y, z = coords[idx]["xyz"]
            m = ATOMIC_MASS[element_of(name)]
            ms += m
            mx += m * x
            my += m * y
            mz += m * z
        beads.append(
            {
                "bead_id": b["bead_id"],
                "mass": ms,
                "xyz": (mx / ms, my / ms, mz / ms),
            }
        )
    return beads


def dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def angle_deg(a, b, c):
    va = [a[i] - b[i] for i in range(3)]
    vc = [c[i] - b[i] for i in range(3)]
    na = math.sqrt(sum(x * x for x in va))
    nc = math.sqrt(sum(x * x for x in vc))
    cos = sum(va[i] * vc[i] for i in range(3)) / (na * nc)
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))


def main():
    mapping = json.loads(MAPPING.read_text())
    coords = parse_gro(GRO)
    beads = project(mapping, coords)

    print(f"AA atoms read from {GRO.name}: {len(coords)}")
    print(f"PEO atoms used (resname=ETHOX): "
          f"{sum(1 for v in coords.values() if v['resname'] == 'ETHOX')}")
    print()
    print(f"{'bead':>4}  {'mass':>7}  {'x':>7}  {'y':>7}  {'z':>7}")
    for b in beads:
        x, y, z = b["xyz"]
        print(f"{b['bead_id']:>4d}  {b['mass']:>7.3f}  {x:>7.3f}  {y:>7.3f}  {z:>7.3f}")

    bonds = [dist(beads[i]["xyz"], beads[i + 1]["xyz"]) for i in range(19)]
    angles = [
        angle_deg(beads[i]["xyz"], beads[i + 1]["xyz"], beads[i + 2]["xyz"])
        for i in range(18)
    ]
    print()
    print(f"1-2 bond lengths (nm): "
          f"mean={sum(bonds) / len(bonds):.3f}  "
          f"min={min(bonds):.3f}  max={max(bonds):.3f}")
    print(f"  Polyply equilibrium: 0.360 nm (single frame, expect spread)")
    print(f"1-2-3 bond angles (deg): "
          f"mean={sum(angles) / len(angles):.1f}  "
          f"min={min(angles):.1f}  max={max(angles):.1f}")
    print(f"  Polyply equilibrium: 123 deg")

    masses = [b["mass"] for b in beads]
    expected = [45.061] + [44.053] * 18 + [45.061]
    bad = [
        (i + 1, m, e)
        for i, (m, e) in enumerate(zip(masses, expected))
        if abs(m - e) > 0.01
    ]
    if bad:
        print()
        print("MASS MISMATCH:")
        for bid, got, exp in bad:
            print(f"  bead {bid}: got {got:.3f}, expected {exp:.3f}")
        raise SystemExit(1)
    print()
    print("OK: all bead masses match expected AA mass sums.")


if __name__ == "__main__":
    main()
