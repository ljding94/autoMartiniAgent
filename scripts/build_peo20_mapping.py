"""Build AA→CG mapping for PEO-20 (CHARMM36 ETHOX × 20 → Polyply 20×SN3r).

Reads the AA topology at reference/gromacs/toppar/S1P1.itp and emits:
  - reference/polyply_PEO20/PEO20_mapping.json  (canonical, atom-index keyed)
  - reference/polyply_PEO20/PEO20.map           (Martini-style, atom-name keyed)

Mapping rule (canonical PEO Martini 3 cold start): one SN3r bead per
ETHOX residue. End-group atoms (HO1 on residue 1; the extra methyl-H on
residue 20) fold into their residue's bead. No splitting, no merging.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
AA_ITP = REPO / "reference/gromacs/toppar/S1P1.itp"
CG_ITP = REPO / "reference/polyply_PEO20/PEO20.itp"
OUT_JSON = REPO / "reference/polyply_PEO20/PEO20_mapping.json"
OUT_MAP = REPO / "reference/polyply_PEO20/PEO20.map"

ATOM_LINE = re.compile(
    r"^\s*(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+\d+\s+\S+\s+(\S+)"
)


def parse_aa_atoms(itp: Path):
    atoms = []
    in_block = False
    for raw in itp.read_text().splitlines():
        line = raw.split(";", 1)[0].rstrip()
        if not line.strip():
            continue
        if line.lstrip().startswith("["):
            in_block = line.strip() == "[ atoms ]"
            continue
        if not in_block:
            continue
        m = ATOM_LINE.match(line)
        if not m:
            continue
        idx, atype, resnr, resname, aname, mass = m.groups()
        atoms.append(
            {
                "index": int(idx),
                "type": atype,
                "resnr": int(resnr),
                "resname": resname,
                "name": aname,
                "mass": float(mass),
            }
        )
    return atoms


def parse_cg_beads(itp: Path):
    beads = []
    in_block = False
    for raw in itp.read_text().splitlines():
        line = raw.split(";", 1)[0].rstrip()
        if not line.strip():
            continue
        if line.lstrip().startswith("["):
            in_block = line.strip() == "[ atoms ]"
            continue
        if not in_block:
            continue
        parts = line.split()
        if len(parts) < 7:
            continue
        beads.append(
            {
                "bead_id": int(parts[0]),
                "bead_type": parts[1],
                "resnr": int(parts[2]),
                "resname": parts[3],
                "bead_name": parts[4],
                "mass": float(parts[6]),
            }
        )
    return beads


def build_mapping(aa_atoms, cg_beads):
    aa_by_resnr = {}
    for a in aa_atoms:
        aa_by_resnr.setdefault(a["resnr"], []).append(a)

    assert len(cg_beads) == len(aa_by_resnr) == 20, (
        f"expected 20 beads and 20 residues; got "
        f"{len(cg_beads)} beads, {len(aa_by_resnr)} residues"
    )

    beads_out = []
    for bead in cg_beads:
        resnr = bead["bead_id"]
        residue_atoms = aa_by_resnr[resnr]
        heavy = [a for a in residue_atoms if not a["name"].startswith("H")]
        bead_mass_aa = sum(a["mass"] for a in residue_atoms)
        comment = {1: "HO terminus", 20: "CH3 terminus"}.get(resnr, "ether monomer")
        beads_out.append(
            {
                "bead_id": bead["bead_id"],
                "bead_name": bead["bead_name"],
                "bead_type": bead["bead_type"],
                "cg_mass": bead["mass"],
                "aa_residue": resnr,
                "aa_resname": residue_atoms[0]["resname"],
                "atom_indices": [a["index"] for a in residue_atoms],
                "atom_names": [a["name"] for a in residue_atoms],
                "heavy_atom_indices": [a["index"] for a in heavy],
                "heavy_atom_count": len(heavy),
                "aa_mass_sum": round(bead_mass_aa, 4),
                "comment": comment,
            }
        )
    return beads_out


def write_json(beads):
    payload = {
        "molecule": "PEO20",
        "aa_topology": "reference/gromacs/toppar/S1P1.itp",
        "aa_residue_name": "ETHOX",
        "aa_total_atoms": 142,
        "cg_topology": "reference/polyply_PEO20/PEO20.itp",
        "cg_bead_count": 20,
        "weighting": "mass",
        "rule": "one SN3r bead per ETHOX residue; bead i ← residue i (all heavy + bonded H)",
        "provenance": {
            "cg_backend": "polyply",
            "cg_command": (
                "polyply gen_params -lib martini3 -seq PEO:20 "
                "-o PEO20.itp -name PEO20"
            ),
            "aa_force_field": "CHARMM36 General (via CHARMM-GUI)",
            "aa_water": "TIP3P",
            "aa_trajectory": "reference/gromacs/step5_200.xtc",
            "mapping_origin": (
                "canonical PEO Martini 3 cold start — residue-to-bead 1:1, "
                "end-group atoms fold into their residue's bead"
            ),
        },
        "beads": beads,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2) + "\n")


def write_map(beads):
    lines = [
        "; PEO-20 AA→CG mapping (CHARMM36 ETHOX → Martini 3 SN3r).",
        "; Generated by scripts/build_peo20_mapping.py — do not hand-edit.",
        "",
        "[ molecule ]",
        "PEO20",
        "",
        "[ martini ]",
        " ".join(b["bead_name"] for b in beads),
        "",
        "[ mapping ]",
        "charmm36",
        "",
        "[ atoms ]",
        ";  idx  name   bead",
    ]
    for b in beads:
        for atom_idx, atom_name in zip(b["atom_indices"], b["atom_names"]):
            lines.append(f"  {atom_idx:4d}  {atom_name:<5s}  {b['bead_id']:>2d}")
    lines.append("")
    OUT_MAP.write_text("\n".join(lines))


def main():
    aa_atoms = parse_aa_atoms(AA_ITP)
    cg_beads = parse_cg_beads(CG_ITP)
    beads = build_mapping(aa_atoms, cg_beads)
    write_json(beads)
    write_map(beads)
    print(f"wrote {OUT_JSON.relative_to(REPO)}")
    print(f"wrote {OUT_MAP.relative_to(REPO)}")
    print()
    print("bead summary:")
    for b in beads:
        print(
            f"  bead {b['bead_id']:>2d} ({b['bead_type']}) "
            f"← residue {b['aa_residue']:>2d}  "
            f"{b['heavy_atom_count']} heavy, {len(b['atom_indices'])} total, "
            f"AA mass {b['aa_mass_sum']:.3f} g/mol  [{b['comment']}]"
        )


if __name__ == "__main__":
    main()
