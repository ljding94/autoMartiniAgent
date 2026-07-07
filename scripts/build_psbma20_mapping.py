"""Derive the PSBMA-20mer AA→CG mapping from the AA .gro + the single-monomer
PSBA.itp template.

Approach
--------
1. **Heavy atoms** map by a static table derived by aligning the SMILES
   ``CC(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]`` (auto-martiniM3's SMILES
   atom indices, per ``tests/fixtures/psbma/PSBA.itp`` bead-atom comments)
   with the AA atom names shipped in ``reference/PSBMA_20mer/*.gro``.

   ``PSBA.itp`` bead partition:
     bead 1  SC1  (C01, "CCC"):        SMILES C0,C1,C2      → AA C1,C2,C3
     bead 2  TN5a (N01, "C=O"):        SMILES C3,O4         → AA C4,O1
     bead 3  TP2a (P01, "CO"):         SMILES O5,C6         → AA O2,C5
     bead 4  Q1   (101, "C[N+](C)C"):  SMILES C7,N8,C9,C10  → AA C6,N,C7,C8
     bead 5  SC1  (C02, "CCC"):        SMILES C11,C12,C13   → AA C9,C10,C11
     bead 6  Q1   (102, "O=[SH](=O)[O-]"): SMILES S14,O15,O16,O17 → AA S,O3,O4,O5

2. **Hydrogens** are assigned by spatial nearest-neighbour on the reference
   ``.gro`` frame: each H is folded into whichever heavy atom it is closest
   to, and inherits that heavy atom's bead. Deterministic and topology-free.

3. **End groups** — residue 1 has an extra ``H6`` (initiator-side proton on
   the terminal backbone C); residue 20 has an extra ``BR1`` (ATRP terminator
   Br). Both are additional atoms in the .gro. They fold into their monomer's
   backbone (bead 1) via the same nearest-heavy rule.

Outputs
-------
- ``reference/PSBMA_20mer/PSBMA20_mapping.json`` — canonical, atom-index keyed
- ``reference/PSBMA_20mer/PSBMA20.map``          — Martini-style human view
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import MDAnalysis as mda
import numpy as np


REPO = Path(__file__).resolve().parents[1]
AA_GRO = REPO / "reference/PSBMA_20mer/PSBMA_20mer_no_water.gro"
OUT_JSON = REPO / "reference/PSBMA_20mer/PSBMA20_mapping.json"
OUT_MAP = REPO / "reference/PSBMA_20mer/PSBMA20.map"


# Static heavy-atom → bead-id (1-based) table, one monomer.
HEAVY_TO_BEAD: dict[str, int] = {
    # bead 1 (SC1, backbone)
    "C1": 1, "C2": 1, "C3": 1,
    # bead 2 (TN5a, ester C=O)
    "C4": 2, "O1": 2,
    # bead 3 (TP2a, ester O-CH2)
    "O2": 3, "C5": 3,
    # bead 4 (Q1+, ammonium head)
    "C6": 4, "N": 4, "C7": 4, "C8": 4,
    # bead 5 (SC1, propyl linker)
    "C9": 5, "C10": 5, "C11": 5,
    # bead 6 (Q1-, sulfonate)
    "S": 6, "O3": 6, "O4": 6, "O5": 6,
}


BEAD_META = [
    {"bead_id": 1, "bead_name": "C01", "bead_type": "SC1", "cg_mass": 54.0, "charge":  0, "comment": "backbone CCC"},
    {"bead_id": 2, "bead_name": "N01", "bead_type": "TN5a", "cg_mass": 36.0, "charge":  0, "comment": "ester C=O"},
    {"bead_id": 3, "bead_name": "P01", "bead_type": "TP2a", "cg_mass": 36.0, "charge":  0, "comment": "ester O-CH2"},
    {"bead_id": 4, "bead_name": "101", "bead_type": "Q1",   "cg_mass": 72.0, "charge":  1, "comment": "ammonium C-N+(C)(C)"},
    {"bead_id": 5, "bead_name": "C02", "bead_type": "SC1",  "cg_mass": 54.0, "charge":  0, "comment": "propyl linker CCC"},
    {"bead_id": 6, "bead_name": "102", "bead_type": "Q1",   "cg_mass": 72.0, "charge": -1, "comment": "sulfonate SO3-"},
]


def _atomic_mass(name: str) -> float:
    """Rough atomic mass from the leading letter of the AA atom name."""
    if name.upper().startswith("BR"):
        return 79.904
    return {"C": 12.011, "H": 1.008, "N": 14.007, "O": 15.999, "S": 32.06}[name[0].upper()]


def build_mapping(u: mda.Universe) -> dict:
    n_monomers = len(u.residues)
    # per-bead accumulator: bead_id → list[atom_index (1-based global)]
    beads_per_res: list[dict[int, list[int]]] = [
        defaultdict(list) for _ in range(n_monomers)
    ]

    # ---- pass 1: heavy atoms via static table
    for res_i, res in enumerate(u.residues):
        heavy_by_name: dict[str, tuple[int, np.ndarray]] = {}
        for atom in res.atoms:
            name = atom.name
            if name in HEAVY_TO_BEAD:
                bead_id = HEAVY_TO_BEAD[name]
                beads_per_res[res_i][bead_id].append(int(atom.index) + 1)
                heavy_by_name[name] = (int(atom.index) + 1, atom.position.copy())

        # ---- pass 2: within this residue, assign each remaining atom to
        # its nearest heavy atom's bead. This covers all H's plus the two
        # end-group extras (H6 on res 1, BR1 on res 20).
        heavy_positions = np.array([xyz for (_, xyz) in heavy_by_name.values()])
        heavy_beads = [HEAVY_TO_BEAD[name] for name in heavy_by_name]
        for atom in res.atoms:
            name = atom.name
            if name in HEAVY_TO_BEAD:
                continue
            d2 = np.sum((heavy_positions - atom.position) ** 2, axis=1)
            nearest = int(np.argmin(d2))
            bead_id = heavy_beads[nearest]
            beads_per_res[res_i][bead_id].append(int(atom.index) + 1)

    # ---- assemble the mapping JSON
    beads_out = []
    for res_i, buckets in enumerate(beads_per_res):
        resid = res_i + 1
        for meta in BEAD_META:
            bid = meta["bead_id"]
            atom_indices = sorted(buckets[bid])
            atom_names = [u.atoms[a - 1].name for a in atom_indices]
            aa_mass_sum = float(sum(_atomic_mass(n) for n in atom_names))
            heavy_indices = [
                a for (a, n) in zip(atom_indices, atom_names)
                if not (n.startswith("H") and not n.startswith("HG"))
            ]
            beads_out.append({
                "bead_id": (resid - 1) * 6 + bid,
                "bead_name": meta["bead_name"],
                "bead_type": meta["bead_type"],
                "cg_mass": meta["cg_mass"],
                "charge": meta["charge"],
                "aa_residue": resid,
                "aa_resname": u.residues[res_i].resname,
                "atom_indices": atom_indices,
                "atom_names": atom_names,
                "heavy_atom_indices": heavy_indices,
                "heavy_atom_count": len(heavy_indices),
                "aa_mass_sum": round(aa_mass_sum, 4),
                "comment": meta["comment"] + (
                    f" [+H6 end-group]" if resid == 1 and bid == 1 and "H6" in atom_names else
                    f" [+BR1 end-group]" if resid == n_monomers and bid == 1 and "BR1" in atom_names else
                    ""
                ),
            })

    return {
        "molecule": "PSBMA20",
        "weighting": "mass",
        "n_beads": len(beads_out),
        "n_monomers": n_monomers,
        "beads_per_monomer": 6,
        "monomer_smiles": "CC(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]",
        "provenance": {
            "aa_source": str(AA_GRO.relative_to(REPO)),
            "cg_template": "tests/fixtures/psbma/PSBA.itp",
            "heavy_atom_rule": "static table (SMILES↔AA name)",
            "hydrogen_rule": "spatial nearest-neighbour to a heavy atom in the reference .gro frame",
            "end_groups": "res 1 extra H6 and res 20 extra BR1 fold into their monomer's bead 1 via the nearest-heavy rule",
        },
        "beads": beads_out,
    }


def write_martini_map(mapping: dict, out_path: Path) -> None:
    lines = [
        "; PSBMA-20 AA→CG mapping (Martini-style mirror of PSBMA20_mapping.json)",
        f"; molecule = {mapping['molecule']}, weighting = {mapping['weighting']}",
        f"; SMILES   = {mapping['monomer_smiles']}",
        "",
    ]
    for bead in mapping["beads"]:
        lines.append(
            f"[ {bead['bead_name']}#{bead['bead_id']:03d}  {bead['bead_type']:<5s} "
            f"res {bead['aa_residue']:>2d} ]"
        )
        lines.append("  " + " ".join(bead["atom_names"]))
        lines.append("")
    out_path.write_text("\n".join(lines))


def main() -> None:
    u = mda.Universe(str(AA_GRO))
    n_atoms = len(u.atoms)
    n_res = len(u.residues)
    print(f"AA: {n_atoms} atoms, {n_res} residues")

    mapping = build_mapping(u)
    n_covered = sum(len(b["atom_indices"]) for b in mapping["beads"])
    print(f"beads: {len(mapping['beads'])}  atoms covered: {n_covered}/{n_atoms}")
    if n_covered != n_atoms:
        missing = set(range(1, n_atoms + 1))
        for b in mapping["beads"]:
            missing -= set(b["atom_indices"])
        raise SystemExit(f"ERROR: {len(missing)} AA atoms unassigned: {sorted(missing)[:20]}...")

    # Bead mass histogram
    mass_by_bead_id: dict[str, list[float]] = {}
    for b in mapping["beads"]:
        mass_by_bead_id.setdefault(b["bead_name"], []).append(b["aa_mass_sum"])
    print("mean AA mass per bead type:")
    for name, ms in mass_by_bead_id.items():
        print(f"  {name:<5s} n={len(ms):>2d}  mean={np.mean(ms):>6.2f}  min={np.min(ms):>6.2f}  max={np.max(ms):>6.2f}")

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(mapping, indent=2))
    print(f"wrote {OUT_JSON.relative_to(REPO)}")
    write_martini_map(mapping, OUT_MAP)
    print(f"wrote {OUT_MAP.relative_to(REPO)}")


if __name__ == "__main__":
    main()
