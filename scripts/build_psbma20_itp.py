"""Tile the single-monomer PSBA.itp into a full 20-mer PSBMA20.itp.

Approach
--------
1. Read the atoms/bonds/angles from the auto-martiniM3 monomer template
   (``tests/fixtures/psbma/PSBA.itp``).
2. Emit 20 copies of the atom block, offsetting bead index and cgnr by
   ``6*(n-1)`` per monomer.
3. Emit 20 copies of the intra-monomer bond and angle blocks with the same
   offset applied to all indices.
4. Add 19 inter-monomer backbone bonds linking bead 1 (SC1 backbone) of
   monomer *n* to bead 1 of monomer *n+1*, with the parameters chosen for
   the Martini-3 methacrylate backbone (b0 = 0.27 nm, kb = 7500).

Output: ``reference/PSBMA_20mer/PSBMA20.itp``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
MONOMER_ITP = REPO / "tests/fixtures/psbma/PSBA.itp"
OUT_ITP = REPO / "reference/PSBMA_20mer/PSBMA20.itp"

N_MONOMERS = 20
BEADS_PER_MONOMER = 6

# Inter-monomer backbone bond (SC1 of monomer N → SC1 of monomer N+1)
INTER_B0 = 0.27
INTER_KB = 7500.0
INTER_FUNC = 1


@dataclass
class Atom:
    idx: int
    type: str
    resnr: int
    resname: str
    name: str
    cgnr: int
    charge: float
    mass: float


@dataclass
class Bond:
    i: int
    j: int
    func: int
    b0: float
    kb: float


@dataclass
class Angle:
    i: int
    j: int
    k: int
    func: int
    theta0: float
    ka: float


def parse_monomer(itp_path: Path) -> tuple[list[Atom], list[Bond], list[Angle]]:
    atoms: list[Atom] = []
    bonds: list[Bond] = []
    angles: list[Angle] = []
    section: str | None = None
    for raw in itp_path.read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        m = re.match(r"\[\s*(\w+)\s*\]", line)
        if m:
            section = m.group(1).lower()
            continue
        toks = line.split()
        try:
            if section == "atoms":
                atoms.append(Atom(
                    idx=int(toks[0]), type=toks[1], resnr=int(toks[2]),
                    resname=toks[3], name=toks[4], cgnr=int(toks[5]),
                    charge=float(toks[6]), mass=float(toks[7]),
                ))
            elif section == "bonds":
                bonds.append(Bond(
                    i=int(toks[0]), j=int(toks[1]), func=int(toks[2]),
                    b0=float(toks[3]), kb=float(toks[4]),
                ))
            elif section == "angles":
                angles.append(Angle(
                    i=int(toks[0]), j=int(toks[1]), k=int(toks[2]),
                    func=int(toks[3]), theta0=float(toks[4]), ka=float(toks[5]),
                ))
        except (ValueError, IndexError):
            continue
    return atoms, bonds, angles


def tile_topology(
    atoms: list[Atom], bonds: list[Bond], angles: list[Angle],
    n_monomers: int,
) -> tuple[list[Atom], list[Bond], list[Angle]]:
    per_mono = len(atoms)
    tiled_atoms: list[Atom] = []
    tiled_bonds: list[Bond] = []
    tiled_angles: list[Angle] = []
    for m in range(n_monomers):
        off = m * per_mono
        for a in atoms:
            tiled_atoms.append(Atom(
                idx=a.idx + off,
                type=a.type,
                resnr=m + 1,
                resname="PSBMA",
                name=a.name,
                cgnr=a.cgnr + off,
                charge=a.charge,
                mass=a.mass,
            ))
        for b in bonds:
            tiled_bonds.append(Bond(
                i=b.i + off, j=b.j + off, func=b.func, b0=b.b0, kb=b.kb,
            ))
        for ang in angles:
            tiled_angles.append(Angle(
                i=ang.i + off, j=ang.j + off, k=ang.k + off,
                func=ang.func, theta0=ang.theta0, ka=ang.ka,
            ))
    # Add inter-monomer backbone bonds (bead 1 of monomer m → bead 1 of m+1)
    for m in range(n_monomers - 1):
        i = m * per_mono + 1
        j = (m + 1) * per_mono + 1
        tiled_bonds.append(Bond(
            i=i, j=j, func=INTER_FUNC, b0=INTER_B0, kb=INTER_KB,
        ))
    return tiled_atoms, tiled_bonds, tiled_angles


def write_itp(
    atoms: list[Atom], bonds: list[Bond], angles: list[Angle],
    out_path: Path,
) -> None:
    lines: list[str] = []
    lines.append("; PSBMA-20 CG topology")
    lines.append("; Tiled from tests/fixtures/psbma/PSBA.itp (single-monomer auto-martiniM3 output)")
    lines.append("; Inter-monomer backbone bond: SC1(m)-SC1(m+1) at b0={} nm, kb={} (Martini 3 methacrylate default)".format(INTER_B0, INTER_KB))
    lines.append(";")
    lines.append("[ moleculetype ]")
    lines.append("PSBMA20   1")
    lines.append("")
    lines.append("[ atoms ]")
    lines.append(";  nr  type  resnr  residue  atom  cgnr  charge  mass")
    for a in atoms:
        lines.append(
            f"  {a.idx:>4d}  {a.type:<5s}  {a.resnr:>3d}  {a.resname:<5s}  "
            f"{a.name:<4s}  {a.cgnr:>4d}  {a.charge:>7.3f}  {a.mass:>7.2f}"
        )
    lines.append("")
    lines.append("[ bonds ]")
    lines.append(";   i    j  func    b0        kb        ; comment")
    n_intra = len([b for b in bonds if b.b0 != INTER_B0 or b.kb != INTER_KB])
    n_intra_per_mono = n_intra // N_MONOMERS
    for k, b in enumerate(bonds):
        is_inter = (k >= N_MONOMERS * n_intra_per_mono)
        comment = "; inter-monomer backbone" if is_inter else "; intra-monomer"
        lines.append(
            f"  {b.i:>4d} {b.j:>4d}  {b.func:>3d}  {b.b0:>8.4f}  {b.kb:>10.2f}  {comment}"
        )
    lines.append("")
    lines.append("[ angles ]")
    lines.append(";   i    j    k  func  theta0    ka")
    for ang in angles:
        lines.append(
            f"  {ang.i:>4d} {ang.j:>4d} {ang.k:>4d}  {ang.func:>3d}  {ang.theta0:>7.2f}  {ang.ka:>8.2f}"
        )
    lines.append("")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines))


def main() -> None:
    atoms, bonds, angles = parse_monomer(MONOMER_ITP)
    print(f"monomer template: {len(atoms)} atoms, {len(bonds)} bonds, {len(angles)} angles")
    tiled_atoms, tiled_bonds, tiled_angles = tile_topology(
        atoms, bonds, angles, N_MONOMERS
    )
    n_inter = N_MONOMERS - 1
    print(f"tiled  20-mer   : {len(tiled_atoms)} atoms, {len(tiled_bonds)} bonds "
          f"({len(tiled_bonds) - n_inter} intra + {n_inter} inter), {len(tiled_angles)} angles")
    write_itp(tiled_atoms, tiled_bonds, tiled_angles, OUT_ITP)
    print(f"wrote {OUT_ITP.relative_to(REPO)}")


if __name__ == "__main__":
    main()
