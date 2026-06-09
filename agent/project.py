"""AA→CG trajectory projector.

Reads an atomistic trajectory + a mapping JSON in the autoMartiniAgent
schema, projects each AA frame through mass-weighted COMs into a
coarse-grained trajectory, and writes a CG ``.xtc`` plus a single-frame
CG ``.gro`` for downstream analysis tools.

Mapping JSON contract (see ``reference/polyply_PEO20/PEO20_mapping.json``):
  - ``beads``: list of {``bead_id``, ``bead_name``, ``bead_type``,
    ``cg_mass``, ``aa_resname``, ``atom_indices`` (1-based GROMACS
    convention), ...}.
  - ``molecule``: short name used for residue label in the CG ``.gro``.

CLI::

  python -m agent.project \\
    --aa-top reference/gromacs/equil.tpr \\
    --aa-traj reference/gromacs/step5_200_center.xtc \\
    --mapping reference/polyply_PEO20/PEO20_mapping.json \\
    --out-dir derived/PEO20

Output filenames default to ``<molecule>_cg.xtc`` and ``<molecule>_cg.gro``
inside ``--out-dir`` (the dir is created if absent).
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import MDAnalysis as mda
import numpy as np


@dataclass(frozen=True)
class ProjectionResult:
    n_frames: int
    n_beads: int
    out_traj: Path
    out_struct: Path


def load_mapping(path: str | Path) -> dict:
    path = Path(path)
    with path.open() as f:
        return json.load(f)


def _build_bead_groups(universe: mda.Universe, mapping: dict) -> list[mda.AtomGroup]:
    """Resolve each bead's 1-based atom_indices to an MDAnalysis AtomGroup.

    GROMACS indices are 1-based; MDAnalysis is 0-based — we shift here.
    """
    groups: list[mda.AtomGroup] = []
    n_atoms = len(universe.atoms)
    for bead in mapping["beads"]:
        idxs = np.asarray(bead["atom_indices"], dtype=int) - 1
        if idxs.min() < 0 or idxs.max() >= n_atoms:
            raise ValueError(
                f"bead {bead['bead_id']} atom_indices out of range "
                f"for universe with {n_atoms} atoms "
                f"(got min={idxs.min() + 1}, max={idxs.max() + 1})"
            )
        groups.append(universe.atoms[idxs])
    return groups


def _validate_mass_sums(
    groups: list[mda.AtomGroup], mapping: dict, tol: float = 0.02
) -> None:
    """Check each bead's AA mass sum agrees with the mapping's recorded value.

    Catches index off-by-ones (the most common projection bug) before any
    expensive trajectory IO.
    """
    for bg, bead in zip(groups, mapping["beads"]):
        if "aa_mass_sum" not in bead:
            continue
        got = float(bg.total_mass())
        expected = float(bead["aa_mass_sum"])
        if abs(got - expected) > tol:
            raise ValueError(
                f"bead {bead['bead_id']} mass mismatch: "
                f"AtomGroup={got:.4f} g/mol, "
                f"mapping={expected:.4f} g/mol (tol {tol})"
            )


def _build_cg_universe(mapping: dict) -> mda.Universe:
    """Construct a CG Universe sized to the bead list, ready to receive frames.

    Topology attributes attached: names, types, resnames, resids, masses,
    segids. One bead per residue per Martini convention.
    """
    beads = mapping["beads"]
    n = len(beads)
    cg_u = mda.Universe.empty(
        n_atoms=n,
        n_residues=n,
        n_segments=1,
        atom_resindex=np.arange(n, dtype=int),
        residue_segindex=np.zeros(n, dtype=int),
        trajectory=True,
    )
    cg_u.add_TopologyAttr("names", [b["bead_name"] for b in beads])
    cg_u.add_TopologyAttr("types", [b["bead_type"] for b in beads])
    cg_u.add_TopologyAttr(
        "resnames", [mapping.get("molecule", "MOL")[:4].upper()] * n
    )
    cg_u.add_TopologyAttr("resids", [b["bead_id"] for b in beads])
    cg_u.add_TopologyAttr("masses", [b["cg_mass"] for b in beads])
    cg_u.add_TopologyAttr("segids", ["A"])
    return cg_u


def project_trajectory(
    aa_top: str | Path,
    aa_traj: str | Path,
    mapping: dict | str | Path,
    out_traj: str | Path,
    out_struct: str | Path,
) -> ProjectionResult:
    """Project an AA trajectory through ``mapping`` into a CG trajectory.

    Parameters
    ----------
    aa_top : path
        AA topology (``.tpr``, ``.gro``, ``.pdb``, ...). Must include masses
        for correct mass-weighted COM.
    aa_traj : path
        AA trajectory (``.xtc``, ``.trr``, ``.dcd``, ...).
    mapping : dict | path
        Mapping in the autoMartiniAgent schema. If a path, it is read first.
    out_traj : path
        Destination CG ``.xtc`` (extension drives the writer).
    out_struct : path
        Destination single-frame CG ``.gro`` (written from frame 0).
    """
    mapping_dict = (
        mapping if isinstance(mapping, dict) else load_mapping(mapping)
    )
    out_traj = Path(out_traj)
    out_struct = Path(out_struct)
    out_traj.parent.mkdir(parents=True, exist_ok=True)
    out_struct.parent.mkdir(parents=True, exist_ok=True)

    u = mda.Universe(str(aa_top), str(aa_traj))
    bead_groups = _build_bead_groups(u, mapping_dict)
    _validate_mass_sums(bead_groups, mapping_dict)

    cg_u = _build_cg_universe(mapping_dict)
    n_beads = len(bead_groups)

    n_frames = 0
    with mda.Writer(str(out_traj), n_atoms=n_beads) as writer:
        for ts in u.trajectory:
            positions = np.vstack(
                [bg.center_of_mass() for bg in bead_groups]
            )
            cg_u.atoms.positions = positions
            if ts.dimensions is not None:
                cg_u.dimensions = ts.dimensions
            writer.write(cg_u.atoms)
            if n_frames == 0:
                cg_u.atoms.write(str(out_struct))
            n_frames += 1

    return ProjectionResult(
        n_frames=n_frames,
        n_beads=n_beads,
        out_traj=out_traj,
        out_struct=out_struct,
    )


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Project an AA trajectory through a Martini mapping JSON"
    )
    p.add_argument("--aa-top", required=True, help="AA topology (.tpr/.gro/.pdb)")
    p.add_argument("--aa-traj", required=True, help="AA trajectory (.xtc/.trr/...)")
    p.add_argument("--mapping", required=True, help="mapping JSON path")
    p.add_argument(
        "--out-dir",
        default=None,
        help="output directory (default: derived/<molecule>)",
    )
    p.add_argument(
        "--out-traj",
        default=None,
        help="override CG trajectory filename (default: <molecule>_cg.xtc)",
    )
    p.add_argument(
        "--out-struct",
        default=None,
        help="override CG structure filename (default: <molecule>_cg.gro)",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    mapping = load_mapping(args.mapping)
    molecule = mapping.get("molecule", "MOL")
    out_dir = Path(args.out_dir or f"derived/{molecule}")
    out_traj = Path(args.out_traj) if args.out_traj else out_dir / f"{molecule}_cg.xtc"
    out_struct = (
        Path(args.out_struct) if args.out_struct else out_dir / f"{molecule}_cg.gro"
    )

    result = project_trajectory(
        aa_top=args.aa_top,
        aa_traj=args.aa_traj,
        mapping=mapping,
        out_traj=out_traj,
        out_struct=out_struct,
    )
    print(
        f"projected {result.n_frames} frame(s) of {Path(args.aa_traj).name} "
        f"through {len(mapping['beads'])}-bead mapping "
        f"({mapping.get('molecule', 'MOL')})"
    )
    print(f"  CG trajectory : {result.out_traj}")
    print(f"  CG structure  : {result.out_struct}")


if __name__ == "__main__":
    main()
