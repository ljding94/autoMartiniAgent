"""Mapping mutation verbs for the agent repair loop (#4).

The agent's optimization variable is the AA-atom -> CG-bead assignment in the
mapping JSON. This module provides the geometry-changing verbs the agent uses to
edit that assignment, plus the IO needed to feed an edited mapping back through
the projector + scorer:

  - ``reassign_atoms`` — move AA atoms between two bonded beads (fixed bead count)
  - ``merge_beads``    — combine beads into one (bead count shrinks)
  - ``split_bead``     — divide a bead's atoms into two (bead count grows)

Every verb takes a ``MappingState`` and returns a *new* one (inputs are never
mutated), re-deriving each touched bead's masses / heavy-atom counts and keeping
the bead-bond graph consistent so a valid parameter-free ``.itp`` can be written.

Why parameter-free? The repair objective is Gaussianity of the AA-projected
distributions (``1 - R^2``), which needs only bead *connectivity* to know which
bond / angle distributions to measure — not force-field targets. So for beads the
agent invents (via merge/split) we carry a bead-bond graph and emit an ``.itp``
with real ``(b0, kb)`` group tags preserved where they exist and unique sentinel
tags for new bonds (so distinct bond families never collapse into one group). See
``agent/evaluate.py`` for the project+score+objective glue that consumes this.
"""

from __future__ import annotations

import copy
import json
import re
from dataclasses import dataclass, field
from pathlib import Path


# ---------- atomic masses (self-contained; matches tests/test_project.py) ----------

_ELEMENT_MASS = {
    "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.9994, "S": 32.06,
    "P": 30.9738, "BR": 79.904, "CL": 35.45, "F": 18.998, "NA": 22.99,
}
_TWO_LETTER = {"BR", "CL", "NA"}


def _element_of(name: str) -> str:
    s = "".join(ch for ch in str(name) if ch.isalpha()).upper()
    if len(s) >= 2 and s[:2] in _TWO_LETTER:
        return s[:2]
    return s[:1] if s else "C"


def _mass_from_name(name: str) -> float:
    return _ELEMENT_MASS.get(_element_of(name), 12.011)


# ---------- bead-bond graph ----------


@dataclass(frozen=True)
class BeadBond:
    """A bond between two 1-based bead indices, carrying its grouping tag.

    ``b0`` / ``kb`` are used purely as a group key by the scorer (bonds sharing a
    ``(b0, kb)`` are one distribution). For real force-field bonds they are the
    itp targets; for bonds the agent creates they are unique sentinels.
    """
    i: int
    j: int
    b0: float | None
    kb: float | None

    def key(self) -> tuple[int, int]:
        return (self.i, self.j) if self.i <= self.j else (self.j, self.i)


def _parse_itp_bonds(itp_path: str | Path) -> list[BeadBond]:
    """Minimal GROMACS ``[ bonds ]`` reader (no mdtraj/scipy dependency)."""
    bonds: list[BeadBond] = []
    section: str | None = None
    for raw in Path(itp_path).read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        m = re.match(r"\[\s*(\w+)\s*\]", line)
        if m:
            section = m.group(1).lower()
            continue
        if section != "bonds":
            continue
        toks = line.split()
        if len(toks) < 2:
            continue
        try:
            i, j = int(toks[0]), int(toks[1])
        except ValueError:
            continue
        b0 = float(toks[3]) if len(toks) > 3 else None
        kb = float(toks[4]) if len(toks) > 4 else None
        bonds.append(BeadBond(i=i, j=j, b0=b0, kb=kb))
    return bonds


# ---------- state ----------


@dataclass
class MappingState:
    """Working state: the mapping JSON plus a bead-bond graph and atom lookups.

    ``atom_name`` / ``heavy_atoms`` / ``atom_masses`` are global (all AA atoms) so
    beads can be re-derived after any edit. ``bead_id`` is kept equal to the bead's
    1-based position, so bond indices, itp ``[atoms]`` nr, and the projected CG
    universe order all agree.
    """
    mapping: dict
    bonds: list[BeadBond]
    atom_name: dict[int, str]
    heavy_atoms: frozenset[int]
    atom_masses: dict[int, float]
    all_atoms: frozenset[int]
    _sentinel_counter: int = field(default=0)

    @property
    def beads(self) -> list[dict]:
        return self.mapping["beads"]

    @property
    def n_beads(self) -> int:
        return len(self.mapping["beads"])

    def bead(self, bead_id: int) -> dict:
        beads = self.beads
        if not 1 <= bead_id <= len(beads):
            raise ValueError(f"bead_id {bead_id} out of range (1..{len(beads)})")
        b = beads[bead_id - 1]
        assert b["bead_id"] == bead_id, "bead_id/position invariant broken"
        return b

    def bead_of_atom(self, atom_index: int) -> int:
        for b in self.beads:
            if atom_index in b["atom_indices"]:
                return b["bead_id"]
        raise ValueError(f"atom {atom_index} not assigned to any bead")

    def neighbors(self, bead_id: int) -> set[int]:
        out: set[int] = set()
        for bd in self.bonds:
            if bd.i == bead_id:
                out.add(bd.j)
            elif bd.j == bead_id:
                out.add(bd.i)
        return out

    def _next_sentinel(self) -> float:
        """A unique negative ``b0`` group tag for an agent-created bond."""
        self._sentinel_counter += 1
        return -float(self._sentinel_counter)


def load_state(
    mapping_path: str | Path,
    itp_path: str | Path,
    atom_masses: dict[int, float] | None = None,
) -> MappingState:
    """Load a mapping JSON + its itp into a MappingState.

    ``atom_masses`` (1-based AA index -> mass) overrides the name-derived masses;
    pass the projector's actual AA masses here so re-derived ``aa_mass_sum`` values
    match what ``agent.project`` validates against.
    """
    mapping = json.loads(Path(mapping_path).read_text())
    bonds = _parse_itp_bonds(itp_path)

    atom_name: dict[int, str] = {}
    heavy: set[int] = set()
    all_atoms: set[int] = set()
    for bd in mapping["beads"]:
        for idx, nm in zip(bd["atom_indices"], bd["atom_names"]):
            atom_name[idx] = nm
            all_atoms.add(idx)
        heavy.update(bd.get("heavy_atom_indices", []))

    if atom_masses is None:
        atom_masses = {idx: _mass_from_name(nm) for idx, nm in atom_name.items()}

    state = MappingState(
        mapping=mapping,
        bonds=bonds,
        atom_name=atom_name,
        heavy_atoms=frozenset(heavy),
        atom_masses=atom_masses,
        all_atoms=frozenset(all_atoms),
    )
    _reindex(state)  # normalize bead_id == position
    validate_state(state)
    return state


# ---------- internals ----------


def _refresh_bead(state: MappingState, bead: dict) -> None:
    """Re-derive a bead's atom_names / heavy counts / masses from atom_indices."""
    idxs = sorted(bead["atom_indices"])
    bead["atom_indices"] = idxs
    bead["atom_names"] = [state.atom_name[i] for i in idxs]
    heavy = [i for i in idxs if i in state.heavy_atoms]
    bead["heavy_atom_indices"] = heavy
    bead["heavy_atom_count"] = len(heavy)
    m = float(sum(state.atom_masses[i] for i in idxs))
    bead["aa_mass_sum"] = round(m, 3)
    bead["cg_mass"] = round(m, 3)


def _reindex(state: MappingState) -> None:
    """Renumber bead_id to 1..N by list position and remap bond indices.

    Bonds reference the *current* ``bead_id`` values before this call; they are
    rewritten to the new positional ids so the invariant bead_id == position holds.
    """
    old_ids = [b["bead_id"] for b in state.beads]
    remap = {old: pos for pos, old in enumerate(old_ids, start=1)}
    for pos, b in enumerate(state.beads, start=1):
        b["bead_id"] = pos
    new_bonds: list[BeadBond] = []
    seen: set[tuple[int, int]] = set()
    for bd in state.bonds:
        if bd.i not in remap or bd.j not in remap:
            continue  # bond to a removed bead
        i, j = remap[bd.i], remap[bd.j]
        if i == j:
            continue  # self-loop (both endpoints merged)
        nb = BeadBond(i=i, j=j, b0=bd.b0, kb=bd.kb)
        if nb.key() in seen:
            continue
        seen.add(nb.key())
        new_bonds.append(nb)
    state.bonds = new_bonds
    state.mapping["n_beads"] = len(state.beads)


def validate_state(state: MappingState) -> None:
    """Assert the atom partition is intact: every AA atom in exactly one bead."""
    seen: dict[int, int] = {}
    for b in state.beads:
        if not b["atom_indices"]:
            raise ValueError(f"bead {b['bead_id']} ({b.get('bead_name')}) is empty")
        for i in b["atom_indices"]:
            if i in seen:
                raise ValueError(
                    f"atom {i} in both bead {seen[i]} and bead {b['bead_id']}"
                )
            seen[i] = b["bead_id"]
    got = frozenset(seen)
    if got != state.all_atoms:
        missing = sorted(state.all_atoms - got)
        extra = sorted(got - state.all_atoms)
        raise ValueError(f"atom set changed: missing={missing} extra={extra}")


# ---------- verbs ----------


def reassign_atoms(
    state: MappingState,
    atom_indices: list[int],
    to_bead_id: int,
    *,
    require_adjacent: bool = True,
) -> MappingState:
    """Move ``atom_indices`` into bead ``to_bead_id`` (fixed bead count).

    The moved atoms may come from any bead(s); by default the target must be a
    bonded neighbour of each source bead (``require_adjacent``) so moves stay local
    and physically meaningful (the classic "shift a boundary atom to the lighter
    neighbour" repair). Bead count and the bond graph are unchanged.
    """
    new = _clone(state)
    target = new.bead(to_bead_id)
    move = set(atom_indices)
    if not move:
        raise ValueError("no atoms to reassign")

    sources: set[int] = set()
    for a in move:
        src = new.bead_of_atom(a)
        if src != to_bead_id:
            sources.add(src)
    if require_adjacent:
        adj = new.neighbors(to_bead_id)
        for s in sources:
            if s not in adj:
                raise ValueError(
                    f"bead {s} is not bonded to target bead {to_bead_id}; "
                    f"pass require_adjacent=False to override"
                )

    touched = {to_bead_id} | sources
    for bd in new.beads:
        if bd["bead_id"] in sources:
            bd["atom_indices"] = [i for i in bd["atom_indices"] if i not in move]
    tgt_atoms = set(target["atom_indices"]) | move
    target["atom_indices"] = sorted(tgt_atoms)
    for bid in touched:
        _refresh_bead(new, new.bead(bid))
    validate_state(new)
    return new


def merge_beads(state: MappingState, bead_ids: list[int]) -> MappingState:
    """Merge ``bead_ids`` into a single bead (bead count shrinks by len-1).

    Atoms are unioned into the lowest-id bead; bonds to any absorbed bead are
    redirected to it (self-loops and duplicates dropped). The merged bead's type /
    name become composites when the inputs disagree. Beads need not be pairwise
    bonded, but a warning-worthy non-adjacent merge is the caller's call.
    """
    ids = sorted(set(bead_ids))
    if len(ids) < 2:
        raise ValueError("merge needs >= 2 distinct beads")
    new = _clone(state)
    for bid in ids:
        new.bead(bid)  # range check

    keep_id = ids[0]
    absorb = set(ids[1:])
    keep = new.bead(keep_id)

    atoms = set(keep["atom_indices"])
    types: list[str] = [keep["bead_type"]]
    names: list[str] = [str(keep.get("bead_name", keep_id))]
    for bid in ids[1:]:
        b = new.bead(bid)
        atoms.update(b["atom_indices"])
        types.append(b["bead_type"])
        names.append(str(b.get("bead_name", bid)))
    keep["atom_indices"] = sorted(atoms)

    uniq_types = list(dict.fromkeys(types))
    keep["bead_type"] = uniq_types[0] if len(uniq_types) == 1 else "/".join(uniq_types)
    keep["bead_name"] = names[0] if len(set(names)) == 1 else "+".join(names)
    keep["comment"] = "merged " + "+".join(names)
    keep["charge"] = int(round(sum(new.bead(b)["charge"] for b in ids)))

    # redirect bonds from absorbed beads to keep_id, then drop absorbed beads
    redirected: list[BeadBond] = []
    for bd in new.bonds:
        i = keep_id if bd.i in absorb else bd.i
        j = keep_id if bd.j in absorb else bd.j
        redirected.append(BeadBond(i=i, j=j, b0=bd.b0, kb=bd.kb))
    new.bonds = redirected
    new.mapping["beads"] = [b for b in new.beads if b["bead_id"] not in absorb]

    _refresh_bead(new, keep)
    _reindex(new)
    new.mapping.pop("beads_per_monomer", None)  # bead/monomer regularity broken
    validate_state(new)
    return new


def split_bead(
    state: MappingState,
    bead_id: int,
    group_a: list[int],
    *,
    ref_positions: dict[int, tuple[float, float, float]] | None = None,
) -> MappingState:
    """Split bead ``bead_id`` into two: ``group_a`` and the remaining atoms.

    A new bead is inserted immediately after the original. External bonds of the
    original bead are re-attached to whichever child is geometrically closer to the
    neighbour (needs ``ref_positions``: 1-based AA index -> xyz, e.g. from the CG
    reference .gro frame); without positions they all stay on child A and a warning
    situation is left to the caller. A fresh A–B bond (unique sentinel tag) is
    always added.
    """
    new = _clone(state)
    orig = new.bead(bead_id)
    a_atoms = sorted(set(group_a))
    if not a_atoms:
        raise ValueError("group_a is empty")
    orig_atoms = set(orig["atom_indices"])
    if not set(a_atoms) <= orig_atoms:
        raise ValueError(f"group_a atoms {sorted(set(a_atoms) - orig_atoms)} not in bead {bead_id}")
    b_atoms = sorted(orig_atoms - set(a_atoms))
    if not b_atoms:
        raise ValueError("split leaves the other child empty (group_a is the whole bead)")

    new_id = max(b["bead_id"] for b in new.beads) + 1  # temporary unique id
    child_b = copy.deepcopy(orig)
    child_b["bead_id"] = new_id
    child_b["bead_name"] = f"{orig.get('bead_name', bead_id)}b"
    child_b["comment"] = f"split-b of {orig.get('bead_name', bead_id)}"
    child_b["atom_indices"] = b_atoms
    orig["bead_name"] = f"{orig.get('bead_name', bead_id)}a"
    orig["comment"] = f"split-a of {orig.get('comment', bead_id)}"
    orig["atom_indices"] = a_atoms

    _refresh_bead(new, orig)
    _refresh_bead(new, child_b)

    # decide, per external neighbour, whether it bonds to child A (orig) or B
    def _com(atoms: list[int]) -> tuple[float, float, float] | None:
        if not ref_positions:
            return None
        pts = [ref_positions[i] for i in atoms if i in ref_positions]
        if not pts:
            return None
        n = len(pts)
        return tuple(sum(p[k] for p in pts) / n for k in range(3))  # type: ignore

    com_a = _com(a_atoms) if ref_positions else None
    com_b = _com(b_atoms) if ref_positions else None

    def _closer_to_b(neighbor_id: int) -> bool:
        if com_a is None or com_b is None:
            return False  # fallback: keep external bonds on child A
        nb = new.bead(neighbor_id)
        cn = _com(nb["atom_indices"])
        if cn is None:
            return False
        da = sum((com_a[k] - cn[k]) ** 2 for k in range(3))
        db = sum((com_b[k] - cn[k]) ** 2 for k in range(3))
        return db < da

    rebuilt: list[BeadBond] = []
    for bd in new.bonds:
        if bd.i == bead_id or bd.j == bead_id:
            other = bd.j if bd.i == bead_id else bd.i
            endpoint = new_id if _closer_to_b(other) else bead_id
            i, j = (endpoint, other) if bd.i == bead_id else (other, endpoint)
            rebuilt.append(BeadBond(i=i, j=j, b0=bd.b0, kb=bd.kb))
        else:
            rebuilt.append(bd)
    rebuilt.append(BeadBond(i=bead_id, j=new_id, b0=new._next_sentinel(), kb=None))
    new.bonds = rebuilt

    # insert child B right after the original in list order
    pos = next(k for k, b in enumerate(new.beads) if b["bead_id"] == bead_id)
    new.mapping["beads"].insert(pos + 1, child_b)

    _reindex(new)
    new.mapping.pop("beads_per_monomer", None)  # bead/monomer regularity broken
    validate_state(new)
    return new


def _clone(state: MappingState) -> MappingState:
    return MappingState(
        mapping=copy.deepcopy(state.mapping),
        bonds=list(state.bonds),
        atom_name=state.atom_name,
        heavy_atoms=state.heavy_atoms,
        atom_masses=state.atom_masses,
        all_atoms=state.all_atoms,
        _sentinel_counter=state._sentinel_counter,
    )


# ---------- output ----------


def write_mapping(state: MappingState, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(state.mapping, indent=2))
    return path


def write_itp(state: MappingState, path: str | Path, mol_name: str | None = None) -> Path:
    """Write a parameter-free ``.itp`` (``[atoms]`` + ``[bonds]``, no ``[angles]``).

    Real ``(b0, kb)`` bond tags are preserved so existing bond families group as
    before; agent-created bonds carry unique sentinel tags. Angles are intentionally
    omitted — the scorer derives every consecutive-bond angle from the bonds when
    run with ``all_bonded_angles=True``.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    mol = mol_name or state.mapping.get("molecule", "MOL")
    lines: list[str] = []
    lines.append(f"; {mol} CG topology — parameter-free (agent-edited mapping)")
    lines.append("; bond b0/kb are group tags only; negative b0 = agent-created bond")
    lines.append("[ moleculetype ]")
    lines.append(f"{mol}   1")
    lines.append("")
    lines.append("[ atoms ]")
    lines.append(";  nr  type  resnr  residue  atom  cgnr  charge  mass")
    for b in state.beads:
        nr = b["bead_id"]
        resnr = b.get("aa_residue", 1)
        resname = b.get("aa_resname", mol)[:5]
        aname = str(b.get("bead_name", nr))
        charge = float(b.get("charge", 0) or 0)
        mass = float(b.get("cg_mass", 0.0))
        lines.append(
            f"{nr:6d} {b['bead_type']:>8} {resnr:6d}  {resname:<5} {aname:>5} "
            f"{nr:6d} {charge:8.3f} {mass:8.2f}"
        )
    lines.append("")
    lines.append("[ bonds ]")
    lines.append(";   i    j  func    b0        kb")
    for bd in sorted(state.bonds, key=lambda x: x.key()):
        i, j = bd.key()
        if bd.b0 is None and bd.kb is None:
            lines.append(f"{i:6d}{j:6d}    1")
        else:
            b0 = bd.b0 if bd.b0 is not None else 0.0
            kb = bd.kb if bd.kb is not None else 0.0
            lines.append(f"{i:6d}{j:6d}    1  {b0:8.4f}  {kb:10.2f}")
    lines.append("")
    path.write_text("\n".join(lines))
    return path


def load_ref_positions(struct_path: str | Path) -> dict[int, tuple[float, float, float]]:
    """1-based AA atom index -> (x, y, z) in nm from a GROMACS ``.gro`` frame.

    Used by ``split_bead`` to route external bonds to the nearer child.
    """
    text = Path(struct_path).read_text().splitlines()
    n = int(text[1].strip())
    out: dict[int, tuple[float, float, float]] = {}
    for k, line in enumerate(text[2 : 2 + n], start=1):
        x = float(line[20:28]); y = float(line[28:36]); z = float(line[36:44])
        out[k] = (x, y, z)
    return out
