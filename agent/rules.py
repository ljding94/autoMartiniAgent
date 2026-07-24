"""Martini 3 rule checker (#3).

Audits a mapping against Martini's hard rules and returns structured violations.
Two severities:

  - **error** — a mapping that cannot be a valid Martini model: a bead with > 4
    heavy atoms (bigger than a regular bead), an empty bead, or a charged bead
    that is not a Q-type. These are the loop's *hard constraints* (pass
    ``hard_violations`` as ``run_loop``'s ``rule_check``).
  - **warn** — fixable by re-labelling, not by moving atoms: a bead whose heavy
    count does not match its size-class prefix (T=2 / S=3 / regular=4), a Q-type
    bead with no charge, an implausible total bead count, or a non-standard
    (agent-merged) composite type. ``relabel_to_fit`` re-types beads to match
    their heavy count — objective-neutral, since the bead type is only a label.

Bead type → size class from the Martini 3 prefix: a leading ``T`` = Tiny (2 heavy),
``S`` = Small (3), otherwise Regular (4). No regular Martini class starts with S or
T, so the prefix is unambiguous.

Not yet enforced: **functional-group integrity** (amide / ester / sulfonate /
ammonium not split across beads) — that needs SMARTS on the AA structure (RDKit)
and is left as a TODO; the sizing + charge + count checks below are structural and
need only the mapping.
"""

from __future__ import annotations

from dataclasses import dataclass

from agent import repair as _repair

_EXPECTED_HEAVY = {"T": 2, "S": 3, "R": 4}


def size_class(bead_type: str) -> str:
    """'T' (tiny), 'S' (small) or 'R' (regular) from the Martini bead-type prefix."""
    if bead_type and bead_type[0] == "T":
        return "T"
    if bead_type and bead_type[0] == "S":
        return "S"
    return "R"


def expected_heavy(bead_type: str) -> int:
    return _EXPECTED_HEAVY[size_class(bead_type)]


@dataclass(frozen=True)
class Violation:
    kind: str          # sizing | oversize | empty | charge_type | bead_count | composite_type
    severity: str      # "error" | "warn"
    message: str
    bead_id: int | None = None


def check_rules(state: _repair.MappingState) -> list[Violation]:
    """All Martini-rule violations in ``state`` (errors first, then warnings)."""
    v: list[Violation] = []
    total_heavy = 0
    for b in state.beads:
        bid, btype, n = b["bead_id"], b["bead_type"], b["heavy_atom_count"]
        charge = float(b.get("charge", 0) or 0)
        total_heavy += n

        if n == 0:
            v.append(Violation("empty", "error", f"bead {bid} has no heavy atoms", bid))
        elif n > 4:
            v.append(Violation("oversize", "error",
                               f"bead {bid} ({btype}) has {n} heavy atoms (> 4, larger than a regular bead)", bid))

        if "/" in btype:
            v.append(Violation("composite_type", "warn",
                               f"bead {bid} has a non-standard merged type {btype!r}; relabel to a Martini type", bid))
        elif 1 <= n <= 4:
            exp = expected_heavy(btype)
            if n != exp:
                v.append(Violation("sizing", "warn",
                                   f"bead {bid} ({btype}, {size_class(btype)}-class) has {n} heavy atoms, "
                                   f"expected {exp}; relabel size class", bid))

        is_q = "Q" in btype.upper()
        if charge != 0 and not is_q:
            v.append(Violation("charge_type", "error",
                               f"bead {bid} ({btype}) carries charge {charge:g} but is not a Q-type", bid))
        elif charge == 0 and is_q:
            v.append(Violation("charge_type", "warn",
                               f"bead {bid} ({btype}) is a Q-type but has no charge", bid))

    n_beads = state.n_beads
    plausible = total_heavy / 4.0
    if plausible > 0 and not (0.75 * plausible <= n_beads <= 1.25 * plausible):
        v.append(Violation("bead_count", "warn",
                           f"{n_beads} beads for {total_heavy} heavy atoms "
                           f"(~{plausible:.0f} expected at 4/bead; outside ±25%)"))

    v.sort(key=lambda x: 0 if x.severity == "error" else 1)
    return v


def hard_violations(state: _repair.MappingState) -> list[str]:
    """Error-severity messages only — the loop's hard constraint (rejects edits that
    make an unphysical mapping, while leaving relabel-fixable sizing to the agent)."""
    return [x.message for x in check_rules(state) if x.severity == "error"]


def relabel_to_fit(state: _repair.MappingState) -> _repair.MappingState:
    """Re-type each bead's size-class prefix to match its heavy-atom count (T=2 /
    S=3 / regular=4). Objective-neutral: the bead type is only a label, so this
    fixes size-class + composite-type warnings without touching any distribution.
    """
    new = _repair._clone(state)
    for b in new.beads:
        n = b["heavy_atom_count"]
        if n < 1 or n > 4:
            continue                       # can't be a valid size; leave for an error
        btype = b["bead_type"]
        base = btype.split("/")[0]         # collapse composite types onto their first class
        base = base[1:] if base[:1] in ("T", "S") else base
        prefix = {2: "T", 3: "S", 4: ""}[n]
        b["bead_type"] = prefix + base
    return new


def format_report(state: _repair.MappingState) -> str:
    vs = check_rules(state)
    errs = [x for x in vs if x.severity == "error"]
    warns = [x for x in vs if x.severity == "warn"]
    lines = [f"Martini rule check: {len(errs)} error(s), {len(warns)} warning(s)"]
    for x in vs:
        lines.append(f"  [{x.severity}] {x.kind}: {x.message}")
    if not vs:
        lines.append("  (clean)")
    return "\n".join(lines)
