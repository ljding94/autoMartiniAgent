"""Tests for the repair-loop controller (``agent.loop``).

Trajectory-free and LLM-free: the controller is exercised with the real PSBMA
mapping fixture, a **fake** ``evaluate_fn`` (a pure function of the mapping, no
projection), and a ``ScriptedPolicy``. This isolates the loop mechanics
(apply → rule-check → evaluate → keep-best-valid → terminate) from the science.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from agent import loop, repair

REPO = Path(__file__).resolve().parents[1]
MAP = REPO / "reference/PSBMA_20mer/PSBMA20_mapping.json"
ITP = REPO / "reference/PSBMA_20mer/PSBMA20.itp"


@pytest.fixture(scope="module")
def state() -> repair.MappingState:
    return repair.load_state(MAP, ITP)


def _fake_eval(objective_fn):
    """Build an evaluate_fn returning .objective + an empty .breakdown.terms."""
    def _eval(st):
        return SimpleNamespace(objective=float(objective_fn(st)),
                               breakdown=SimpleNamespace(terms=[]))
    return _eval


# distance from the all-4-heavy ideal; a clean, mapping-only synthetic objective
def _sq_dev(st):
    return sum((b["heavy_atom_count"] - 4) ** 2 for b in st.beads)


# ---------- apply_ops (tiled) ----------


def test_apply_ops_tiled_reassign_all_monomers(state):
    # move C9 (+2H) from the propyl bead (role 5) into the ammonium bead (role 4),
    # tiled across all 20 monomers
    out = loop.apply_ops(state, [{"op": "reassign", "atom_names": ["C9", "H18", "H19"],
                                  "from_role": 5, "to_role": 4}])
    assert out.n_beads == 120
    assert sum(len(b["atom_indices"]) for b in out.beads) == 782  # partition intact
    # every monomer's ammonium bead gained a heavy atom (4 -> 5), propyl lost one
    for m in range(20):
        assert out.bead(6 * m + 4)["heavy_atom_count"] == 5
        assert out.bead(6 * m + 5)["heavy_atom_count"] == 2


def test_apply_ops_tiled_merge_shrinks_every_monomer(state):
    out = loop.apply_ops(state, [{"op": "merge", "roles": [2, 3]}])
    assert out.n_beads == 100                       # 120 - 20
    assert sum(len(b["atom_indices"]) for b in out.beads) == 782
    # merged ester bead is 2+2 = 4 heavy in every (now 5-bead) monomer
    for m in range(20):
        assert out.bead(5 * m + 2)["heavy_atom_count"] == 4


def test_apply_ops_rejects_bad_role(state):
    with pytest.raises(ValueError, match="role"):
        loop.apply_ops(state, [{"op": "reassign", "atom_names": ["C9"],
                                "from_role": 5, "to_role": 99}])


def test_apply_ops_auto_carries_hydrogens(state):
    ref = repair.load_ref_positions(REPO / "reference/PSBMA_20mer/PSBMA_20mer_no_water.gro")
    base_n3 = len(state.bead(3)["atom_indices"])          # ester bead: O2,C5,H8,H9 (4)
    # name only the heavy carbon C6 (a CH2); its 2 hydrogens should follow via ref
    out = loop.apply_ops(state, [{"op": "reassign", "atom_names": ["C6"],
                                  "from_role": 4, "to_role": 3}], ref_positions=ref)
    b3 = out.bead(3)
    assert "C6" in b3["atom_names"] and b3["heavy_atom_count"] == 3   # O2,C5,C6
    assert len(b3["atom_indices"]) == base_n3 + 3                     # C6 + its 2 H
    assert out.bead(4)["heavy_atom_count"] == 3                       # N,C7,C8 left
    assert sum(len(b["atom_indices"]) for b in out.beads) == 782      # partition intact
    # without a reference frame, only the named atom moves (no auto-carry)
    bare = loop.apply_ops(state, [{"op": "reassign", "atom_names": ["C6"],
                                   "from_role": 4, "to_role": 3}])
    assert len(bare.bead(3)["atom_indices"]) == base_n3 + 1           # just C6


# ---------- rule check ----------


def test_rule_check_flags_oversized_bead(state):
    good = loop.default_rule_check(state)
    assert good == []                                # baseline: all beads <= 4 heavy
    over = loop.apply_ops(state, [{"op": "reassign", "atom_names": ["C9", "H18", "H19"],
                                   "from_role": 5, "to_role": 4}])
    bad = loop.default_rule_check(over)
    assert len(bad) == 20 and "heavy atoms (>4)" in bad[0]


# ---------- controller ----------


def test_run_loop_accepts_improving_merge(state):
    policy = loop.ScriptedPolicy([
        loop.Action(thought="merge the two over-split ester T-beads", ops=[{"op": "merge", "roles": [2, 3]}]),
    ])
    best, result = loop.run_loop(state, policy, _fake_eval(_sq_dev), max_iters=6)
    # baseline dev^2 = 199 (monomer 20's backbone bead is 4-heavy: extra BR1 terminator)
    assert result.initial_objective == 199.0
    assert result.best_objective == 39.0            # after merging the two 2-heavy ester beads
    assert best.n_beads == 100
    assert result.stop_reason == "submit"
    step = result.trajectory[0]
    assert step.accepted and step.improved
    assert result.trajectory[-1].note == "submit"


def test_run_loop_rejects_rule_violation_keeps_best(state):
    policy = loop.ScriptedPolicy([
        loop.Action(thought="stuff C9 into the ammonium bead", ops=[{"op": "reassign",
            "atom_names": ["C9", "H18", "H19"], "from_role": 5, "to_role": 4}]),
    ])
    best, result = loop.run_loop(state, policy, _fake_eval(_sq_dev), max_iters=6)
    assert result.best_objective == 199.0           # rejected → best unchanged
    assert best.n_beads == 120                       # working never advanced
    assert result.trajectory[0].accepted is False
    assert "REJECTED (Martini rule)" in result.trajectory[0].note
    assert result.stop_reason == "submit"            # scripted policy then submits


def test_run_loop_hillclimb_discards_regression(state):
    # objective worsens as beads shrink, so a merge (120->100) is non-improving.
    # hill-climb must DISCARD it and propose the next edit from the pristine 120-bead
    # state — so step 2's merge again scores the 120->100 value (0.7), not 100->80 (0.9).
    obj = _fake_eval(lambda st: {120: 0.5, 100: 0.7, 80: 0.9}.get(st.n_beads, 1.0))
    policy = loop.ScriptedPolicy([
        loop.Action(ops=[{"op": "merge", "roles": [2, 3]}]),
        loop.Action(ops=[{"op": "merge", "roles": [2, 3]}]),
    ])
    best, result = loop.run_loop(state, policy, obj, max_iters=6, plateau_k=2)
    assert result.best_objective == 0.5 and best.n_beads == 120     # never left the best
    assert result.stop_reason == "plateau" and result.n_iterations == 2
    for step in result.trajectory:
        assert step.accepted is False and step.improved is False    # discarded, not adopted
        assert step.objective == 0.7                                # both scored 120->100 => reverted


def test_run_loop_plateau_stops(state):
    # constant objective → every valid step is non-improving → plateau after K
    op = {"op": "reassign", "atom_names": ["C2"], "from_role": 1, "to_role": 2}
    policy = loop.ScriptedPolicy([loop.Action(ops=[op]) for _ in range(6)])
    _, result = loop.run_loop(state, policy, _fake_eval(lambda st: 1.0),
                              max_iters=20, plateau_k=3)
    assert result.stop_reason == "plateau"
    assert result.n_iterations == 3


def test_run_loop_budget_stops(state):
    op = {"op": "reassign", "atom_names": ["C2"], "from_role": 1, "to_role": 2}
    policy = loop.ScriptedPolicy([loop.Action(ops=[op]) for _ in range(10)])
    _, result = loop.run_loop(state, policy, _fake_eval(lambda st: 1.0),
                              max_iters=2, plateau_k=99)
    assert result.stop_reason == "budget"
    assert result.n_iterations == 2
