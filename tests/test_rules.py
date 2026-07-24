"""Tests for the Martini rule checker (``agent.rules``). Trajectory-free."""

from __future__ import annotations

from pathlib import Path

import pytest

from agent import loop, repair, rules

REPO = Path(__file__).resolve().parents[1]
MAP = REPO / "reference/PSBMA_20mer/PSBMA20_mapping.json"
ITP = REPO / "reference/PSBMA_20mer/PSBMA20.itp"
GRO = REPO / "reference/PSBMA_20mer/PSBMA_20mer_no_water.gro"


@pytest.fixture(scope="module")
def state() -> repair.MappingState:
    return repair.load_state(MAP, ITP)


@pytest.fixture(scope="module")
def ref():
    return repair.load_ref_positions(GRO)


def test_size_class_and_expected_heavy():
    assert rules.size_class("TN5a") == "T" and rules.expected_heavy("TN5a") == 2
    assert rules.size_class("SC1") == "S" and rules.expected_heavy("SC1") == 3
    assert rules.size_class("Q1") == "R" and rules.expected_heavy("Q1") == 4


def test_baseline_has_no_hard_errors(state):
    assert rules.hard_violations(state) == []          # the good mapping is Martini-error-free


def test_baseline_flags_the_br1_end_bead(state):
    # the ATRP Br terminator makes monomer 20's backbone bead 4-heavy (SC1 = S expects 3)
    sizing = [v for v in rules.check_rules(state) if v.kind == "sizing"]
    assert sizing, "expected a size-class warning from the 4-heavy Br-terminated end bead"


def test_oversize_edit_is_a_hard_error(state, ref):
    over = loop.apply_ops(state, [{"op": "reassign", "atom_names": ["C9"],
                                   "from_role": 5, "to_role": 4}], ref_positions=ref)
    errs = rules.hard_violations(over)
    assert errs and any("> 4" in e for e in errs)      # the 5-heavy ammonium bead


def test_relabel_to_fit_clears_sizing_warnings_objective_neutrally(state):
    before = sum(1 for v in rules.check_rules(state) if v.kind == "sizing")
    fixed = rules.relabel_to_fit(state)
    after = sum(1 for v in rules.check_rules(fixed) if v.kind == "sizing")
    assert after < before                              # relabelling fixes size-class mismatches
    # objective-neutral: the atom partition is untouched (only labels change)
    assert [b["atom_indices"] for b in fixed.beads] == [b["atom_indices"] for b in state.beads]


def test_fg_integrity_flags_split_sulfonate(state, ref):
    # the good mapping keeps functional groups intact; moving a sulfonate O out of
    # the SO3 bead splits the group (this is exactly greedy's objective-gaming edit)
    assert rules.fg_violations(state, rules.FUNCTIONAL_GROUPS_PSBMA) == []
    split = loop.apply_ops(state, [{"op": "reassign", "atom_names": ["O3"],
                                    "from_role": 6, "to_role": 5}], ref_positions=ref)
    bad = rules.fg_violations(split, rules.FUNCTIONAL_GROUPS_PSBMA)
    assert bad and all(v.kind == "fg_split" and v.severity == "error" for v in bad)
    assert any("sulfonate" in v.message for v in bad)
    # and it becomes a hard-constraint error when FGs are supplied
    assert rules.hard_violations(split, rules.FUNCTIONAL_GROUPS_PSBMA)


def test_charge_on_non_q_bead_is_an_error(state):
    bad = repair._clone(state)
    bad.bead(2)["charge"] = 1                           # TN5a is not a Q-type
    errs = [v for v in rules.check_rules(bad) if v.kind == "charge_type" and v.severity == "error"]
    assert errs
