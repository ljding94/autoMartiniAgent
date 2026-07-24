"""Tests for the rank validation harness (``agent.validate``).

Trajectory-free: the perturbations (named + random scramble) are pure mapping
transforms, so they are checked for validity (partition intact, deterministic)
without projecting; the rank driver is exercised with a fake evaluator.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from agent import loop, repair, validate

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


def _atoms(st):
    return sum(len(b["atom_indices"]) for b in st.beads)


def test_scramble_preserves_partition(state, ref):
    s = validate.scramble(state, 5, seed=0, ref_positions=ref)
    assert s.n_beads == 120                 # reassignments keep bead count fixed
    assert _atoms(s) == 782
    repair.validate_state(s)                # no lost atoms, no empty beads


def test_scramble_is_deterministic(state):
    a = validate.scramble(state, 4, seed=3)
    b = validate.scramble(state, 4, seed=3)
    assert [bd["atom_indices"] for bd in a.beads] == [bd["atom_indices"] for bd in b.beads]


def test_scramble_actually_moves_atoms(state, ref):
    s = validate.scramble(state, 6, seed=1, ref_positions=ref)
    changed = any(a["atom_indices"] != b["atom_indices"]
                  for a, b in zip(state.beads, s.beads))
    assert changed


def test_named_perturbations_are_valid(state, ref):
    for name, ops in validate.NAMED_PERTURBATIONS_PSBMA.items():
        s = loop.apply_ops(state, ops, ref_positions=ref)
        assert _atoms(s) == 782, name
        repair.validate_state(s)


def test_run_recover_mechanics(state, ref):
    # perturbation: split_ammonium (C7 role4->role5). A scripted policy that moves
    # C7 back (role5->role4) should fully recover. Fake objective = # atoms off their
    # good-mapping role (0 = perfect), so we can check recovery_fraction/grouping_match
    # without projecting or calling an LLM.
    good_roles = validate._atom_roles(state)

    def fake(st):
        rr = validate._atom_roles(st)
        off = sum(1 for a, r in good_roles.items() if rr.get(a) != r)
        return SimpleNamespace(objective=float(off), breakdown=SimpleNamespace(terms=[]))

    def make_policy():
        return loop.ScriptedPolicy([loop.Action(
            ops=[{"op": "reassign", "atom_names": ["C7"], "from_role": 5, "to_role": 4}])])

    rows = validate.run_recover(
        state, make_policy, fake, ref_positions=ref,
        perturbations={"split_ammonium": validate.NAMED_PERTURBATIONS_PSBMA["split_ammonium"]},
        restarts=1, max_iters=3, plateau_k=2)
    r = rows[0]
    assert r.good == 0.0                 # good mapping matches itself
    assert r.degraded > 0.0             # the split moved atoms off-role
    assert r.recovered == 0.0          # scripted inverse fully recovers
    assert r.recovery_fraction == pytest.approx(1.0)
    assert r.grouping_match == pytest.approx(1.0)


def test_run_rank_shape(state):
    # fake objective: total atoms in role-1 beads (varies as scramble moves atoms)
    def fake(st):
        obj = sum(len(st.bead(6 * m + 1)["atom_indices"]) for m in range(20))
        return SimpleNamespace(objective=float(obj), breakdown=SimpleNamespace(terms=[]))

    res = validate.run_rank(state, fake, severities=(1, 2), seeds=(0, 1))
    assert res.baseline > 0
    n_named = len(validate.NAMED_PERTURBATIONS_PSBMA)
    assert len(res.rows) == n_named + 2 * 2          # named + severities×seeds
    assert all(r.delta == pytest.approx(r.objective - res.baseline) for r in res.rows)
