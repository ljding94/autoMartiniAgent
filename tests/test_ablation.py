"""Tests for the ablation baselines (``agent.ablation``).

Trajectory-free / LLM-free: the greedy and random searches are exercised with a
fake evaluator; the LLM-agent leg is not unit-tested (it needs the API).
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from agent import ablation, loop, repair, validate

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


def test_candidate_moves_are_bonded_single_reassigns(state):
    moves = ablation.candidate_moves(state)
    assert moves and all(m["op"] == "reassign" and len(m["atom_names"]) == 1 for m in moves)
    _n, bpm = loop._beads_per_monomer(state)
    adj = validate._bonded_role_adjacency(state, bpm)
    assert all(m["to_role"] in adj[m["from_role"]] for m in moves)


def _fake_role1_heavy_target(target_heavy=2):
    # objective = sum over monomers of |role-1 heavy count - target|
    def fake(st):
        obj = sum(abs(st.bead(6 * m + 1)["heavy_atom_count"] - target_heavy) for m in range(20))
        return SimpleNamespace(objective=float(obj), breakdown=SimpleNamespace(terms=[]))
    return fake


def test_greedy_takes_an_improving_move(state, ref):
    fake = _fake_role1_heavy_target(2)          # backbone role-1 is 3-4 heavy → shrinking it improves
    start = fake(state).objective
    r, best = ablation.greedy_search(state, fake, ref_positions=ref, max_steps=2)
    assert r.best_objective < start             # found a strictly-improving single move
    assert r.n_evaluations > 1
    assert best.n_beads == state.n_beads


def test_greedy_stops_when_stuck(state):
    fake = lambda st: SimpleNamespace(objective=1.0, breakdown=SimpleNamespace(terms=[]))
    r, _ = ablation.greedy_search(state, fake, max_steps=5)
    assert r.best_objective == 1.0              # nothing improves a constant objective
    assert r.n_evaluations >= 1


def test_random_search_respects_budget(state, ref):
    fake = lambda st: SimpleNamespace(objective=1.0, breakdown=SimpleNamespace(terms=[]))
    r, _ = ablation.random_search(state, fake, ref_positions=ref, budget=12, seed=0)
    assert r.n_evaluations <= 12
