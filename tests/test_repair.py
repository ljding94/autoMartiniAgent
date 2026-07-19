"""Tests for the mapping mutation verbs (``agent.repair``) and the repair
objective (``agent.evaluate.mapping_error``).

These are trajectory-free: the verbs and the objective are pure functions of the
mapping JSON + itp, so the suite runs without the (gitignored) AA trajectory. The
project+score+objective integration is exercised interactively in the repair loop,
not here.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from agent import repair
from agent.evaluate import mapping_error
from agent.score import ScoreReport, TermStats

REPO = Path(__file__).resolve().parents[1]
MAP = REPO / "reference/PSBMA_20mer/PSBMA20_mapping.json"
ITP = REPO / "reference/PSBMA_20mer/PSBMA20.itp"
REPAIRED = REPO / "derived/PSBMA20/repair/PSBMA20_repaired_mapping.json"
REPAIRED_ITP = REPO / "derived/PSBMA20/repair/PSBMA20_repaired.itp"


@pytest.fixture(scope="module")
def state() -> repair.MappingState:
    return repair.load_state(MAP, ITP)


def _total_mass(st: repair.MappingState) -> float:
    return sum(b["aa_mass_sum"] for b in st.beads)


def _total_atoms(st: repair.MappingState) -> int:
    return sum(len(b["atom_indices"]) for b in st.beads)


def test_load_shape(state):
    assert state.n_beads == 120
    assert len(state.bonds) == 119          # 100 intra + 19 inter-monomer backbone
    assert len(state.all_atoms) == 782
    assert _total_atoms(state) == 782
    # bead_id == position invariant
    assert [b["bead_id"] for b in state.beads] == list(range(1, 121))


def test_reassign_conserves_atoms_and_mass(state):
    m0 = _total_mass(state)
    # C6(8) + its two H (27,28) : ammonium bead 4 -> ester bead 3 (bonded neighbour)
    r = repair.reassign_atoms(state, [8, 27, 28], 3)
    assert _total_atoms(r) == 782
    assert set(r.all_atoms) == set(state.all_atoms)
    assert math.isclose(_total_mass(r), m0, abs_tol=1e-6)   # atoms only moved
    assert 8 in r.bead(3)["atom_indices"] and 8 not in r.bead(4)["atom_indices"]
    assert r.n_beads == 120                                  # fixed bead count
    # input state untouched
    assert 8 in state.bead(4)["atom_indices"]


def test_reassign_requires_adjacency(state):
    # bead 4 (ammonium) is not bonded to bead 6 (sulfonate) -> must refuse
    assert 6 not in state.neighbors(4)
    with pytest.raises(ValueError, match="not bonded"):
        repair.reassign_atoms(state, [8], 6)
    # ... unless explicitly overridden
    repair.reassign_atoms(state, [8, 27, 28], 6, require_adjacent=False)


def test_merge_beads(state):
    m0 = _total_mass(state)
    m = repair.merge_beads(state, [3, 4])
    assert m.n_beads == 119
    assert _total_atoms(m) == 782
    assert math.isclose(_total_mass(m), m0, abs_tol=1e-6)
    # merged bead carries the union and a composite type
    assert set(m.bead(3)["atom_indices"]) == set(
        state.bead(3)["atom_indices"] + state.bead(4)["atom_indices"]
    )
    assert m.bead(3)["bead_type"] == "TP2a/Q1"
    # bond between the merged pair is gone; their external neighbours survive.
    # after reindex old bead5 (propyl) -> id 4, so merged bead 3 neighbours {2, 4}
    assert len(m.bonds) == 118
    assert m.neighbors(3) == {2, 4}     # old bead2 (carbonyl) and old bead5 (propyl)


def test_split_bead_partitions_and_bonds(state):
    ref = repair.load_ref_positions(REPO / "reference/PSBMA_20mer/PSBMA_20mer_no_water.gro")
    s = repair.split_bead(state, 4, [8, 27, 28], ref_positions=ref)  # C6 | N,C7,C8
    assert s.n_beads == 121
    assert _total_atoms(s) == 782
    # the two children partition the original bead's atoms
    child_a = set(s.bead(4)["atom_indices"])
    child_b = set(s.bead(5)["atom_indices"])
    assert child_a.isdisjoint(child_b)
    assert child_a | child_b == set(state.bead(4)["atom_indices"])
    # a fresh A-B bond exists and each child keeps one external neighbour
    assert 5 in s.neighbors(4)          # A-B bond
    assert 3 in s.neighbors(4)          # C6 -> ester (bead3)
    assert 6 in s.neighbors(5)          # N-side -> propyl (now bead6)


def test_validate_catches_lost_atom(state):
    bad = repair._clone(state)
    bad.bead(1)["atom_indices"] = bad.bead(1)["atom_indices"][1:]  # drop an atom
    with pytest.raises(ValueError, match="atom set changed"):
        repair.validate_state(bad)


def test_validate_catches_empty_bead(state):
    bad = repair._clone(state)
    # move all of bead 2's atoms into bead 1 by hand, leaving bead 2 empty
    bad.bead(1)["atom_indices"] += bad.bead(2)["atom_indices"]
    bad.bead(2)["atom_indices"] = []
    with pytest.raises(ValueError, match="empty"):
        repair.validate_state(bad)


def test_itp_and_mapping_roundtrip(state, tmp_path):
    m = repair.merge_beads(state, [3, 4])
    mp = repair.write_mapping(m, tmp_path / "m.json")
    ip = repair.write_itp(m, tmp_path / "m.itp")
    reloaded = repair.load_state(mp, ip)
    assert reloaded.n_beads == 119
    assert len(reloaded.bonds) == 118
    assert set(reloaded.all_atoms) == set(state.all_atoms)


def _term(kind_label, r2, n=1000):
    return TermStats(
        group_id=0, label=kind_label, bead_pattern=kind_label, n_members=1,
        n_observations=n, measured_mu=0.0, measured_sigma=0.0, fit_mu=0.0,
        fit_sigma=0.0, fit_amp=0.0, fit_r2=r2, target_mu=None, target_k=None,
        delta_vs_target=None,
    )


def test_mapping_error_is_mean_one_minus_r2():
    report = ScoreReport(
        molecule="X", n_frames=1, n_beads=2, end_exclude=0,
        bond_terms=[_term("b1", 1.0), _term("b2", 0.90)],
        angle_terms=[_term("a1", 0.80)],
    )
    b = mapping_error(report)
    assert b.n_terms == 3
    # mean of (0.0, 0.10, 0.20)
    assert math.isclose(b.objective, 0.10, abs_tol=1e-9)
    assert math.isclose(b.bond_mean, 0.05, abs_tol=1e-9)
    assert math.isclose(b.angle_mean, 0.20, abs_tol=1e-9)
    assert b.worst[0].label == "a1"          # largest error first


def test_mapping_error_skips_nan():
    report = ScoreReport(
        molecule="X", n_frames=1, n_beads=2, end_exclude=0,
        bond_terms=[_term("b1", float("nan")), _term("b2", 0.90)],
    )
    b = mapping_error(report)
    assert b.n_terms == 1
    assert math.isclose(b.objective, 0.10, abs_tol=1e-9)


@pytest.mark.skipif(not REPAIRED.exists(), reason="repaired PSBMA deliverable not generated")
def test_repaired_deliverable_is_valid():
    st = repair.load_state(REPAIRED, REPAIRED_ITP)
    assert st.n_beads == 120                 # W2 keeps the bead count
    assert len(st.all_atoms) == 782
    # every bead stays within Martini regular-bead sizing (<= 4 heavy atoms)
    assert all(b["heavy_atom_count"] <= 4 for b in st.beads)
