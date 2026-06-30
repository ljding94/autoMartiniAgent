"""Smoke tests for ``agent.score``.

Strategy: read the canonical PEO-20 ``.itp``, score the committed CG
trajectory snapshot in ``derived/PEO20_solu``, and assert that
  - the bonded skeleton parsed out of the ``.itp`` has the expected
    cardinality (19 bonds, 18 angles for a 20-bead chain),
  - the fitted means agree with the values reported in the README
    table (sanity-check against drift in the Gaussian fitter),
  - target μ and Δ-vs-target are wired correctly.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from agent.score import (
    load_topology,
    score_mapping,
    write_json_report,
)

REPO = Path(__file__).resolve().parents[1]
ITP = REPO / "reference/polyply_PEO20/PEO20.itp"
CG_GRO = REPO / "derived/PEO20_solu/PEO20_cg.gro"
CG_XTC = REPO / "derived/PEO20_solu/PEO20_cg.xtc"


pytestmark = pytest.mark.skipif(
    not (CG_GRO.exists() and CG_XTC.exists()),
    reason="derived/PEO20_solu CG trajectory not present (run agent.project first)",
)


@pytest.fixture(scope="module")
def topology():
    return load_topology(ITP)


@pytest.fixture(scope="module")
def report():
    return score_mapping(
        itp=ITP,
        cg_struct=CG_GRO,
        cg_traj=CG_XTC,
        molecule="PEO20",
        end_exclude=2,
    )


def test_topology_shape(topology):
    assert len(topology.bonds) == 19
    assert len(topology.angles) == 18
    assert topology.bonds[0].b0 == pytest.approx(0.36)
    assert topology.bonds[0].kb == pytest.approx(7000.0)
    assert topology.angles[0].theta0 == pytest.approx(123.0)
    assert topology.angles[0].ka == pytest.approx(80.0)


def test_report_dimensions(report):
    assert report.molecule == "PEO20"
    assert report.n_beads == 20
    assert report.n_frames == 2000
    assert report.end_exclude == 2
    assert report.bond is not None
    assert report.angle is not None


def test_end_exclude_applied(report):
    # 19 raw bonds, drop those touching beads 0,1,18,19 → keep 15
    assert report.bond.n_members == 15
    # 18 raw angles, drop those touching beads 0,1,18,19 → keep 14
    assert report.angle.n_members == 14
    assert report.bond.n_observations == report.n_frames * report.bond.n_members
    assert report.angle.n_observations == report.n_frames * report.angle.n_members


def test_targets_pulled_from_itp(report):
    assert report.bond.target_mu == pytest.approx(0.36)
    assert report.angle.target_mu == pytest.approx(123.0)


def test_fit_means_in_expected_band(report):
    # Locked from the README PEO-20 first-numbers table; tolerances are
    # generous so a Gaussian-fitter tweak doesn't fail this trivially.
    assert 0.32 < report.bond.fit_mu < 0.34
    assert 125 < report.angle.fit_mu < 135
    # Δ vs target should have the right sign + ballpark
    assert -0.05 < report.bond.delta_vs_target < -0.01
    assert 5 < report.angle.delta_vs_target < 11


def test_json_round_trip(report, tmp_path):
    import json
    out = tmp_path / "score_report.json"
    write_json_report(report, out)
    obj = json.loads(out.read_text())
    assert obj["molecule"] == "PEO20"
    assert obj["bond"]["target_mu"] == pytest.approx(0.36)
    assert obj["angle"]["target_mu"] == pytest.approx(123.0)
