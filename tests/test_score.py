"""Smoke tests for ``agent.score``.

Strategy: read the canonical PEO-20 ``.itp``, score the committed CG
trajectory snapshot in ``derived/PEO20_solu``, and assert that
  - the topology parses ([atoms] + [bonds] + [angles]) with the right shapes,
  - term groups (by target μ/k) come out correctly (PEO has 1 group each),
  - the fitted means agree with the values reported in the README table,
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
    assert len(topology.atoms) == 20
    assert all(a.type == "SN3r" for a in topology.atoms)
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
    assert len(report.bond_terms) == 1
    assert len(report.angle_terms) == 1


def test_end_exclude_applied(report):
    bond = report.bond_terms[0]
    angle = report.angle_terms[0]
    # 19 raw bonds, drop those touching beads 0,1,18,19 → keep 15
    assert bond.n_members == 15
    # 18 raw angles, drop those touching beads 0,1,18,19 → keep 14
    assert angle.n_members == 14
    assert bond.n_observations == report.n_frames * bond.n_members
    assert angle.n_observations == report.n_frames * angle.n_members


def test_group_labels(report):
    assert report.bond_terms[0].label == "SN3r-SN3r"
    assert report.bond_terms[0].bead_pattern == "SN3r-SN3r"
    assert report.angle_terms[0].label == "SN3r-SN3r-SN3r"


def test_targets_pulled_from_itp(report):
    assert report.bond_terms[0].target_mu == pytest.approx(0.36)
    assert report.bond_terms[0].target_k == pytest.approx(7000.0)
    assert report.angle_terms[0].target_mu == pytest.approx(123.0)
    assert report.angle_terms[0].target_k == pytest.approx(80.0)


def test_fit_means_in_expected_band(report):
    # Locked from the README PEO-20 first-numbers table; tolerances are
    # generous so a Gaussian-fitter tweak doesn't fail this trivially.
    bond = report.bond_terms[0]
    angle = report.angle_terms[0]
    assert 0.32 < bond.fit_mu < 0.34
    assert 125 < angle.fit_mu < 135
    assert -0.05 < bond.delta_vs_target < -0.01
    assert 5 < angle.delta_vs_target < 11


def test_fit_r2_is_reasonable(report):
    """PEO bond and angle distributions are near-Gaussian, so the fitted-Gaussian
    R² (in density space) should be high — 1.0 is a perfect trace and clean
    single-Gaussian fits sit well above 0.9, versus the much lower R² of the
    badly non-Gaussian (bimodal / skewed) PSBMA terms. If a scoring change
    degrades this we want to notice."""
    import math
    for t in [report.bond_terms[0], report.angle_terms[0]]:
        assert math.isfinite(t.fit_r2)
        assert t.fit_r2 <= 1.0  # R² is bounded above by 1 (perfect fit)
        # PEO fits are clean, so R² should be close to 1.
        assert t.fit_r2 > 0.9


def test_angle_coverage_always_populated(report):
    """Every score run audits which consecutive-bond angles the itp defines. PEO
    is a linear chain whose one angle type (SN3r-SN3r-SN3r) IS defined, so its
    coverage is a single in-itp entry with nothing missing. The always-on audit
    is how the workflow flags bonded angles a force field omits."""
    cov = report.angle_coverage
    assert len(cov) == 1
    assert cov[0].bead_pattern == "SN3r-SN3r-SN3r"
    assert cov[0].in_itp is True
    assert [c for c in cov if not c.in_itp] == []  # nothing missing for PEO


def test_json_round_trip(report, tmp_path):
    import json
    out = tmp_path / "score_report.json"
    write_json_report(report, out)
    obj = json.loads(out.read_text())
    assert obj["molecule"] == "PEO20"
    assert obj["bond_terms"][0]["target_mu"] == pytest.approx(0.36)
    assert obj["angle_terms"][0]["target_mu"] == pytest.approx(123.0)
    assert obj["bond_terms"][0]["label"] == "SN3r-SN3r"
