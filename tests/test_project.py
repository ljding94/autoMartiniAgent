"""Smoke tests for ``agent.project``.

Strategy: project the equilibrated single-frame ``.gro`` through the
PEO-20 mapping and assert (a) the output shape is right, (b) the bead
positions match a from-scratch mass-weighted COM computed in plain
Python with no MDAnalysis. Catches index off-by-ones and mass-weighting
regressions — the two ways the projector can silently lie.
"""

from __future__ import annotations

import math
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pytest

from agent.project import (
    _build_bead_groups,
    _build_cg_universe,
    load_mapping,
    project_trajectory,
)

REPO = Path(__file__).resolve().parents[1]
AA_GRO = REPO / "reference/gromacs/equil.gro"
AA_TPR = REPO / "reference/gromacs/equil.tpr"
MAPPING = REPO / "reference/polyply_PEO20/PEO20_mapping.json"

ATOMIC_MASS = {"O": 15.9994, "C": 12.011, "H": 1.008, "N": 14.007, "S": 32.06}


@pytest.fixture(scope="module")
def mapping() -> dict:
    return load_mapping(MAPPING)


@pytest.fixture(scope="module")
def projected(tmp_path_factory) -> dict:
    out_dir = tmp_path_factory.mktemp("PEO20")
    out_traj = out_dir / "PEO20_cg.xtc"
    out_struct = out_dir / "PEO20_cg.gro"
    result = project_trajectory(
        aa_top=AA_TPR,
        aa_traj=AA_GRO,
        mapping=MAPPING,
        out_traj=out_traj,
        out_struct=out_struct,
    )
    return {"result": result, "traj": out_traj, "struct": out_struct}


def _parse_gro_peo(gro: Path) -> dict[int, dict]:
    """Pull only ETHOX atoms from a .gro file, keyed by 1-based atom index."""
    text = gro.read_text().splitlines()
    n = int(text[1].strip())
    coords = {}
    for line in text[2 : 2 + n]:
        resname = line[5:10].strip()
        if resname != "ETHOX":
            continue
        idx = int(line[15:20])
        aname = line[10:15].strip()
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        coords[idx] = {"name": aname, "xyz": (x, y, z)}
    return coords


def _python_com(mapping: dict, coords: dict[int, dict]) -> list[tuple[float, float, float]]:
    out = []
    for bead in mapping["beads"]:
        ms = mx = my = mz = 0.0
        for idx, name in zip(bead["atom_indices"], bead["atom_names"]):
            x, y, z = coords[idx]["xyz"]
            m = ATOMIC_MASS[name[0].upper()]
            ms += m
            mx += m * x
            my += m * y
            mz += m * z
        out.append((mx / ms, my / ms, mz / ms))
    return out


def test_load_mapping_has_expected_shape(mapping):
    assert mapping["molecule"] == "PEO20"
    assert len(mapping["beads"]) == 20
    assert all(b["bead_type"] == "SN3r" for b in mapping["beads"])
    assert mapping["beads"][0]["aa_residue"] == 1
    assert mapping["beads"][-1]["aa_residue"] == 20


def test_bead_groups_match_mapping_atom_counts(mapping):
    u = mda.Universe(str(AA_TPR), str(AA_GRO))
    groups = _build_bead_groups(u, mapping)
    assert len(groups) == 20
    for bg, bead in zip(groups, mapping["beads"]):
        assert len(bg) == len(bead["atom_indices"])


def test_cg_universe_metadata(mapping):
    cg_u = _build_cg_universe(mapping)
    assert len(cg_u.atoms) == 20
    assert list(cg_u.atoms.names) == [b["bead_name"] for b in mapping["beads"]]
    assert list(cg_u.atoms.types) == [b["bead_type"] for b in mapping["beads"]]
    np.testing.assert_allclose(
        cg_u.atoms.masses,
        [b["cg_mass"] for b in mapping["beads"]],
    )


def test_projected_output_shape(projected):
    cg_u = mda.Universe(str(projected["struct"]), str(projected["traj"]))
    assert len(cg_u.atoms) == 20
    assert len(cg_u.trajectory) == 1
    assert projected["result"].n_frames == 1
    assert projected["result"].n_beads == 20


def test_projected_positions_match_python_com(projected, mapping):
    cg_u = mda.Universe(str(projected["struct"]), str(projected["traj"]))
    aa_coords = _parse_gro_peo(AA_GRO)
    expected = _python_com(mapping, aa_coords)

    cg_u.trajectory[0]
    got = cg_u.atoms.positions / 10.0  # Å → nm

    for i, (gx, gy, gz) in enumerate(expected):
        dx, dy, dz = got[i] - np.array([gx, gy, gz])
        d = math.sqrt(dx * dx + dy * dy + dz * dz)
        assert d < 1e-3, f"bead {i + 1} drift {d:.4f} nm exceeds tolerance"
