"""Repair-loop evaluator: mapping edit -> Gaussianity error.

Glue that closes ingredient #4's inner loop. Given a ``MappingState`` (or a
mapping+itp on disk), it:

  1. writes the mapping JSON + a parameter-free ``.itp`` (``agent.repair``),
  2. projects the AA trajectory through it (``agent.project``, optional stride),
  3. scores every bond + every bonded angle (``agent.score`` all-bonded mode),
  4. reduces the per-term R^2 to a single scalar **error** the agent minimizes.

The objective is deliberately *target-free*: it measures how Gaussian (harmonic-
consistent) the AA-projected distributions are, not how close they sit to any
force-field value. That is exactly "make the measured CG distributions as Gaussian
as possible".

    error = mean over (all bonds + all bonded angles) of (1 - R^2_fit)

R^2 = 1 is a perfect single-Gaussian trace, so error = 0 is the ideal and larger
is worse (a bimodal term contributes > 1). Results are cached by a content hash of
(mapping, itp, stride, trajectory identity) so re-evaluating an identical state is
free.

CLI::

  python -m agent.evaluate \\
    --mapping reference/PSBMA_20mer/PSBMA20_mapping.json \\
    --itp     reference/PSBMA_20mer/PSBMA20.itp \\
    --aa-top  reference/PSBMA_20mer/PSBMA_20mer_no_water.gro \\
    --aa-traj reference/PSBMA_20mer/PSBMA_20mer_no_water_skip10.xtc \\
    --frame-stride 25
"""

from __future__ import annotations

import argparse
import glob as _glob
import hashlib
import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Sequence

from agent import repair as _repair
from agent.project import project_trajectory
from agent.score import ScoreReport, score_mapping, write_plots


@dataclass(frozen=True)
class TermError:
    kind: str          # "bond" | "angle"
    label: str
    r2: float
    error: float       # 1 - r2
    n_obs: int


@dataclass(frozen=True)
class ObjectiveBreakdown:
    """The scalar the agent minimizes, plus per-category means and worst terms."""
    objective: float               # mean (1-R^2) over all scored terms
    n_terms: int
    bond_mean: float
    angle_mean: float
    terms: list[TermError] = field(default_factory=list)

    @property
    def worst(self) -> list[TermError]:
        return sorted(self.terms, key=lambda t: t.error, reverse=True)


def mapping_error(report: ScoreReport) -> ObjectiveBreakdown:
    """Reduce a ScoreReport to the scalar Gaussianity error (lower = better)."""
    terms: list[TermError] = []
    for t in report.bond_terms:
        if not math.isnan(t.fit_r2):
            terms.append(TermError("bond", t.label, t.fit_r2, 1.0 - t.fit_r2, t.n_observations))
    for t in report.angle_terms:
        if not math.isnan(t.fit_r2):
            terms.append(TermError("angle", t.label, t.fit_r2, 1.0 - t.fit_r2, t.n_observations))

    bond_errs = [t.error for t in terms if t.kind == "bond"]
    angle_errs = [t.error for t in terms if t.kind == "angle"]
    all_errs = [t.error for t in terms]
    return ObjectiveBreakdown(
        objective=(sum(all_errs) / len(all_errs)) if all_errs else float("nan"),
        n_terms=len(all_errs),
        bond_mean=(sum(bond_errs) / len(bond_errs)) if bond_errs else float("nan"),
        angle_mean=(sum(angle_errs) / len(angle_errs)) if angle_errs else float("nan"),
        terms=terms,
    )


@dataclass
class EvalResult:
    objective: float
    breakdown: ObjectiveBreakdown
    report: ScoreReport
    n_frames: int
    n_beads: int
    work_dir: Path
    mapping_path: Path
    itp_path: Path
    plots: dict[str, Path] = field(default_factory=dict)


def _resolve_traj(aa_traj: str | Path | Sequence[str | Path]) -> list[str]:
    if isinstance(aa_traj, (str, Path)):
        s = str(aa_traj)
        if any(c in s for c in "*?["):
            return sorted(_glob.glob(s))
        return [s]
    return [str(p) for p in aa_traj]


def _state_hash(state: _repair.MappingState, frame_stride: int, traj: list[str]) -> str:
    h = hashlib.sha1()
    h.update(json.dumps(state.mapping, sort_keys=True).encode())
    bonds = sorted((b.key()[0], b.key()[1], b.b0, b.kb) for b in state.bonds)
    h.update(json.dumps(bonds).encode())
    h.update(str(frame_stride).encode())
    h.update("|".join(traj).encode())
    return h.hexdigest()[:12]


def evaluate_state(
    state: _repair.MappingState,
    *,
    aa_top: str | Path,
    aa_traj: str | Path | Sequence[str | Path],
    work_root: str | Path,
    molecule: str | None = None,
    frame_stride: int = 1,
    make_plots: bool = False,
    use_cache: bool = True,
) -> EvalResult:
    """Project + score a MappingState and return its Gaussianity error.

    Artifacts (mapping.json, topo.itp, cg.xtc/gro, score_report.json, optional
    PDFs) land in ``work_root/<hash>/`` so distinct edits never clobber each other
    and identical states are served from cache.
    """
    traj = _resolve_traj(aa_traj)
    mol = molecule or state.mapping.get("molecule", "MOL")
    tag = _state_hash(state, frame_stride, traj)
    work_dir = Path(work_root) / tag
    work_dir.mkdir(parents=True, exist_ok=True)

    mapping_path = _repair.write_mapping(state, work_dir / "mapping.json")
    itp_path = _repair.write_itp(state, work_dir / "topo.itp", mol_name=mol)
    report_path = work_dir / "score_report.json"
    cg_xtc = work_dir / "cg.xtc"
    cg_gro = work_dir / "cg.gro"

    if use_cache and report_path.exists() and cg_gro.exists():
        report = _load_report(report_path)
    else:
        project_trajectory(
            aa_top=aa_top,
            aa_traj=traj if len(traj) > 1 else traj[0],
            mapping=json.loads(mapping_path.read_text()),
            out_traj=cg_xtc,
            out_struct=cg_gro,
            frame_stride=frame_stride,
        )
        report = score_mapping(
            itp=itp_path,
            cg_struct=cg_gro,
            cg_traj=cg_xtc,
            molecule=mol,
            all_bonded_angles=True,
        )
        report_path.write_text(json.dumps(asdict(report), indent=2))

    breakdown = mapping_error(report)
    plots: dict[str, Path] = {}
    if make_plots:
        plots = write_plots(
            report, itp=itp_path, cg_struct=cg_gro, cg_traj=cg_xtc,
            out_dir=work_dir, all_bonded_angles=True,
        )

    return EvalResult(
        objective=breakdown.objective,
        breakdown=breakdown,
        report=report,
        n_frames=report.n_frames,
        n_beads=report.n_beads,
        work_dir=work_dir,
        mapping_path=mapping_path,
        itp_path=itp_path,
        plots=plots,
    )


def evaluate_paths(
    mapping_path: str | Path,
    itp_path: str | Path,
    *,
    aa_top: str | Path,
    aa_traj: str | Path | Sequence[str | Path],
    work_root: str | Path,
    molecule: str | None = None,
    frame_stride: int = 1,
    make_plots: bool = False,
    atom_masses: dict[int, float] | None = None,
) -> EvalResult:
    """Convenience: load a mapping+itp from disk and evaluate it."""
    state = _repair.load_state(mapping_path, itp_path, atom_masses=atom_masses)
    return evaluate_state(
        state, aa_top=aa_top, aa_traj=aa_traj, work_root=work_root,
        molecule=molecule, frame_stride=frame_stride, make_plots=make_plots,
    )


def _load_report(path: str | Path) -> ScoreReport:
    from agent.score import AngleCoverage, TermStats
    d = json.loads(Path(path).read_text())
    return ScoreReport(
        molecule=d["molecule"], n_frames=d["n_frames"], n_beads=d["n_beads"],
        end_exclude=d["end_exclude"],
        bond_terms=[TermStats(**t) for t in d["bond_terms"]],
        angle_terms=[TermStats(**t) for t in d["angle_terms"]],
        angle_coverage=[AngleCoverage(**c) for c in d.get("angle_coverage", [])],
    )


def write_comparison_plots(
    *,
    itp: str | Path,
    struct_a: str | Path,
    traj_a: str | Path,
    struct_b: str | Path,
    traj_b: str | Path,
    out_dir: str | Path,
    label_a: str = "baseline",
    label_b: str = "repaired",
    bond_binwidth_nm: float = 0.001,
    angle_binwidth_deg: float = 1.0,
    name_suffix: str = "_compare",
) -> dict[str, Path]:
    """Overlay two mappings' measured distributions per bond / angle group.

    Both mappings must share the same bead-type / bond topology (true for any
    ``reassign_atoms`` edit — bead count, types and bonds are unchanged, only atom
    membership differs), so groups line up one-to-one and the same CG bead index is
    the same conceptual bead in both trajectories. Each panel draws mapping A vs
    mapping B (measured histogram + its Gaussian fit) with both R² values, so the
    before/after change in Gaussianity is visible at a glance.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mdtraj as md
    import numpy as np
    from agent.score import (
        _all_bonded_angle_groups, _drop_end_excluded, _fit_gaussian, _gaussian,
        _group_bonds, _hist, load_topology,
    )

    topology = load_topology(itp)
    ta = md.load(str(traj_a), top=str(struct_a))
    tb = md.load(str(traj_b), top=str(struct_b))
    n_beads = ta.n_atoms
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}

    specs = {
        "bond": dict(
            groups=_group_bonds(topology),
            idx=lambda ijk: [(i - 1, j - 1) for (i, j) in ijk],
            measure=lambda tr, kept: md.compute_distances(tr, kept, opt=True).ravel(),
            binwidth=bond_binwidth_nm, lo=0.0, hi=0.8, xlim=(0.15, 0.42),
            xlabel="bond distance (nm)", unit="nm", fmt=".4f",
        ),
        "angle": dict(
            groups=_all_bonded_angle_groups(topology),
            idx=lambda ijk: [(i - 1, j - 1, k - 1) for (i, j, k) in ijk],
            measure=lambda tr, kept: (180.0 / np.pi) * md.compute_angles(tr, kept, opt=True).ravel(),
            binwidth=angle_binwidth_deg, lo=0.0, hi=180.0, xlim=(40.0, 180.0),
            xlabel="bond angle (deg)", unit="°", fmt=".1f",
        ),
    }

    for kind, sp in specs.items():
        groups = sp["groups"]
        if not groups:
            continue
        n = len(groups)
        ncol = 2 if n > 1 else 1
        nrow = (n + ncol - 1) // ncol
        fig, axes = plt.subplots(nrow, ncol, figsize=(6.5 * ncol, 3.8 * nrow), squeeze=False)
        for gi, (_key, ijk_1, label, _pat, target_mu, _tk) in enumerate(groups):
            ax = axes[gi // ncol][gi % ncol]
            kept = _drop_end_excluded(sp["idx"](ijk_1), n_beads, 0)
            for tr, lab, color in ((ta, label_a, "0.45"), (tb, label_b, "tab:red")):
                vals = sp["measure"](tr, kept)
                centers, density = _hist(vals, sp["binwidth"], sp["lo"], sp["hi"])
                popt, _fit, r2 = _fit_gaussian(centers, density)
                fitted = _gaussian(centers, *popt)
                ax.plot(centers, density, lw=1.6, color=color,
                        label=f"{lab}: μ={format(popt[0], sp['fmt'])}{sp['unit']}  R²={r2:.3f}")
                ax.plot(centers, fitted, ls="--", lw=1.0, color=color, alpha=0.7)
            if target_mu is not None:
                ax.axvline(target_mu, color="green", ls=":", lw=1.1, label=f"target={target_mu:g}")
            ax.set_title(f"{label}  (n_obs={len(kept) * ta.n_frames})", fontsize=10)
            ax.set_xlabel(sp["xlabel"], fontsize=9)
            ax.set_ylabel("probability density", fontsize=9)
            ax.set_xlim(sp["xlim"])
            ax.legend(loc="upper right", fontsize=7, frameon=False)
        for j in range(n, nrow * ncol):
            axes[j // ncol][j % ncol].axis("off")
        fig.suptitle(f"{label_a} vs {label_b} — {kind} distributions", fontsize=12)
        plt.tight_layout(rect=(0, 0, 1, 0.98))
        path = out_dir / f"{kind}_hists{name_suffix}.pdf"
        fig.savefig(path)
        plt.close(fig)
        written[kind] = path
    return written


def format_breakdown(b: ObjectiveBreakdown) -> str:
    lines = [
        f"objective (mean 1-R^2) = {b.objective:.4f}  over {b.n_terms} terms "
        f"(bonds {b.bond_mean:.4f}, angles {b.angle_mean:.4f})",
        "  worst terms:",
    ]
    for t in b.worst[:6]:
        lines.append(f"    {t.kind:<5} {t.label:<20} R²={t.r2:.3f}  err={t.error:.3f}")
    return "\n".join(lines)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Evaluate a CG mapping's Gaussianity error")
    p.add_argument("--mapping", required=True)
    p.add_argument("--itp", required=True)
    p.add_argument("--aa-top", required=True)
    p.add_argument("--aa-traj", required=True, nargs="+")
    p.add_argument("--work-root", default=None, help="artifact dir (default: derived/<mol>/repair)")
    p.add_argument("--molecule", default=None)
    p.add_argument("--frame-stride", type=int, default=1)
    p.add_argument("--plots", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    state = _repair.load_state(args.mapping, args.itp)
    mol = args.molecule or state.mapping.get("molecule", "MOL")
    work_root = Path(args.work_root or f"derived/{mol}/repair")
    result = evaluate_state(
        state, aa_top=args.aa_top, aa_traj=args.aa_traj, work_root=work_root,
        molecule=mol, frame_stride=args.frame_stride, make_plots=args.plots,
    )
    print(f"evaluated {mol}: {result.n_frames} frames, {result.n_beads} beads "
          f"(stride={args.frame_stride})")
    print(format_breakdown(result.breakdown))
    print(f"  artifacts: {result.work_dir}")


if __name__ == "__main__":
    main()
