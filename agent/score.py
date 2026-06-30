"""CG mapping quality scorer.

Reads a CG trajectory + the canonical Polyply ``.itp`` for the molecule,
measures bonded distance and bond-angle distributions over all frames,
fits a Gaussian per term-type, and emits a score report comparing the
measured mean to the ``.itp`` target.

Torsions are intentionally not scored here — Martini polymer torsions are
typically weak or non-Gaussian and don't give a clean single-mode signal.

CLI::

  python -m agent.score \\
    --itp        reference/polyply_PEO20/PEO20.itp \\
    --cg-struct  derived/PEO20_solu/PEO20_cg.gro \\
    --cg-traj    derived/PEO20_solu/PEO20_cg.xtc \\
    --out-dir    derived/PEO20_solu

Outputs (in ``--out-dir``):
  - ``score_report.json``  : scalar stats per term (measured μ/σ, fit μ/σ,
                             fit RMSE, target μ, Δ vs target)
  - ``bond_hists.pdf``     : measured distribution + Gaussian fit + target line
  - ``angle_hists.pdf``    : same for the angle term
"""

from __future__ import annotations

import argparse
import glob as _glob
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import mdtraj as md
import numpy as np
from scipy.optimize import curve_fit


# ---------- helpers ----------


def _natural_key(s: str) -> list:
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", str(s))]


def _gaussian(x, mu, sigma, amp):
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def _fit_gaussian(centers: np.ndarray, density: np.ndarray):
    """Fit a single Gaussian to (centers, density). Returns (popt, fitted_y, rmse)."""
    x = np.asarray(centers, dtype=float)
    y = np.asarray(density, dtype=float)
    s = float(y.sum())
    if s > 0:
        mu0 = float(np.average(x, weights=y))
        var0 = float(np.average((x - mu0) ** 2, weights=y))
        sigma0 = max(np.sqrt(max(var0, 0.0)), float(x[1] - x[0]))
    else:
        mu0 = float(x.mean())
        sigma0 = float(x.std()) or 1.0
    amp0 = float(max(y.max(), 1e-12))
    try:
        popt, _ = curve_fit(_gaussian, x, y, p0=[mu0, sigma0, amp0], maxfev=10000)
    except Exception:
        popt = np.array([mu0, sigma0, amp0])
    fitted = _gaussian(x, *popt)
    rmse = float(np.sqrt(np.mean((y - fitted) ** 2)))
    return popt, fitted, rmse


def _moments(centers: np.ndarray, density: np.ndarray) -> tuple[float, float]:
    x = np.asarray(centers, dtype=float)
    y = np.asarray(density, dtype=float)
    if float(y.sum()) == 0:
        return float(x.mean()), float(x.std())
    mu = float(np.average(x, weights=y))
    var = float(np.average((x - mu) ** 2, weights=y))
    return mu, float(np.sqrt(max(var, 0.0)))


# ---------- .itp parsing ----------


@dataclass(frozen=True)
class BondTerm:
    i: int
    j: int
    b0: float | None
    kb: float | None
    func: int | None


@dataclass(frozen=True)
class AngleTerm:
    i: int
    j: int
    k: int
    theta0: float | None
    ka: float | None
    func: int | None


@dataclass(frozen=True)
class Topology:
    bonds: list[BondTerm]
    angles: list[AngleTerm]


def load_topology(itp_path: str | Path) -> Topology:
    """Parse [ bonds ] and [ angles ] sections from a GROMACS-style .itp file."""
    bonds: list[BondTerm] = []
    angles: list[AngleTerm] = []
    section: str | None = None
    for raw in Path(itp_path).read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        m = re.match(r"\[\s*(\w+)\s*\]", line)
        if m:
            section = m.group(1).lower()
            continue
        toks = line.split()
        try:
            if section == "bonds" and len(toks) >= 2:
                i, j = int(toks[0]), int(toks[1])
                func = int(toks[2]) if len(toks) > 2 else None
                b0 = float(toks[3]) if len(toks) > 3 else None
                kb = float(toks[4]) if len(toks) > 4 else None
                bonds.append(BondTerm(i=i, j=j, b0=b0, kb=kb, func=func))
            elif section == "angles" and len(toks) >= 3:
                i, j, k = int(toks[0]), int(toks[1]), int(toks[2])
                func = int(toks[3]) if len(toks) > 3 else None
                t0 = float(toks[4]) if len(toks) > 4 else None
                ka = float(toks[5]) if len(toks) > 5 else None
                angles.append(AngleTerm(i=i, j=j, k=k, theta0=t0, ka=ka, func=func))
        except (ValueError, IndexError):
            continue
    return Topology(bonds=bonds, angles=angles)


# ---------- core scoring ----------


@dataclass(frozen=True)
class TermStats:
    """Aggregated stats for one bonded-term type."""
    n_members: int
    n_observations: int
    measured_mu: float
    measured_sigma: float
    fit_mu: float
    fit_sigma: float
    fit_amp: float
    fit_rmse: float
    target_mu: float | None
    delta_vs_target: float | None


@dataclass(frozen=True)
class ScoreReport:
    molecule: str
    n_frames: int
    n_beads: int
    end_exclude: int
    bond: TermStats | None
    angle: TermStats | None


def _drop_end_excluded(tuples, n_beads: int, end_exclude: int):
    """Drop any term whose 0-based bead indices touch the end exclusion zone."""
    if end_exclude <= 0:
        return list(tuples)
    lo, hi = end_exclude, n_beads - end_exclude
    return [t for t in tuples if all(lo <= i < hi for i in t)]


def _hist(values: np.ndarray, binwidth: float, lo: float, hi: float):
    edges = np.arange(lo, hi + binwidth, binwidth)
    centers = 0.5 * (edges[:-1] + edges[1:])
    counts, _ = np.histogram(values, bins=edges, density=True)
    return centers, counts


def _stats_from_values(
    values: np.ndarray,
    binwidth: float,
    lo: float,
    hi: float,
    target_mu: float | None,
    n_members: int,
) -> TermStats:
    centers, density = _hist(values, binwidth, lo, hi)
    measured_mu, measured_sigma = _moments(centers, density)
    popt, _fitted, rmse = _fit_gaussian(centers, density)
    delta = float(measured_mu - target_mu) if target_mu is not None else None
    return TermStats(
        n_members=n_members,
        n_observations=int(values.size),
        measured_mu=float(measured_mu),
        measured_sigma=float(measured_sigma),
        fit_mu=float(popt[0]),
        fit_sigma=float(abs(popt[1])),
        fit_amp=float(popt[2]),
        fit_rmse=float(rmse),
        target_mu=float(target_mu) if target_mu is not None else None,
        delta_vs_target=delta,
    )


def score_mapping(
    *,
    itp: str | Path,
    cg_struct: str | Path,
    cg_traj: str | Path | Sequence[str | Path],
    molecule: str | None = None,
    end_exclude: int = 2,
    bond_binwidth_nm: float = 0.001,
    angle_binwidth_deg: float = 1.0,
) -> ScoreReport:
    """Score a CG mapping by comparing bonded distributions to ``.itp`` targets.

    Parameters
    ----------
    itp
        Canonical Polyply (or hand-written) ``.itp`` file. Provides both the
        bonded skeleton (which beads form bonds/angles) and the target values.
    cg_struct
        CG structure (``.gro``/``.pdb``) — used as topology for mdtraj.
    cg_traj
        CG trajectory; a list/glob is loaded as a single chained trajectory.
    end_exclude
        Drop any term whose bead indices touch the first/last ``end_exclude``
        beads on either end. Matches the convention used in the reference
        ``validation/PEO20/analyze_peo_model.py``.
    """
    topology = load_topology(itp)

    if isinstance(cg_traj, (str, Path)):
        traj_paths: list[str] | str = str(cg_traj)
    else:
        traj_paths = [str(p) for p in cg_traj]
    traj = md.load(traj_paths, top=str(cg_struct))
    n_beads = traj.n_atoms
    n_frames = traj.n_frames

    bond_pairs_0 = [(b.i - 1, b.j - 1) for b in topology.bonds]
    kept_bonds = _drop_end_excluded(bond_pairs_0, n_beads, end_exclude)
    bond_stats: TermStats | None = None
    if kept_bonds:
        target_b0 = topology.bonds[0].b0 if topology.bonds else None
        distances_nm = md.compute_distances(traj, kept_bonds, opt=True).ravel()
        bond_stats = _stats_from_values(
            distances_nm, bond_binwidth_nm, 0.0, 0.8, target_b0, len(kept_bonds)
        )

    angle_triplets_0 = [(a.i - 1, a.j - 1, a.k - 1) for a in topology.angles]
    kept_angles = _drop_end_excluded(angle_triplets_0, n_beads, end_exclude)
    angle_stats: TermStats | None = None
    if kept_angles:
        target_theta0 = topology.angles[0].theta0 if topology.angles else None
        angles_deg = (
            (180.0 / np.pi)
            * md.compute_angles(traj, kept_angles, opt=True).ravel()
        )
        angle_stats = _stats_from_values(
            angles_deg,
            angle_binwidth_deg,
            0.0,
            180.0,
            target_theta0,
            len(kept_angles),
        )

    return ScoreReport(
        molecule=molecule or "MOL",
        n_frames=n_frames,
        n_beads=n_beads,
        end_exclude=end_exclude,
        bond=bond_stats,
        angle=angle_stats,
    )


# ---------- output writers ----------


def write_json_report(report: ScoreReport, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(asdict(report), indent=2))


def _plot_term(
    *,
    itp: str | Path,
    cg_struct: str | Path,
    cg_traj: str | Path | Sequence[str | Path],
    stats: TermStats,
    kind: str,
    out_path: Path,
    end_exclude: int,
    binwidth: float,
    lo: float,
    hi: float,
) -> None:
    """Re-histogram for plotting (avoids carrying arrays through the JSON report)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    topology = load_topology(itp)
    if isinstance(cg_traj, (str, Path)):
        traj_paths: list[str] | str = str(cg_traj)
    else:
        traj_paths = [str(p) for p in cg_traj]
    traj = md.load(traj_paths, top=str(cg_struct))
    n_beads = traj.n_atoms

    if kind == "bond":
        terms = _drop_end_excluded(
            [(b.i - 1, b.j - 1) for b in topology.bonds], n_beads, end_exclude
        )
        values = md.compute_distances(traj, terms, opt=True).ravel()
        xlabel = "bond distance (nm)"
    elif kind == "angle":
        terms = _drop_end_excluded(
            [(a.i - 1, a.j - 1, a.k - 1) for a in topology.angles],
            n_beads,
            end_exclude,
        )
        values = (180.0 / np.pi) * md.compute_angles(traj, terms, opt=True).ravel()
        xlabel = "bond angle (deg)"
    else:
        raise ValueError(f"unknown kind {kind!r}")

    centers, density = _hist(values, binwidth, lo, hi)
    fitted = _gaussian(centers, stats.fit_mu, stats.fit_sigma, stats.fit_amp)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(centers, density, lw=1.5, label="measured")
    ax.plot(centers, fitted, "r--", lw=1.5, label="Gaussian fit")
    if stats.target_mu is not None:
        ax.axvline(
            stats.target_mu,
            color="green",
            ls=":",
            lw=1.5,
            label=f"target μ = {stats.target_mu:g}",
        )

    is_bond = kind == "bond"
    unit = "nm" if is_bond else "°"
    mu_fmt = ".4f" if is_bond else ".2f"
    annot = [
        f"μ = {stats.fit_mu:{mu_fmt}} {unit}",
        f"σ = {stats.fit_sigma:{mu_fmt}} {unit}",
        f"RMSE = {stats.fit_rmse:.4f}",
    ]
    if stats.delta_vs_target is not None:
        delta_fmt = "+" + mu_fmt
        annot.append(f"Δ vs target = {format(stats.delta_vs_target, delta_fmt)} {unit}")

    ax.text(
        0.97, 0.95,
        "\n".join(annot),
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=10,
        bbox=dict(facecolor="white", alpha=0.85, edgecolor="0.7"),
    )
    ax.set_title(f"{kind} distribution  (n_obs = {stats.n_observations})")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("probability density")
    if is_bond:
        ax.set_xlim(0.25, 0.45)
    else:
        ax.set_xlim(0.0, 180.0)
    ax.legend(loc="upper left", fontsize=9, frameon=False)
    plt.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def write_plots(
    report: ScoreReport,
    *,
    itp: str | Path,
    cg_struct: str | Path,
    cg_traj: str | Path | Sequence[str | Path],
    out_dir: str | Path,
    bond_binwidth_nm: float = 0.001,
    angle_binwidth_deg: float = 1.0,
) -> dict[str, Path]:
    """Render bond_hists.pdf and angle_hists.pdf into ``out_dir``."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    if report.bond is not None:
        path = out_dir / "bond_hists.pdf"
        _plot_term(
            itp=itp, cg_struct=cg_struct, cg_traj=cg_traj,
            stats=report.bond, kind="bond", out_path=path,
            end_exclude=report.end_exclude,
            binwidth=bond_binwidth_nm, lo=0.0, hi=0.8,
        )
        written["bond"] = path
    if report.angle is not None:
        path = out_dir / "angle_hists.pdf"
        _plot_term(
            itp=itp, cg_struct=cg_struct, cg_traj=cg_traj,
            stats=report.angle, kind="angle", out_path=path,
            end_exclude=report.end_exclude,
            binwidth=angle_binwidth_deg, lo=0.0, hi=180.0,
        )
        written["angle"] = path
    return written


# ---------- CLI ----------


def _resolve_traj_args(traj_args: list[str]) -> list[str]:
    if len(traj_args) == 1 and any(c in traj_args[0] for c in "*?["):
        expanded = sorted(_glob.glob(traj_args[0]), key=_natural_key)
        if not expanded:
            raise SystemExit(f"glob matched no files: {traj_args[0]}")
        return expanded
    return traj_args


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Score a CG mapping by measuring bonded distributions vs .itp targets"
    )
    p.add_argument("--itp", required=True, help="canonical Polyply .itp")
    p.add_argument("--cg-struct", required=True, help="CG structure (.gro/.pdb)")
    p.add_argument(
        "--cg-traj",
        required=True,
        nargs="+",
        help="CG trajectory file(s); multiple files chain in time order, single quoted glob is natural-sorted",
    )
    p.add_argument(
        "--out-dir",
        default=None,
        help="output directory (default: same dir as --cg-struct)",
    )
    p.add_argument(
        "--end-exclude", type=int, default=2,
        help="drop terms touching the first/last N beads (default: 2)",
    )
    p.add_argument(
        "--molecule", default=None,
        help="molecule label for the report (default: derived from --itp filename)",
    )
    p.add_argument(
        "--no-plots", action="store_true",
        help="skip PDF plot generation",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    traj_paths = _resolve_traj_args(args.cg_traj)
    out_dir = Path(args.out_dir or Path(args.cg_struct).parent)
    molecule = args.molecule or Path(args.itp).stem

    report = score_mapping(
        itp=args.itp,
        cg_struct=args.cg_struct,
        cg_traj=traj_paths if len(traj_paths) > 1 else traj_paths[0],
        molecule=molecule,
        end_exclude=args.end_exclude,
    )

    json_path = out_dir / "score_report.json"
    write_json_report(report, json_path)

    print(f"scored {report.n_frames} frame(s) of {molecule} ({report.n_beads} beads)")
    if report.bond is not None:
        b = report.bond
        print(
            f"  bond  : μ_fit={b.fit_mu:.4f} nm  σ={b.fit_sigma:.4f} nm  "
            f"RMSE={b.fit_rmse:.4f}  target={b.target_mu:g}  "
            f"Δ={b.delta_vs_target:+.4f} nm  (n_obs={b.n_observations})"
        )
    if report.angle is not None:
        a = report.angle
        print(
            f"  angle : μ_fit={a.fit_mu:.2f}°   σ={a.fit_sigma:.2f}°  "
            f"RMSE={a.fit_rmse:.4f}  target={a.target_mu:g}  "
            f"Δ={a.delta_vs_target:+.2f}°  (n_obs={a.n_observations})"
        )
    print(f"  report : {json_path}")

    if not args.no_plots:
        plot_paths = write_plots(
            report,
            itp=args.itp,
            cg_struct=args.cg_struct,
            cg_traj=traj_paths if len(traj_paths) > 1 else traj_paths[0],
            out_dir=out_dir,
        )
        for kind, p in plot_paths.items():
            print(f"  {kind:<7}: {p}")


if __name__ == "__main__":
    main()
