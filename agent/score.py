"""CG mapping quality scorer.

Reads a CG trajectory + the canonical Polyply (or auto-martini) ``.itp``
for the molecule, measures bonded distance and bond-angle distributions
over all frames, fits a Gaussian per **term group** (bonds grouped by
``(b0, kb)``; angles by ``(theta0, ka)``), and emits a score report
comparing each group's measured mean to the ``.itp`` target.

Torsions are intentionally not scored — Martini polymer torsions are
typically weak or non-Gaussian.

CLI::

  python -m agent.score \\
    --itp        reference/polyply_PEO20/PEO20.itp \\
    --cg-struct  derived/PEO20_solu/PEO20_cg.gro \\
    --cg-traj    derived/PEO20_solu/PEO20_cg.xtc \\
    --out-dir    derived/PEO20_solu

Outputs (in ``--out-dir``):
  - ``score_report.json``  : per-group stats (measured μ/σ, fit μ/σ,
                             fit RMSE, target μ, Δ vs target)
  - ``bond_hists.pdf``     : one panel per bond group
  - ``angle_hists.pdf``    : one panel per angle group
"""

from __future__ import annotations

import argparse
import glob as _glob
import json
import re
from dataclasses import asdict, dataclass, field
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
class AtomEntry:
    idx: int         # 1-based bead index
    type: str        # Martini bead type, e.g. "SN3r"
    resnr: int | None
    resname: str | None
    name: str | None
    charge: float | None
    mass: float | None


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
    atoms: list[AtomEntry]
    bonds: list[BondTerm]
    angles: list[AngleTerm]

    def bead_type(self, idx1: int) -> str:
        """1-based lookup of Martini bead type; falls back to '?' if unknown."""
        if 1 <= idx1 <= len(self.atoms):
            return self.atoms[idx1 - 1].type
        return "?"


def load_topology(itp_path: str | Path) -> Topology:
    """Parse [ atoms ], [ bonds ], [ angles ] from a GROMACS-style .itp file."""
    atoms: list[AtomEntry] = []
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
            if section == "atoms" and len(toks) >= 2:
                idx = int(toks[0])
                bead_type = toks[1]
                resnr = int(toks[2]) if len(toks) > 2 else None
                resname = toks[3] if len(toks) > 3 else None
                name = toks[4] if len(toks) > 4 else None
                charge = float(toks[6]) if len(toks) > 6 else None
                mass = float(toks[7]) if len(toks) > 7 else None
                atoms.append(AtomEntry(
                    idx=idx, type=bead_type, resnr=resnr, resname=resname,
                    name=name, charge=charge, mass=mass,
                ))
            elif section == "bonds" and len(toks) >= 2:
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
    return Topology(atoms=atoms, bonds=bonds, angles=angles)


# ---------- grouping helpers ----------


def _bond_group_key(b: BondTerm) -> tuple[float | None, float | None]:
    """Group bonds by physical target: (b0, kb). Rounded to catch FP noise."""
    b0 = round(b.b0, 4) if b.b0 is not None else None
    kb = round(b.kb, 1) if b.kb is not None else None
    return (b0, kb)


def _angle_group_key(a: AngleTerm) -> tuple[float | None, float | None]:
    t0 = round(a.theta0, 2) if a.theta0 is not None else None
    ka = round(a.ka, 2) if a.ka is not None else None
    return (t0, ka)


def _pair_label(ti: str, tj: str) -> str:
    """Order-invariant bead-type pair, hyphenated."""
    a, b = sorted([ti, tj])
    return f"{a}-{b}"


def _triple_label(ti: str, tj: str, tk: str) -> str:
    """Reverse-invariant triple label; middle bead is anchored."""
    a, c = sorted([ti, tk])
    return f"{a}-{tj}-{c}"


# ---------- core scoring ----------


@dataclass(frozen=True)
class TermStats:
    """Aggregated stats for one bonded-term group (one (b0,kb) or (θ0,ka) set)."""
    group_id: int              # index within its kind (bond or angle)
    label: str                 # human-readable, e.g. "SC1-TN5a"
    bead_pattern: str          # bead-type pair or triple
    n_members: int             # number of terms of this type in the topology
    n_observations: int        # total measurements across trajectory
    measured_mu: float
    measured_sigma: float
    fit_mu: float
    fit_sigma: float
    fit_amp: float
    fit_rmse: float
    target_mu: float | None
    target_k: float | None
    delta_vs_target: float | None


@dataclass(frozen=True)
class ScoreReport:
    molecule: str
    n_frames: int
    n_beads: int
    end_exclude: int
    bond_terms: list[TermStats] = field(default_factory=list)
    angle_terms: list[TermStats] = field(default_factory=list)


def _drop_end_excluded(tuples, n_beads: int, end_exclude: int):
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
    *,
    group_id: int,
    label: str,
    bead_pattern: str,
    values: np.ndarray,
    binwidth: float,
    lo: float,
    hi: float,
    target_mu: float | None,
    target_k: float | None,
    n_members: int,
) -> TermStats:
    centers, density = _hist(values, binwidth, lo, hi)
    measured_mu, measured_sigma = _moments(centers, density)
    popt, _fitted, rmse = _fit_gaussian(centers, density)
    delta = float(measured_mu - target_mu) if target_mu is not None else None
    return TermStats(
        group_id=group_id,
        label=label,
        bead_pattern=bead_pattern,
        n_members=n_members,
        n_observations=int(values.size),
        measured_mu=float(measured_mu),
        measured_sigma=float(measured_sigma),
        fit_mu=float(popt[0]),
        fit_sigma=float(abs(popt[1])),
        fit_amp=float(popt[2]),
        fit_rmse=float(rmse),
        target_mu=float(target_mu) if target_mu is not None else None,
        target_k=float(target_k) if target_k is not None else None,
        delta_vs_target=delta,
    )


def _group_bonds(topology: Topology):
    """Return list of (group_key, [(i,j)...], label, target_mu, target_k) sorted by target_mu.

    ``label`` disambiguates when multiple groups share the same bead-type pair.
    """
    buckets: dict[tuple, dict] = {}
    for b in topology.bonds:
        key = _bond_group_key(b)
        pair = _pair_label(topology.bead_type(b.i), topology.bead_type(b.j))
        buckets.setdefault(key, {"pairs": [], "pattern": pair, "target_mu": b.b0, "target_k": b.kb})
        buckets[key]["pairs"].append((b.i, b.j))
        if buckets[key]["pattern"] != pair:
            buckets[key]["pattern"] = pair  # keep last; label collisions are handled below
    ordered = sorted(buckets.items(), key=lambda kv: (kv[0][0] if kv[0][0] is not None else -1.0))
    pattern_counts: dict[str, int] = {}
    for _, v in ordered:
        pattern_counts[v["pattern"]] = pattern_counts.get(v["pattern"], 0) + 1
    seen: dict[str, int] = {}
    out = []
    for key, v in ordered:
        p = v["pattern"]
        if pattern_counts[p] > 1:
            seen[p] = seen.get(p, 0) + 1
            label = f"{p} #{seen[p]}"
        else:
            label = p
        out.append((key, v["pairs"], label, p, v["target_mu"], v["target_k"]))
    return out


def _group_angles(topology: Topology):
    buckets: dict[tuple, dict] = {}
    for a in topology.angles:
        key = _angle_group_key(a)
        triple = _triple_label(
            topology.bead_type(a.i),
            topology.bead_type(a.j),
            topology.bead_type(a.k),
        )
        buckets.setdefault(key, {"triples": [], "pattern": triple, "target_mu": a.theta0, "target_k": a.ka})
        buckets[key]["triples"].append((a.i, a.j, a.k))
    ordered = sorted(buckets.items(), key=lambda kv: (kv[0][0] if kv[0][0] is not None else -1.0))
    pattern_counts: dict[str, int] = {}
    for _, v in ordered:
        pattern_counts[v["pattern"]] = pattern_counts.get(v["pattern"], 0) + 1
    seen: dict[str, int] = {}
    out = []
    for key, v in ordered:
        p = v["pattern"]
        if pattern_counts[p] > 1:
            seen[p] = seen.get(p, 0) + 1
            label = f"{p} #{seen[p]}"
        else:
            label = p
        out.append((key, v["triples"], label, p, v["target_mu"], v["target_k"]))
    return out


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
    """Score a CG mapping by grouping bonded terms per (target μ, target k) and
    reporting per-group Gaussian-fit statistics vs the ``.itp`` target.
    """
    topology = load_topology(itp)

    if isinstance(cg_traj, (str, Path)):
        traj_arg: str | list[str] = str(cg_traj)
    else:
        traj_arg = [str(p) for p in cg_traj]
    traj = md.load(traj_arg, top=str(cg_struct))
    n_beads = traj.n_atoms
    n_frames = traj.n_frames

    bond_terms: list[TermStats] = []
    for group_id, (_key, pairs_1, label, pattern, target_mu, target_k) in enumerate(_group_bonds(topology)):
        pairs_0 = [(i - 1, j - 1) for (i, j) in pairs_1]
        kept = _drop_end_excluded(pairs_0, n_beads, end_exclude)
        if not kept:
            continue
        distances_nm = md.compute_distances(traj, kept, opt=True).ravel()
        bond_terms.append(_stats_from_values(
            group_id=group_id,
            label=label,
            bead_pattern=pattern,
            values=distances_nm,
            binwidth=bond_binwidth_nm,
            lo=0.0,
            hi=0.8,
            target_mu=target_mu,
            target_k=target_k,
            n_members=len(kept),
        ))

    angle_terms: list[TermStats] = []
    for group_id, (_key, trips_1, label, pattern, target_mu, target_k) in enumerate(_group_angles(topology)):
        trips_0 = [(i - 1, j - 1, k - 1) for (i, j, k) in trips_1]
        kept = _drop_end_excluded(trips_0, n_beads, end_exclude)
        if not kept:
            continue
        angles_deg = (180.0 / np.pi) * md.compute_angles(traj, kept, opt=True).ravel()
        angle_terms.append(_stats_from_values(
            group_id=group_id,
            label=label,
            bead_pattern=pattern,
            values=angles_deg,
            binwidth=angle_binwidth_deg,
            lo=0.0,
            hi=180.0,
            target_mu=target_mu,
            target_k=target_k,
            n_members=len(kept),
        ))

    return ScoreReport(
        molecule=molecule or "MOL",
        n_frames=n_frames,
        n_beads=n_beads,
        end_exclude=end_exclude,
        bond_terms=bond_terms,
        angle_terms=angle_terms,
    )


# ---------- output writers ----------


def write_json_report(report: ScoreReport, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(asdict(report), indent=2))


def _plot_groups(
    *,
    itp: str | Path,
    cg_struct: str | Path,
    cg_traj: str | Path | Sequence[str | Path],
    terms: list[TermStats],
    kind: str,
    out_path: Path,
    end_exclude: int,
    binwidth: float,
    lo: float,
    hi: float,
) -> None:
    """One PDF, one subplot per group, with measured + Gaussian fit + target."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not terms:
        return

    topology = load_topology(itp)
    if isinstance(cg_traj, (str, Path)):
        traj_arg: str | list[str] = str(cg_traj)
    else:
        traj_arg = [str(p) for p in cg_traj]
    traj = md.load(traj_arg, top=str(cg_struct))
    n_beads = traj.n_atoms

    if kind == "bond":
        groups = _group_bonds(topology)
        get_indices = lambda ijk: [(i - 1, j - 1) for (i, j) in ijk]
        measure = lambda kept: md.compute_distances(traj, kept, opt=True).ravel()
        xlabel = "bond distance (nm)"
        xlim = (0.20, 0.50)
        unit = "nm"
        val_fmt = ".4f"
    elif kind == "angle":
        groups = _group_angles(topology)
        get_indices = lambda ijk: [(i - 1, j - 1, k - 1) for (i, j, k) in ijk]
        measure = lambda kept: (180.0 / np.pi) * md.compute_angles(traj, kept, opt=True).ravel()
        xlabel = "bond angle (deg)"
        xlim = (0.0, 180.0)
        unit = "°"
        val_fmt = ".2f"
    else:
        raise ValueError(f"unknown kind {kind!r}")

    stats_by_gid = {t.group_id: t for t in terms}
    plotted = [(gid, g) for (gid, g) in enumerate(groups) if gid in stats_by_gid]
    n = len(plotted)
    ncol = 2 if n > 1 else 1
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.5 * ncol, 3.8 * nrow), squeeze=False)

    for i, (gid, (_key, ijk_1, label, _pattern, _target_mu, _target_k)) in enumerate(plotted):
        ax = axes[i // ncol][i % ncol]
        stats = stats_by_gid[gid]
        kept = _drop_end_excluded(get_indices(ijk_1), n_beads, end_exclude)
        values = measure(kept)
        centers, density = _hist(values, binwidth, lo, hi)
        fitted = _gaussian(centers, stats.fit_mu, stats.fit_sigma, stats.fit_amp)
        ax.plot(centers, density, lw=1.4, label="measured")
        ax.plot(centers, fitted, "r--", lw=1.2, label="Gaussian fit")
        if stats.target_mu is not None:
            ax.axvline(
                stats.target_mu, color="green", ls=":", lw=1.2,
                label=f"target = {stats.target_mu:g}",
            )
        annot = [
            f"μ = {format(stats.fit_mu, val_fmt)} {unit}",
            f"σ = {format(stats.fit_sigma, val_fmt)} {unit}",
            f"RMSE = {stats.fit_rmse:.4f}",
        ]
        if stats.delta_vs_target is not None:
            delta_fmt = "+" + val_fmt
            annot.append(f"Δ = {format(stats.delta_vs_target, delta_fmt)} {unit}")
        ax.text(
            0.97, 0.95, "\n".join(annot),
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="0.7"),
        )
        ax.set_title(f"{label}  (n_obs={stats.n_observations})", fontsize=10)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_ylabel("probability density", fontsize=9)
        ax.set_xlim(xlim)
        ax.legend(loc="upper left", fontsize=7, frameon=False)

    for j in range(n, nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")

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
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    if report.bond_terms:
        path = out_dir / "bond_hists.pdf"
        _plot_groups(
            itp=itp, cg_struct=cg_struct, cg_traj=cg_traj,
            terms=report.bond_terms, kind="bond", out_path=path,
            end_exclude=report.end_exclude,
            binwidth=bond_binwidth_nm, lo=0.0, hi=0.8,
        )
        written["bond"] = path
    if report.angle_terms:
        path = out_dir / "angle_hists.pdf"
        _plot_groups(
            itp=itp, cg_struct=cg_struct, cg_traj=cg_traj,
            terms=report.angle_terms, kind="angle", out_path=path,
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
    p.add_argument("--itp", required=True, help="canonical CG .itp")
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

    print(f"scored {report.n_frames} frame(s) of {molecule} ({report.n_beads} beads, end_exclude={report.end_exclude})")
    print(f"  bond groups  ({len(report.bond_terms)}):")
    for t in report.bond_terms:
        tgt = f"target={t.target_mu:g}" if t.target_mu is not None else "target=—"
        d = f"Δ={t.delta_vs_target:+.4f} nm" if t.delta_vs_target is not None else "Δ=—"
        print(
            f"    {t.label:<22} μ_fit={t.fit_mu:.4f} σ={t.fit_sigma:.4f}  "
            f"RMSE={t.fit_rmse:.4f}  {tgt}  {d}  (n_obs={t.n_observations})"
        )
    print(f"  angle groups ({len(report.angle_terms)}):")
    for t in report.angle_terms:
        tgt = f"target={t.target_mu:g}" if t.target_mu is not None else "target=—"
        d = f"Δ={t.delta_vs_target:+.2f}°" if t.delta_vs_target is not None else "Δ=—"
        print(
            f"    {t.label:<22} μ_fit={t.fit_mu:.2f}° σ={t.fit_sigma:.2f}°  "
            f"RMSE={t.fit_rmse:.4f}  {tgt}  {d}  (n_obs={t.n_observations})"
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
