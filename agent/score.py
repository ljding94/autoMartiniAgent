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
    """Fit a single Gaussian to (centers, density). Returns (popt, fitted_y, r2).

    Fit parameters come from least-squares (``curve_fit``). Goodness-of-fit is
    reported as the **coefficient of determination** R² = 1 − SS_res / SS_tot
    in density space, i.e. how well the fitted Gaussian traces the measured
    histogram curve. R² = 1 is a perfect fit; lower is worse and it can go
    negative for a fit worse than a flat mean line (e.g. a single Gaussian over
    a clearly bimodal distribution).

    R² is used in preference to information-theoretic scores (cross-entropy,
    KL, JS) because those are mass/overlap-based and are dominated by tail
    behaviour: a sharp, well-fit peak with rare non-Gaussian tails gets a large
    KL even though the harmonic model captures it well, while a broad bimodal
    blob that the Gaussian traces poorly gets a small KL because the fitted
    curve still covers the occupied range. R² weights every bin's density
    residual equally, so it tracks visual / physical fit quality (does the peak
    match) instead of exploding on low-density tail bins.
    """
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
    r2 = _r_squared(y, fitted)
    return popt, fitted, r2


def _r_squared(measured: np.ndarray, fitted: np.ndarray) -> float:
    """Coefficient of determination R² = 1 − SS_res / SS_tot in density space.

    ``measured`` and ``fitted`` are the histogram density and the fitted
    Gaussian evaluated at the same bin centers (both in the same density
    units, so no normalization is needed). R² = 1 is a perfect trace; it falls
    toward 0 as the fit worsens and can go negative when the Gaussian is worse
    than a flat mean line (a strongly bimodal / skewed distribution). Unlike
    KL / cross-entropy this weights each bin's residual equally, so it is not
    dominated by rare non-Gaussian tails.
    """
    measured = np.asarray(measured, dtype=float)
    fitted = np.asarray(fitted, dtype=float)
    ss_tot = float(np.sum((measured - measured.mean()) ** 2))
    if ss_tot <= 0.0:
        return float("nan")
    ss_res = float(np.sum((measured - fitted) ** 2))
    return 1.0 - ss_res / ss_tot


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
    fit_r2: float              # R² of fitted Gaussian vs measured density (1 = perfect)
    target_mu: float | None
    target_k: float | None
    delta_vs_target: float | None


@dataclass(frozen=True)
class AngleCoverage:
    """One consecutive-bond angle bead-type triple present in the topology, and
    whether the itp defines a harmonic bonded angle term for it.

    Populated on every score run (topology-only, no trajectory) so the workflow
    always surfaces bonded angles the force field omits. ``in_itp=False`` means a
    real bond-bond angle exists geometrically but has no [angles] entry.
    """
    bead_pattern: str                      # canonical bead-type triple, e.g. "Q1-SC1-Q1"
    example_triple: tuple[int, int, int]   # a representative 1-based (i, j, k)
    in_itp: bool


@dataclass(frozen=True)
class ScoreReport:
    molecule: str
    n_frames: int
    n_beads: int
    end_exclude: int
    bond_terms: list[TermStats] = field(default_factory=list)
    angle_terms: list[TermStats] = field(default_factory=list)
    angle_coverage: list[AngleCoverage] = field(default_factory=list)


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
    popt, _fitted, r2 = _fit_gaussian(centers, density)
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
        fit_r2=float(r2),
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


def _bond_pairs(topology: Topology) -> set[tuple[int, int]]:
    """Undirected set of bonded index pairs (both orderings) for O(1) lookup."""
    pairs: set[tuple[int, int]] = set()
    for b in topology.bonds:
        pairs.add((b.i, b.j))
        pairs.add((b.j, b.i))
    return pairs


def _angle_is_bonded(a: AngleTerm, bonded: set[tuple[int, int]]) -> bool:
    """True iff the angle's vertex (middle atom j) is bonded to both arms i, k.

    GROMACS ``[angles]`` entries are 3-body potentials that do NOT require the
    arms to be bonded to the vertex; some CG force fields add non-adjacent
    "structural" angle restraints (spanning several bonds) to lock a group's
    geometry. Those are not classical harmonic bond-bond angles and a single
    Gaussian is a poor model for them, so they are excluded from scoring.
    """
    return (a.j, a.i) in bonded and (a.j, a.k) in bonded


def _group_angles(topology: Topology):
    bonded = _bond_pairs(topology)
    buckets: dict[tuple, dict] = {}
    for a in topology.angles:
        if not _angle_is_bonded(a, bonded):
            continue
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


def _all_bonded_angle_groups(topology: Topology):
    """Enumerate EVERY consecutive-bond angle (vertex bonded to both arms) from
    the bonds, grouped by canonical bead-type triple — regardless of whether the
    itp defines an [angles] term for it.

    This is the "measure-only" angle set: it surfaces bonded angles the force
    field omits (e.g. sidechain 2-3-4 / 3-4-5 / 4-5-6 or backbone-bend angles)
    so their measured distributions can be inspected. Where the itp *does* define
    the angle, its (theta0, ka) target is attached so a Δ-vs-target is still
    reported; otherwise the target is None (measured-only, no Δ). Groups are keyed
    purely by bead-type pattern, so distinct angles sharing a pattern would merge.
    """
    adj: dict[int, set[int]] = {}
    for b in topology.bonds:
        adj.setdefault(b.i, set()).add(b.j)
        adj.setdefault(b.j, set()).add(b.i)

    # itp targets keyed by canonical triple (min arm, vertex, max arm)
    itp_target: dict[tuple[int, int, int], tuple[float | None, float | None]] = {}
    for a in topology.angles:
        itp_target[(min(a.i, a.k), a.j, max(a.i, a.k))] = (a.theta0, a.ka)

    # every vertex-bonded triple, canonicalized
    triples: set[tuple[int, int, int]] = set()
    for j, nbrs in adj.items():
        nl = sorted(nbrs)
        for x in range(len(nl)):
            for y in range(x + 1, len(nl)):
                triples.add((min(nl[x], nl[y]), j, max(nl[x], nl[y])))

    buckets: dict[str, dict] = {}
    for (i, j, k) in sorted(triples):
        lbl = _triple_label(
            topology.bead_type(i), topology.bead_type(j), topology.bead_type(k)
        )
        tmu, tk = itp_target.get((i, j, k), (None, None))
        v = buckets.setdefault(
            lbl, {"triples": [], "pattern": lbl, "target_mu": None, "target_k": None}
        )
        v["triples"].append((i, j, k))
        if v["target_mu"] is None and tmu is not None:  # adopt a target if any member has one
            v["target_mu"], v["target_k"] = tmu, tk

    # itp-defined groups first, then by measured-agnostic target/label
    ordered = sorted(
        buckets.items(),
        key=lambda kv: (kv[1]["target_mu"] is None, kv[1]["target_mu"] or 0.0, kv[0]),
    )
    return [
        ((v["target_mu"], v["target_k"]), v["triples"], lbl, v["pattern"], v["target_mu"], v["target_k"])
        for lbl, v in ordered
    ]


def _angle_coverage(topology: Topology) -> list[AngleCoverage]:
    """Topology-only audit: every consecutive-bond angle bead-type triple and
    whether the itp defines a harmonic bonded angle for it. No trajectory needed."""
    return [
        AngleCoverage(
            bead_pattern=pattern,
            example_triple=tuple(triples[0]),
            in_itp=target_mu is not None,
        )
        for (_key, triples, _label, pattern, target_mu, _tk) in _all_bonded_angle_groups(topology)
    ]


def score_mapping(
    *,
    itp: str | Path,
    cg_struct: str | Path,
    cg_traj: str | Path | Sequence[str | Path],
    molecule: str | None = None,
    end_exclude: int = 0,
    bond_binwidth_nm: float = 0.001,
    angle_binwidth_deg: float = 1.0,
    all_bonded_angles: bool = False,
) -> ScoreReport:
    """Score a CG mapping by grouping bonded terms per (target μ, target k) and
    reporting per-group Gaussian-fit statistics vs the ``.itp`` target.

    By default only angles the itp defines (and whose vertex is bonded to both
    arms) are scored. With ``all_bonded_angles=True`` the scorer instead measures
    EVERY consecutive-bond angle in the topology — including bonded angles the itp
    omits — reporting their fitted μ/σ/R² and a Δ-vs-target only where the itp
    supplies one. Useful for spotting bonded angles missing from the force field.
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

    angle_groups = _all_bonded_angle_groups(topology) if all_bonded_angles else _group_angles(topology)
    angle_terms: list[TermStats] = []
    for group_id, (_key, trips_1, label, pattern, target_mu, target_k) in enumerate(angle_groups):
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
        angle_coverage=_angle_coverage(topology),
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
    all_bonded_angles: bool = False,
) -> None:
    """One PDF, one subplot per group, with measured + Gaussian fit + target.

    ``all_bonded_angles`` must match how the report's angle terms were grouped so
    that group ids and triples line up (else panels go missing / mislabelled).
    """
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
        groups = _all_bonded_angle_groups(topology) if all_bonded_angles else _group_angles(topology)
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
            f"R² = {stats.fit_r2:.3f}",
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
    name_suffix: str = "",
    all_bonded_angles: bool = False,
) -> dict[str, Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    if report.bond_terms:
        path = out_dir / f"bond_hists{name_suffix}.pdf"
        _plot_groups(
            itp=itp, cg_struct=cg_struct, cg_traj=cg_traj,
            terms=report.bond_terms, kind="bond", out_path=path,
            end_exclude=report.end_exclude,
            binwidth=bond_binwidth_nm, lo=0.0, hi=0.8,
        )
        written["bond"] = path
    if report.angle_terms:
        path = out_dir / f"angle_hists{name_suffix}.pdf"
        _plot_groups(
            itp=itp, cg_struct=cg_struct, cg_traj=cg_traj,
            terms=report.angle_terms, kind="angle", out_path=path,
            end_exclude=report.end_exclude,
            binwidth=angle_binwidth_deg, lo=0.0, hi=180.0,
            all_bonded_angles=all_bonded_angles,
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
        "--end-exclude", type=int, default=0,
        help="drop terms touching the first/last N *beads* (default: 0 = keep all "
             "monomers). Note this counts beads, not monomers: for a molecule with "
             "K beads per monomer, use a multiple of K to trim whole end monomers.",
    )
    p.add_argument(
        "--molecule", default=None,
        help="molecule label for the report (default: derived from --itp filename)",
    )
    p.add_argument(
        "--no-plots", action="store_true",
        help="skip PDF plot generation",
    )
    p.add_argument(
        "--all-bonded-angles", action="store_true",
        help="measure EVERY consecutive-bond angle, including those the itp omits "
             "(Δ-vs-target only where the itp defines one)",
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
        all_bonded_angles=args.all_bonded_angles,
    )

    suffix = "_all_angles" if args.all_bonded_angles else ""
    json_path = out_dir / f"score_report{suffix}.json"
    write_json_report(report, json_path)

    print(f"scored {report.n_frames} frame(s) of {molecule} ({report.n_beads} beads, end_exclude={report.end_exclude})")
    if args.all_bonded_angles:
        print("  [all-bonded-angles mode] angles = every consecutive-bond angle; Δ shown only where the itp defines a target")
    print(f"  bond groups  ({len(report.bond_terms)}):")
    for t in report.bond_terms:
        tgt = f"target={t.target_mu:g}" if t.target_mu is not None else "target=—"
        d = f"Δ={t.delta_vs_target:+.4f} nm" if t.delta_vs_target is not None else "Δ=—"
        print(
            f"    {t.label:<22} μ_fit={t.fit_mu:.4f} σ={t.fit_sigma:.4f}  "
            f"R²={t.fit_r2:.3f}  {tgt}  {d}  (n_obs={t.n_observations})"
        )
    print(f"  angle groups ({len(report.angle_terms)}):")
    for t in report.angle_terms:
        tgt = f"target={t.target_mu:g}" if t.target_mu is not None else "target=—"
        d = f"Δ={t.delta_vs_target:+.2f}°" if t.delta_vs_target is not None else "Δ=—"
        print(
            f"    {t.label:<22} μ_fit={t.fit_mu:.2f}° σ={t.fit_sigma:.2f}°  "
            f"R²={t.fit_r2:.3f}  {tgt}  {d}  (n_obs={t.n_observations})"
        )
    cov = report.angle_coverage
    missing = [c.bead_pattern for c in cov if not c.in_itp]
    print(
        f"  angle coverage: {len(cov)} bonded-angle type(s) in topology; "
        f"itp defines {len(cov) - len(missing)}, missing {len(missing)}"
    )
    if missing:
        print(f"    missing from itp: {', '.join(missing)}")
        if not args.all_bonded_angles:
            print("    → re-run with --all-bonded-angles to measure their distributions")
    print(f"  report : {json_path}")

    if not args.no_plots:
        plot_paths = write_plots(
            report,
            itp=args.itp,
            cg_struct=args.cg_struct,
            cg_traj=traj_paths if len(traj_paths) > 1 else traj_paths[0],
            out_dir=out_dir,
            name_suffix=suffix,
            all_bonded_angles=args.all_bonded_angles,
        )
        for kind, p in plot_paths.items():
            print(f"  {kind:<7}: {p}")


if __name__ == "__main__":
    main()
