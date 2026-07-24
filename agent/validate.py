"""Validation harness — does the Gaussianity objective track mapping quality?

Two sim-free experiments against a good reference mapping (the paper's core
validation, using Martini library / accepted mappings as a relative ground truth):

  - **rank**: perturb the good mapping in chemically-wrong ways of increasing
    severity and show the objective monotonically worsens → the objective ORDERS
    mappings the way a chemist would (a valid quality signal, no simulation).
  - **recover** (see ``recover`` below): degrade the good mapping, run the agent
    loop, and show it climbs back toward the good objective.

Everything reuses what we already built — ``repair.py`` verbs (via
``agent.loop.apply_ops``) to perturb, and ``agent.evaluate`` to score — so this is
an experiment harness, not new core machinery. Perturbations keep the **bead count
fixed** (reassignments only) so the mean-``(1-R^2)`` objective stays directly
comparable across the ladder.

CLI::

  python -m agent.validate rank \\
    --mapping reference/PSBMA_20mer/PSBMA20_mapping.json \\
    --itp     reference/PSBMA_20mer/PSBMA20.itp \\
    --aa-top  reference/PSBMA_20mer/PSBMA_20mer_no_water.gro \\
    --aa-traj reference/PSBMA_20mer/PSBMA_20mer_no_water_skip10.xtc \\
    --cg-struct reference/PSBMA_20mer/PSBMA_20mer_no_water.gro \\
    --frame-stride 25
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import dataclass, field
from pathlib import Path

from agent import loop as _loop
from agent import repair as _repair


# ---------- perturbations (chemically-wrong, fixed bead count) ----------

# Named PSBMA perturbations: each splits a functional group across a bonded bead
# boundary or mis-assigns a backbone atom — all things a good objective should
# penalise. Roles are within one monomer (tiled across all 20).
NAMED_PERTURBATIONS_PSBMA: dict[str, list[dict]] = {
    "split_sulfonate":  [{"op": "reassign", "atom_names": ["O3"], "from_role": 6, "to_role": 5}],
    "split_ammonium":   [{"op": "reassign", "atom_names": ["C7"], "from_role": 4, "to_role": 5}],
    "split_ester_OCH2": [{"op": "reassign", "atom_names": ["C5"], "from_role": 3, "to_role": 2}],
    "mislabel_backbone":[{"op": "reassign", "atom_names": ["C2"], "from_role": 1, "to_role": 2}],
}


def _bonded_role_adjacency(state: _repair.MappingState, bpm: int) -> dict[int, set[int]]:
    adj: dict[int, set[int]] = {}
    for b in state.bonds:
        if b.i <= bpm and b.j <= bpm:
            adj.setdefault(b.i, set()).add(b.j)
            adj.setdefault(b.j, set()).add(b.i)
    return adj


def scramble(state: _repair.MappingState, n_moves: int, seed: int,
             ref_positions=None) -> _repair.MappingState:
    """Apply ``n_moves`` random valid role-based reassignments (tiled across all
    monomers). Each move picks a random heavy atom in a random role and sends it
    (with its hydrogens) to a random bonded-neighbour role; moves that would empty
    a bead of heavy atoms or break adjacency are retried. Deterministic per seed."""
    rng = random.Random(seed)
    _, bpm = _loop._beads_per_monomer(state)
    adj = _bonded_role_adjacency(state, bpm)
    cur = state
    moves = attempts = 0
    while moves < n_moves and attempts < 60 * max(1, n_moves):
        attempts += 1
        fr = rng.randint(1, bpm)
        nbrs = list(adj.get(fr, ()))
        if not nbrs:
            continue
        bead = cur.bead(fr)
        heavy = [nm for i, nm in zip(bead["atom_indices"], bead["atom_names"])
                 if i in cur.heavy_atoms]
        if len(heavy) <= 1:                       # don't strip a bead's last heavy atom
            continue
        op = {"op": "reassign", "atom_names": [rng.choice(heavy)],
              "from_role": fr, "to_role": rng.choice(nbrs)}
        try:
            cur = _loop.apply_ops(cur, [op], ref_positions=ref_positions)
            moves += 1
        except Exception:
            continue
    return cur


# ---------- rank experiment ----------


@dataclass
class RankRow:
    label: str
    severity: int | None      # None for named perturbations; move-count for scrambles
    objective: float
    delta: float              # objective - baseline


@dataclass
class RankResult:
    baseline: float
    rows: list[RankRow] = field(default_factory=list)


def run_rank(
    state: _repair.MappingState,
    evaluate_fn,
    *,
    ref_positions=None,
    severities: tuple[int, ...] = (1, 2, 4, 8),
    seeds: tuple[int, ...] = tuple(range(8)),
    named: dict[str, list[dict]] | None = NAMED_PERTURBATIONS_PSBMA,
) -> RankResult:
    """Score the baseline, each named perturbation, and a random-scramble ladder."""
    base = float(evaluate_fn(state).objective)
    rows: list[RankRow] = []
    for label, ops in (named or {}).items():
        s = _loop.apply_ops(state, ops, ref_positions=ref_positions)
        obj = float(evaluate_fn(s).objective)
        rows.append(RankRow(label, None, obj, obj - base))
    for sev in severities:
        for seed in seeds:
            s = scramble(state, sev, seed, ref_positions=ref_positions)
            obj = float(evaluate_fn(s).objective)
            rows.append(RankRow(f"scramble×{sev}", sev, obj, obj - base))
    return RankResult(baseline=base, rows=rows)


def plot_rank(result: RankResult, out_path: str | Path,
              severities: tuple[int, ...] = (1, 2, 4, 8)) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))

    # left: scramble ladder (objective vs # random moves), mean +/- range
    xs, means, los, his = [], [], [], []
    for sev in severities:
        objs = [r.objective for r in result.rows if r.severity == sev]
        if not objs:
            continue
        xs.append(sev)
        means.append(sum(objs) / len(objs)); los.append(min(objs)); his.append(max(objs))
    ax1.axhline(result.baseline, color="tab:green", ls="--", lw=1.4,
                label=f"good mapping = {result.baseline:.4f}")
    if xs:
        ax1.fill_between(xs, los, his, color="tab:red", alpha=0.15)
        ax1.plot(xs, means, "o-", color="tab:red", label="random scramble (mean, range)")
    ax1.set_xlabel("perturbation severity (# random atom moves)")
    ax1.set_ylabel("Gaussianity error  (mean 1 − R²)")
    ax1.set_title("Objective degrades with scrambling")
    ax1.legend(fontsize=8, frameon=False)

    # right: named chemically-wrong perturbations vs baseline
    named = [r for r in result.rows if r.severity is None]
    labels = ["good\nmapping"] + [r.label for r in named]
    vals = [result.baseline] + [r.objective for r in named]
    colors = ["tab:green"] + ["tab:red"] * len(named)
    ax2.bar(range(len(vals)), vals, color=colors, alpha=0.8)
    ax2.axhline(result.baseline, color="tab:green", ls="--", lw=1.0)
    ax2.set_xticks(range(len(vals)))
    ax2.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax2.set_ylabel("Gaussianity error  (mean 1 − R²)")
    ax2.set_title("Chemically-wrong edits score worse")

    fig.suptitle("Rank validation — the objective orders mappings by quality", fontsize=12)
    plt.tight_layout(rect=(0, 0, 1, 0.96))
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


# ---------- recover experiment ----------


def _atom_roles(state: _repair.MappingState) -> dict[int, int]:
    """Map each AA atom index to its bead's role (1-based position within a monomer)."""
    _, bpm = _loop._beads_per_monomer(state)
    roles: dict[int, int] = {}
    for b in state.beads:
        role = ((b["bead_id"] - 1) % bpm) + 1
        for a in b["atom_indices"]:
            roles[a] = role
    return roles


def _grouping_match(good_roles: dict[int, int], other_roles: dict[int, int]) -> float:
    """Fraction of atoms assigned to the same role as in the good mapping (1 = identical grouping)."""
    if not good_roles:
        return float("nan")
    same = sum(1 for a, r in good_roles.items() if other_roles.get(a) == r)
    return same / len(good_roles)


@dataclass
class RecoverRow:
    label: str
    good: float
    degraded: float
    recovered: float
    recovery_fraction: float   # (degraded - recovered) / (degraded - good); 1 = back to good
    grouping_match: float      # fraction of atoms back in their good-mapping role


def run_recover(
    state: _repair.MappingState,
    make_policy,
    evaluate_fn,
    *,
    ref_positions=None,
    perturbations: dict[str, list[dict]] | None = NAMED_PERTURBATIONS_PSBMA,
    restarts: int = 2,
    max_iters: int = 8,
    plateau_k: int = 4,
    on_step=None,
) -> list[RecoverRow]:
    """For each perturbation: degrade the good mapping, run the loop, and measure how
    far it climbs back (objective recovery + how much of the original grouping it
    reconstructs)."""
    good = float(evaluate_fn(state).objective)
    good_roles = _atom_roles(state)
    rows: list[RecoverRow] = []
    for name, ops in (perturbations or {}).items():
        degraded = _loop.apply_ops(state, ops, ref_positions=ref_positions)
        deg = float(evaluate_fn(degraded).objective)
        best, result = _loop.run_loop_restarts(
            degraded, make_policy, evaluate_fn, restarts=restarts,
            ref_positions=ref_positions, max_iters=max_iters, plateau_k=plateau_k,
        )
        rec = float(result.best_objective)
        frac = (deg - rec) / (deg - good) if deg > good else float("nan")
        row = RecoverRow(name, good, deg, rec, frac, _grouping_match(good_roles, _atom_roles(best)))
        rows.append(row)
        if on_step:
            on_step(row)
    return rows


def plot_recover(rows: list[RecoverRow], good: float, out_path: str | Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(1.6 * len(rows) + 3, 4.8))
    x = range(len(rows))
    ax.axhline(good, color="tab:green", ls="--", lw=1.4, label=f"good mapping = {good:.4f}")
    ax.bar([i - 0.2 for i in x], [r.degraded for r in rows], width=0.38,
           color="tab:red", alpha=0.8, label="degraded")
    ax.bar([i + 0.2 for i in x], [r.recovered for r in rows], width=0.38,
           color="tab:blue", alpha=0.8, label="recovered by loop")
    for i, r in enumerate(rows):
        ax.annotate(f"{100*r.recovery_fraction:.0f}%\nrecovered\n(grp {100*r.grouping_match:.0f}%)",
                    (i, r.recovered), textcoords="offset points", xytext=(0, 6),
                    ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(list(x))
    ax.set_xticklabels([r.label for r in rows], rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("Gaussianity error  (mean 1 − R²)")
    ax.set_title("Recover validation — the agent climbs damaged mappings back toward good")
    ax.legend(fontsize=8, frameon=False)
    plt.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


# ---------- CLI ----------


def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--mapping", required=True)
    p.add_argument("--itp", required=True)
    p.add_argument("--aa-top", required=True)
    p.add_argument("--aa-traj", required=True, nargs="+")
    p.add_argument("--cg-struct", default=None, help="AA .gro for hydrogen auto-carry")
    p.add_argument("--frame-stride", type=int, default=25)
    p.add_argument("--work-root", default=None)
    p.add_argument("--out-dir", default=None)


def _evaluate_fn(args, state):
    from agent.evaluate import evaluate_state
    mol = state.mapping.get("molecule", "MOL")
    work_root = Path(args.work_root or f"derived/{mol}/validate")
    def fn(s):
        return evaluate_state(s, aa_top=args.aa_top, aa_traj=args.aa_traj,
                              work_root=work_root, molecule=mol,
                              frame_stride=args.frame_stride)
    return fn, mol


def _cmd_rank(args) -> None:
    state = _repair.load_state(args.mapping, args.itp)
    ref = _repair.load_ref_positions(args.cg_struct) if args.cg_struct else None
    evaluate_fn, mol = _evaluate_fn(args, state)
    result = run_rank(state, evaluate_fn, ref_positions=ref)

    out_dir = Path(args.out_dir or f"derived/{mol}")
    print(f"rank validation for {mol}  (baseline objective = {result.baseline:.4f})\n")
    print(f"  {'perturbation':<20} {'objective':>10} {'Δ vs good':>10}")
    for r in sorted(result.rows, key=lambda r: r.objective):
        print(f"  {r.label:<20} {r.objective:>10.4f} {r.delta:>+10.4f}")
    worse = sum(1 for r in result.rows if r.delta > 0)
    print(f"\n  {worse}/{len(result.rows)} perturbations scored WORSE than the good mapping")

    from dataclasses import asdict
    (out_dir).mkdir(parents=True, exist_ok=True)
    (out_dir / "rank_report.json").write_text(json.dumps(asdict(result), indent=2))
    plot = plot_rank(result, out_dir / "rank_validation.pdf")
    print(f"  report: {out_dir / 'rank_report.json'}")
    print(f"  plot  : {plot}")


def _add_llm_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--model", default="anthropic/claude-opus-4.8",
                   help="OpenRouter/OpenAI-compatible model slug")
    p.add_argument("--base-url", default=None)
    p.add_argument("--api-key-env", default="OPENROUTER_API_KEY")
    p.add_argument("--temperature", type=float, default=0.6)
    p.add_argument("--restarts", type=int, default=2)
    p.add_argument("--max-iters", type=int, default=8)
    p.add_argument("--plateau-k", type=int, default=4)


def _cmd_recover(args) -> None:
    from agent.loop import LLMPolicy
    state = _repair.load_state(args.mapping, args.itp)
    ref = _repair.load_ref_positions(args.cg_struct) if args.cg_struct else None
    evaluate_fn, mol = _evaluate_fn(args, state)

    def make_policy():
        return LLMPolicy(model=args.model, base_url=args.base_url,
                         api_key_env=args.api_key_env, temperature=args.temperature)

    print(f"recover validation for {mol}  (model={args.model}, restarts={args.restarts})\n")

    def on_step(r: RecoverRow) -> None:
        print(f"  {r.label:<20} degraded {r.degraded:.4f} -> recovered {r.recovered:.4f}  "
              f"({100*r.recovery_fraction:.0f}% of the gap closed, grouping {100*r.grouping_match:.0f}%)")

    rows = run_recover(state, make_policy, evaluate_fn, ref_positions=ref,
                       restarts=args.restarts, max_iters=args.max_iters,
                       plateau_k=args.plateau_k, on_step=on_step)
    good = rows[0].good if rows else float("nan")
    from dataclasses import asdict
    out_dir = Path(args.out_dir or f"derived/{mol}")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "recover_report.json").write_text(
        json.dumps({"good": good, "rows": [asdict(r) for r in rows]}, indent=2))
    plot = plot_recover(rows, good, out_dir / "recover_validation.pdf")
    mean_frac = sum(r.recovery_fraction for r in rows) / len(rows) if rows else float("nan")
    print(f"\n  mean recovery: {100*mean_frac:.0f}% of the degradation gap closed")
    print(f"  report: {out_dir / 'recover_report.json'}")
    print(f"  plot  : {plot}")


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Validation harness for the mapping objective")
    sub = p.add_subparsers(dest="cmd", required=True)

    rank = sub.add_parser("rank", help="perturb a good mapping, show the objective degrades")
    _add_common(rank)
    rank.set_defaults(func=_cmd_rank)

    recover = sub.add_parser("recover", help="degrade a good mapping, run the loop, show it climbs back")
    _add_common(recover)
    _add_llm_args(recover)
    recover.set_defaults(func=_cmd_recover)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
