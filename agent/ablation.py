"""Agentic-vs-deterministic ablation — the LLM's chemistry keeps the search honest.

A blind optimizer of the Gaussianity objective can *game* it with chemically-
invalid edits. On PSBMA, steepest-ascent single-move greedy drives the objective to
**0.0217** — lower than the good mapping — by **splitting the sulfonate** (S in one
bead, its three O's in another): lower error, but a −1 charge divorced from its
sulfur. The objective alone does not forbid it.

So the ablation compares three strategies over the same action space (repair verbs)
+ same objective, and tracks **chemical validity** (functional-group integrity):

  - **greedy (objective-only)**   — free to game it → low error but INVALID;
  - **greedy (chem-constrained)**  — FG-integrity enforced → the best VALID mapping;
  - **agent (LLM, no FG rule)**    — reaches a valid repair (W2) from chemical
                                     reasoning alone, without being told not to
                                     split functional groups.

The result: a blind optimizer needs explicit chemical rules to stay valid; the LLM
agent gets that for free (and is also more evaluation-efficient, though stochastic —
hence best-of-N). This *motivates* the functional-group-integrity rule in
``agent.rules`` and is the paper's "why chemical reasoning matters" figure.

CLI::

  python -m agent.ablation \\
    --mapping reference/PSBMA_20mer/PSBMA20_mapping.json \\
    --itp reference/PSBMA_20mer/PSBMA20.itp \\
    --aa-top reference/PSBMA_20mer/PSBMA_20mer_no_water.gro \\
    --aa-traj reference/PSBMA_20mer/PSBMA_20mer_no_water_skip10.xtc \\
    --cg-struct reference/PSBMA_20mer/PSBMA_20mer_no_water.gro --frame-stride 25
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import asdict, dataclass, field
from pathlib import Path

from agent import loop as _loop
from agent import repair as _repair
from agent import validate as _validate


# ---------- shared single-move action space ----------


def candidate_moves(state: _repair.MappingState) -> list[dict]:
    """Every single-role reassignment op: send a heavy atom (with its Hs) to a bonded
    neighbour role. This is the local-search action space greedy/random draw from."""
    _n, bpm = _loop._beads_per_monomer(state)
    adj = _validate._bonded_role_adjacency(state, bpm)
    moves: list[dict] = []
    for fr in range(1, bpm + 1):
        bead = state.bead(fr)
        heavy = [nm for i, nm in zip(bead["atom_indices"], bead["atom_names"])
                 if i in state.heavy_atoms]
        for nm in heavy:
            for to in sorted(adj.get(fr, ())):
                moves.append({"op": "reassign", "atom_names": [nm], "from_role": fr, "to_role": to})
    return moves


def _apply_valid(state, op, ref, rule_check):
    """Apply one op; return the new state if it is valid + rule-clean, else None."""
    try:
        cand = _loop.apply_ops(state, [op], ref_positions=ref)
    except Exception:
        return None
    if rule_check(cand):
        return None
    return cand


# ---------- deterministic baselines ----------


@dataclass
class SearchResult:
    method: str
    best_objective: float
    n_evaluations: int
    found_target: bool
    fg_valid: bool | None = None        # does the best mapping keep functional groups intact?


def greedy_search(state, evaluate_fn, *, ref_positions=None, rule_check=None,
                  max_steps=8, target: float | None = None):
    """Steepest-ascent single-move hill-climb under ``rule_check``: each step
    evaluate all rule-clean single-move candidates and take the best strictly-
    improving one; stop when none improves. Returns (SearchResult, best_state)."""
    rc = rule_check or _loop.default_rule_check
    best = state
    best_obj = float(evaluate_fn(state).objective)
    n_evals = 1
    for _ in range(max_steps):
        step_best, step_obj = None, best_obj
        for op in candidate_moves(best):
            cand = _apply_valid(best, op, ref_positions, rc)
            if cand is None:
                continue
            obj = float(evaluate_fn(cand).objective)
            n_evals += 1
            if obj < step_obj - 1e-12:
                step_best, step_obj = cand, obj
        if step_best is None:
            break                       # no single move improves → stuck
        best, best_obj = step_best, step_obj
    return (SearchResult("greedy", best_obj, n_evals,
                         target is not None and best_obj <= target), best)


def random_search(state, evaluate_fn, *, ref_positions=None, rule_check=None,
                  budget=40, seed=0, target: float | None = None):
    """Random walk over single moves (keep-best) for a fixed evaluation budget.
    Returns (SearchResult, best_state)."""
    rc = rule_check or _loop.default_rule_check
    rng = random.Random(seed)
    cur = best = state
    best_obj = float(evaluate_fn(state).objective)
    n_evals = 1
    while n_evals < budget:
        cand = _apply_valid(cur, rng.choice(candidate_moves(cur)), ref_positions, rc)
        if cand is None:
            continue
        obj = float(evaluate_fn(cand).objective)
        n_evals += 1
        if obj < best_obj:
            best, best_obj = cand, obj
        cur = cand                       # random walk
    return (SearchResult(f"random (seed {seed})", best_obj, n_evals,
                         target is not None and best_obj <= target), best)


# ---------- the ablation ----------


@dataclass
class AblationResult:
    baseline: float
    target: float
    rows: list[SearchResult] = field(default_factory=list)


def run_ablation(
    state, make_llm_policy, evaluate_fn,
    *,
    ref_positions=None,
    functional_groups=None,
    target: float = 0.0375,
    greedy_max_steps: int = 6,
    llm_restarts: int = 3,
    llm_max_iters: int = 10,
    llm_plateau_k: int = 5,
    on_result=None,
) -> AblationResult:
    """Three strategies over the same action space + objective, tracking chemical
    validity (functional-group integrity):

      1. greedy, **objective-only** (no FG constraint) — free to game the objective;
      2. greedy, **chem-constrained** (FG-integrity enforced) — must stay valid;
      3. the **LLM agent** (given NO FG rule) — does it stay valid on its own?
    """
    from agent import rules as _rules

    baseline = float(evaluate_fn(state).objective)
    rows: list[SearchResult] = []

    def fg_ok(st) -> bool:
        return not _rules.fg_violations(st, functional_groups)

    def _emit(r):
        rows.append(r)
        if on_result:
            on_result(r)

    def fg_rule_check(st):
        return _loop.default_rule_check(st) + [v.message for v in _rules.fg_violations(st, functional_groups)]

    # 1. greedy, objective-only — may split a functional group to lower the objective
    r, s = greedy_search(state, evaluate_fn, ref_positions=ref_positions,
                         rule_check=_loop.default_rule_check, max_steps=greedy_max_steps, target=target)
    _emit(SearchResult("greedy (objective-only)", r.best_objective, r.n_evaluations,
                       r.found_target, fg_ok(s)))

    # 2. greedy, chem-constrained — FG-integrity enforced
    r2, s2 = greedy_search(state, evaluate_fn, ref_positions=ref_positions,
                          rule_check=fg_rule_check, max_steps=greedy_max_steps, target=target)
    _emit(SearchResult("greedy (chem-constrained)", r2.best_objective, r2.n_evaluations,
                       r2.found_target, fg_ok(s2)))

    # 3. the LLM agent — given no FG rule; does chemical reasoning keep it valid?
    best, result = _loop.run_loop_restarts(
        state, make_llm_policy, evaluate_fn, restarts=llm_restarts,
        ref_positions=ref_positions, max_iters=llm_max_iters, plateau_k=llm_plateau_k)
    _emit(SearchResult("agent (LLM, no FG rule)", float(result.best_objective),
                       max(1, result.n_iterations), result.best_objective <= target, fg_ok(best)))

    return AblationResult(baseline=baseline, target=target, rows=rows)


def plot_ablation(result: AblationResult, out_path: str | Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(1.9 * len(result.rows) + 3, 4.8))
    labels = [r.method for r in result.rows]
    objs = [r.best_objective for r in result.rows]
    # colour by chemical validity: red = invalid (gamed the objective), green = valid
    colors = ["tab:red" if r.fg_valid is False else "tab:green" for r in result.rows]
    ax.axhline(result.baseline, color="0.4", ls="--", lw=1.2, label=f"start = {result.baseline:.4f}")
    ax.bar(range(len(objs)), objs, color=colors, alpha=0.8)
    for i, r in enumerate(result.rows):
        valid = "✗ INVALID\n(splits a group)" if r.fg_valid is False else "✓ chemically valid"
        ax.annotate(f"{r.best_objective:.4f}\n{r.n_evaluations} evals\n{valid}",
                    (i, r.best_objective), textcoords="offset points", xytext=(0, 6),
                    ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=15, ha="right", fontsize=8.5)
    ax.set_ylim(0, max(objs) * 1.42)        # headroom so annotations clear the title
    ax.set_ylabel("best Gaussianity error reached  (mean 1 − R²)")
    ax.set_title("Ablation — pure optimization games the objective; the agent's chemistry keeps it valid",
                 fontsize=11)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    plt.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


# ---------- CLI ----------


def main(argv: list[str] | None = None) -> None:
    from agent.evaluate import evaluate_state
    from agent.loop import LLMPolicy

    p = argparse.ArgumentParser(description="Agentic-vs-deterministic ablation")
    p.add_argument("--mapping", required=True)
    p.add_argument("--itp", required=True)
    p.add_argument("--aa-top", required=True)
    p.add_argument("--aa-traj", required=True, nargs="+")
    p.add_argument("--cg-struct", default=None)
    p.add_argument("--frame-stride", type=int, default=25)
    p.add_argument("--model", default="anthropic/claude-opus-4.8")
    p.add_argument("--restarts", type=int, default=3)
    p.add_argument("--target", type=float, default=0.0375)
    p.add_argument("--work-root", default=None)
    p.add_argument("--out-dir", default=None)
    args = p.parse_args(argv)

    state = _repair.load_state(args.mapping, args.itp)
    ref = _repair.load_ref_positions(args.cg_struct) if args.cg_struct else None
    mol = state.mapping.get("molecule", "MOL")
    work_root = Path(args.work_root or f"derived/{mol}/ablation")

    def evaluate_fn(s):
        return evaluate_state(s, aa_top=args.aa_top, aa_traj=args.aa_traj,
                              work_root=work_root, molecule=mol, frame_stride=args.frame_stride)

    from agent import rules

    def make_policy():
        return LLMPolicy(model=args.model, temperature=0.6)

    fgs = rules.FUNCTIONAL_GROUPS_PSBMA if mol.upper().startswith("PSBMA") else None
    print(f"ablation for {mol}  (functional-group integrity: {'on' if fgs else 'off'})\n")
    result = run_ablation(state, make_policy, evaluate_fn, ref_positions=ref,
                          functional_groups=fgs, target=args.target, llm_restarts=args.restarts,
                          on_result=lambda r: print(
                              f"  {r.method:<26} best {r.best_objective:.4f}  "
                              f"{r.n_evaluations:>4} evals  "
                              f"{'✓ chemically valid' if r.fg_valid else '✗ INVALID (splits a group)'}"))

    out_dir = Path(args.out_dir or f"derived/{mol}")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "ablation_report.json").write_text(json.dumps(asdict(result), indent=2))
    plot = plot_ablation(result, out_dir / "ablation.pdf")
    print(f"\n  report: {out_dir / 'ablation_report.json'}")
    print(f"  plot  : {plot}")


if __name__ == "__main__":
    main()
