"""SWE-agent-style repair loop controller (#4).

The improvement loop is deliberately minimal: the intelligence lives in the
*action space* and the *feedback observation*, not the loop (the SWE-agent
Agent-Computer-Interface lesson). Both hard parts already exist —

  - **action space (ACI)** = the mapping-edit verbs in ``agent.repair``
    (``reassign_atoms`` / ``merge_beads`` / ``split_bead``), here exposed as
    monomer-tiled *ops* so one call edits every identical monomer at once;
  - **feedback observation** = the per-term Gaussianity report from
    ``agent.evaluate`` (per-term R² + the scalar objective).

so ``run_loop`` is a short ``while``:

    obs = evaluate(mapping)
    while budget left and not converged:
        action = policy.propose(obs)          # one edit (>=1 tiled ops), or submit
        cand   = apply_ops(mapping, action)    # atomic; reject on invalid / rule break
        obs    = evaluate(cand); keep_best_valid(cand)
        log(trajectory)

The **policy is swappable**: ``LLMPolicy`` is the agentic driver (chemistry
reasoning), ``ScriptedPolicy`` replays a fixed edit list (tests + a trivial
deterministic baseline). Same controller, swap the policy → the
agentic-vs-deterministic ablation falls out for free.

Design choices (see PROGRESS.md §4):
  - **hill-climb from best**: an edit that lowers the objective is adopted; a
    rule-valid edit that does not is evaluated (shown to the policy as feedback)
    then discarded, so the next proposal always starts from the best mapping — no
    steps wasted on manual reverts. Coupled multi-move plans stay reachable because
    a single action can carry several ops applied atomically. The best valid
    mapping is returned.
  - **Martini rules as hard constraints**: an edit that breaks a rule (or the
    atom partition) is rejected and the *reason* is fed back as the next
    observation — SWE-agent's lint-error trick — so the agent self-corrects.
  - **termination**: budget (max_iters) ∨ plateau (no best improvement over K) ∨
    the policy calling ``submit``. The full trajectory is the provenance +
    reproducibility artifact.
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Protocol

from agent import repair as _repair


# ---------- action space (monomer-tiled ops over the repair verbs) ----------


@dataclass
class Action:
    """One policy step: apply ``ops`` atomically, or ``submit`` to stop."""
    thought: str = ""
    ops: list[dict] = field(default_factory=list)
    submit: bool = False


def _check_role(role: int, bpm: int) -> None:
    if not 1 <= role <= bpm:
        raise ValueError(f"role {role} out of range 1..{bpm} for the current monomer")


def _beads_per_monomer(state: _repair.MappingState) -> tuple[int, int]:
    n_mon = int(state.mapping.get("n_monomers", 1) or 1)
    if state.n_beads % n_mon != 0:
        raise ValueError(
            f"bead count {state.n_beads} not divisible by n_monomers {n_mon}; "
            f"cannot tile a role-based op"
        )
    return n_mon, state.n_beads // n_mon


def _d2(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def _atoms_for_names(state, bead, names, ref_positions) -> list[int]:
    """Resolve heavy-atom ``names`` in ``bead`` to indices, **auto-carrying each
    heavy atom's hydrogens** — a hydrogen follows if its nearest in-bead heavy atom
    (by the reference frame) is one of the named atoms. So the policy names only
    heavy atoms and moves whole CHₙ groups. Falls back to exactly the named atoms
    when no reference positions are supplied; explicitly-named atoms always move.
    """
    idx_name = dict(zip(bead["atom_indices"], bead["atom_names"]))
    named = [i for i in bead["atom_indices"] if idx_name[i] in names]
    result = set(named)
    heavy_named = [i for i in named if i in state.heavy_atoms]
    if ref_positions and heavy_named:
        heavy_in_bead = [i for i in bead["atom_indices"] if i in state.heavy_atoms]
        for i in bead["atom_indices"]:
            if i in state.heavy_atoms or i in result:
                continue
            pi = ref_positions.get(i)
            if pi is None:
                continue
            nearest = min(heavy_in_bead, key=lambda h: _d2(pi, ref_positions.get(h, pi)))
            if nearest in heavy_named:
                result.add(i)
    return sorted(result)


def apply_ops(
    state: _repair.MappingState,
    ops: list[dict],
    ref_positions: dict[int, tuple[float, float, float]] | None = None,
) -> _repair.MappingState:
    """Apply a list of role-based ops, each **tiled across every monomer**.

    Ops (roles are 1-based *within a monomer*; the same edit is applied to all
    ``n_monomers`` identical repeat units):
      - ``{"op": "reassign", "atom_names": [...], "from_role": r, "to_role": r}``
      - ``{"op": "merge", "roles": [r, ...]}``
      - ``{"op": "split", "role": r, "group_a_atom_names": [...]}``

    Returns a new state (input untouched). Raises ``ValueError`` on any invalid op
    — the controller turns that into a rejected step + feedback. Merges/splits are
    processed high-monomer-first so lower monomers' bead indices stay valid across
    the reindexing.
    """
    cur = state
    for op in ops:
        n_mon, bpm = _beads_per_monomer(cur)
        kind = op.get("op")
        if kind == "reassign":
            names = set(op["atom_names"])
            fr, to = int(op["from_role"]), int(op["to_role"])
            _check_role(fr, bpm)
            _check_role(to, bpm)
            for m in range(n_mon):
                bead = cur.bead(m * bpm + fr)
                atoms = _atoms_for_names(cur, bead, names, ref_positions)
                if atoms:
                    cur = _repair.reassign_atoms(cur, atoms, m * bpm + to)
        elif kind == "merge":
            roles = sorted(int(r) for r in op["roles"])
            if len(roles) < 2:
                raise ValueError("merge needs >= 2 roles")
            for r in roles:
                _check_role(r, bpm)
            for m in range(n_mon - 1, -1, -1):
                cur = _repair.merge_beads(cur, [m * bpm + r for r in roles])
        elif kind == "split":
            role = int(op["role"])
            names = set(op["group_a_atom_names"])
            _check_role(role, bpm)
            for m in range(n_mon - 1, -1, -1):
                bead = cur.bead(m * bpm + role)
                ga = _atoms_for_names(cur, bead, names, ref_positions)
                cur = _repair.split_bead(cur, m * bpm + role, ga, ref_positions=ref_positions)
        else:
            raise ValueError(f"unknown op {kind!r} (expected reassign/merge/split)")
    return cur


# ---------- rule check (hard constraints) ----------


def default_rule_check(state: _repair.MappingState) -> list[str]:
    """v1 Martini hard constraints: no bead over 4 heavy atoms (regular-bead max).

    (Empty beads and a broken atom partition are already rejected inside the
    ``repair`` verbs.) Functional-group integrity + R/S/T exact sizing land with
    the dedicated rule checker (#3); this is the constraint that guards the loop
    from the size-violating Gaussianity optimum, e.g. the 5-heavy PSBMA case.
    """
    bad = []
    for b in state.beads:
        n = b["heavy_atom_count"]
        if n > 4:
            bad.append(f"bead {b['bead_id']} ({b.get('bead_name')}) has {n} heavy atoms (>4)")
    return bad


# ---------- observation ----------


def format_observation(
    state: _repair.MappingState,
    ev,
    *,
    best_obj: float,
    it: int,
    max_iters: int,
    last_result: str | None,
) -> str:
    """Render the feedback the policy sees each step: objective, worst terms,
    and one monomer's bead layout (roles + atom names) so it can pick an edit."""
    lines = [
        f"iteration {it}/{max_iters} | current objective {ev.objective:.4f} "
        f"(lower is better) | best so far {best_obj:.4f}"
    ]
    if last_result:
        lines.append(f"result of last action: {last_result}")
    terms = sorted(getattr(ev.breakdown, "terms", []), key=lambda t: t.error, reverse=True)
    if terms:
        lines.append("per-term Gaussianity, worst first (err = 1 - R², 0 = perfect Gaussian):")
        for t in terms[:8]:
            lines.append(f"    {t.kind:<5} {t.label:<18} R²={t.r2:.3f}  err={t.error:.3f}")
    try:
        n_mon, bpm = _beads_per_monomer(state)
        lines.append(f"bead layout (one of {n_mon} identical monomers; roles 1..{bpm}; "
                     f"heavy atoms only — their H's move with them):")
        for r in range(1, bpm + 1):
            b = state.bead(r)
            heavy = [nm for i, nm in zip(b["atom_indices"], b["atom_names"])
                     if i in state.heavy_atoms]
            lines.append(f"    role {r}: {b['bead_type']:<5} heavy_atoms={heavy}")
        intra = sorted({(min(bd.i, bd.j), max(bd.i, bd.j)) for bd in state.bonds
                        if bd.i <= bpm and bd.j <= bpm})
        if intra:
            lines.append("  bonded role pairs (reassign/merge only work between bonded "
                         "roles): " + ", ".join(f"{i}-{j}" for i, j in intra))
    except ValueError:
        lines.append(f"(irregular monomer structure: {state.n_beads} beads)")
    return "\n".join(lines)


# ---------- policies ----------


class _EvalResult(Protocol):
    """Structural type for whatever ``evaluate_fn`` returns (real EvalResult or a
    test fake): a scalar ``objective`` plus a ``breakdown`` with per-term stats."""
    objective: float
    breakdown: object


class Policy(Protocol):
    def propose(self, observation: str) -> Action: ...


class ScriptedPolicy:
    """Replay a fixed list of ``Action``s (then ``submit``). Deterministic — used
    for tests and as a trivial baseline in the agentic-vs-deterministic ablation."""

    def __init__(self, actions: list[Action]):
        self._actions = list(actions)
        self._i = 0

    def propose(self, observation: str) -> Action:
        if self._i >= len(self._actions):
            return Action(thought="scripted: done", submit=True)
        a = self._actions[self._i]
        self._i += 1
        return a


_SYSTEM_PROMPT = """\
You refine a Martini 3 coarse-grained MAPPING of a polymer by regrouping which \
all-atom atoms belong to which CG bead. Objective: make every bead bond-length \
and bead angle distribution — measured by projecting the AA trajectory through \
the mapping — as close to a single Gaussian as possible, i.e. MAXIMIZE each \
term's R² (minimize err = 1 - R²), lowering the mean objective.

Physical intuition:
- A broad or bimodal bond/angle means a rotatable/floppy degree of freedom is \
buried inside a bead or split awkwardly across a bead boundary. Fix it by moving \
the boundary atom(s) to a neighbouring bead, or by merging/splitting so the \
rotatable bond becomes a clean inter-bead bond.
- Re-centring a charged bead on its charged atom often needs a COUPLED SWAP in one \
atomic edit: push one substituent OUT to a neighbour and pull an atom from the OTHER \
neighbour IN, so the bead stays <= 4 heavy but its centroid shifts onto the charge. \
Consider every bonded boundary, including the aliphatic linker/tail beads, not just \
the polar groups.

Hard constraints (Martini): no bead may exceed 4 heavy atoms; do not split a \
chemical functional group (ester, sulfonate, ammonium) across beads. Rule-breaking \
edits are rejected and returned to you — fix and retry.

Name ONLY heavy atoms (C, N, O, S); each heavy atom's hydrogens move with it \
automatically. The bead layout lists heavy atoms only.

Every op is TILED automatically across all identical monomers: give roles within \
ONE monomer (1-based) and heavy-atom names. You may pass several ops in one \
edit_mapping call; they apply atomically, so a two-move swap that keeps every bead \
<= 4 heavy is allowed even if either move alone would break the rule — this is often \
how to re-centre a charged bead (move one heavy atom in and another out together).

You ALWAYS work from the best mapping found so far. An edit that lowers the \
objective is kept; an edit that does NOT is automatically discarded (you still see \
its score as feedback, then you are back at the best mapping). So never revert by \
hand, and never repeat an edit that just failed — each turn propose a genuinely \
different idea (a different bead boundary, a coupled two-move swap, a merge or \
split), or submit when you have run out of chemically-motivated moves.

Each turn: reason briefly, then call edit_mapping with one or more ops, or submit \
when no further chemically-motivated improvement is available. Prefer small, \
justified edits."""

_EDIT_TOOL = {
    "name": "edit_mapping",
    "description": "Apply one or more atom-regrouping ops (tiled across all monomers), "
                   "atomically; the mapping is then re-scored and returned to you.",
    "input_schema": {
        "type": "object",
        "properties": {
            "thought": {"type": "string", "description": "brief chemical rationale"},
            "ops": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "op": {"type": "string", "enum": ["reassign", "merge", "split"]},
                        "atom_names": {"type": "array", "items": {"type": "string"}},
                        "from_role": {"type": "integer"},
                        "to_role": {"type": "integer"},
                        "roles": {"type": "array", "items": {"type": "integer"}},
                        "role": {"type": "integer"},
                        "group_a_atom_names": {"type": "array", "items": {"type": "string"}},
                    },
                    "required": ["op"],
                },
            },
        },
        "required": ["thought", "ops"],
    },
}

_SUBMIT_TOOL = {
    "name": "submit",
    "description": "Stop refining; the best valid mapping seen is returned.",
    "input_schema": {
        "type": "object",
        "properties": {"thought": {"type": "string"}},
    },
}


def _openai_tools() -> list[dict]:
    """The two tools in OpenAI function-calling format (from the schemas above)."""
    return [
        {"type": "function", "function": {"name": t["name"], "description": t["description"],
                                          "parameters": t["input_schema"]}}
        for t in (_EDIT_TOOL, _SUBMIT_TOOL)
    ]


class LLMPolicy:
    """Agentic driver: an LLM proposes edits via tool-use over the observation.

    Speaks the **OpenAI-compatible** chat/tool protocol (via ``requests``), so it
    works with OpenRouter (the default), a native OpenAI endpoint, or any compatible
    gateway — no vendor SDK needed. Defaults to OpenRouter + ``OPENROUTER_API_KEY``
    and ``anthropic/claude-sonnet-4.6``. The conversation is stateful: each fresh
    observation is delivered as the ``tool`` result for the model's previous tool
    call, so it sees the consequence (accepted Δ / rejection reason) of its last edit.
    """

    DEFAULT_BASE_URL = "https://openrouter.ai/api/v1"

    def __init__(
        self,
        model: str = "anthropic/claude-sonnet-4.6",
        *,
        base_url: str | None = None,
        api_key_env: str = "OPENROUTER_API_KEY",
        max_tokens: int = 1500,
        temperature: float = 0.0,
        timeout: int = 180,
    ):
        try:
            import requests
        except ImportError as e:
            raise ImportError("LLMPolicy needs `requests` (pip install requests)") from e
        self._requests = requests
        self._model = model
        self._base_url = (base_url or os.environ.get("OPENAI_BASE_URL")
                          or self.DEFAULT_BASE_URL).rstrip("/")
        self._api_key = os.environ.get(api_key_env) or os.environ.get("OPENAI_API_KEY")
        if not self._api_key:
            raise RuntimeError(f"no API key found: set ${api_key_env}")
        self._max_tokens = max_tokens
        self._temperature = temperature
        self._timeout = timeout
        self._tools = _openai_tools()
        self._messages: list[dict] = [{"role": "system", "content": _SYSTEM_PROMPT}]
        self._pending_tool_call_id: str | None = None

    def propose(self, observation: str) -> Action:
        if self._pending_tool_call_id is None:
            self._messages.append({"role": "user", "content": observation})
        else:
            self._messages.append({"role": "tool", "tool_call_id": self._pending_tool_call_id,
                                   "content": observation})
        resp = self._requests.post(
            self._base_url + "/chat/completions",
            headers={"Authorization": f"Bearer {self._api_key}",
                     "Content-Type": "application/json",
                     "X-Title": "autoMartiniAgent repair loop"},
            json={"model": self._model, "max_tokens": self._max_tokens,
                  "temperature": self._temperature, "messages": self._messages,
                  "tools": self._tools, "tool_choice": "auto"},
            timeout=self._timeout,
        )
        resp.raise_for_status()
        data = resp.json()
        if "choices" not in data:
            raise RuntimeError(f"unexpected API response: {json.dumps(data)[:400]}")
        msg = data["choices"][0]["message"]
        tool_calls = msg.get("tool_calls") or []
        if not tool_calls:  # model chose not to call a tool → treat as submit
            self._messages.append({"role": "assistant", "content": msg.get("content") or ""})
            self._pending_tool_call_id = None
            return Action(thought=(msg.get("content") or "(no tool call)"), submit=True)
        tc = tool_calls[0]
        # keep exactly one tool_call in history so we owe exactly one tool response
        self._messages.append({
            "role": "assistant", "content": msg.get("content"),
            "tool_calls": [{"id": tc["id"], "type": "function",
                            "function": {"name": tc["function"]["name"],
                                         "arguments": tc["function"]["arguments"]}}],
        })
        self._pending_tool_call_id = tc["id"]
        try:
            args = json.loads(tc["function"]["arguments"] or "{}")
        except json.JSONDecodeError:
            args = {}
        if tc["function"]["name"] == "submit":
            return Action(thought=args.get("thought", ""), submit=True)
        return Action(thought=args.get("thought", ""), ops=args.get("ops", []) or [])


# ---------- controller ----------


@dataclass
class TrajectoryStep:
    iteration: int
    thought: str
    ops: list[dict]
    accepted: bool          # rule-valid and applied (working state advanced)
    improved: bool          # produced a new best
    objective: float | None  # objective after the edit (None if rejected)
    best_objective: float
    note: str               # human-readable result / rejection reason


@dataclass
class LoopResult:
    initial_objective: float
    best_objective: float
    n_iterations: int
    stop_reason: str        # "submit" | "plateau" | "budget"
    trajectory: list[TrajectoryStep] = field(default_factory=list)


def run_loop(
    state: _repair.MappingState,
    policy: Policy,
    evaluate_fn: Callable[[_repair.MappingState], _EvalResult],
    *,
    ref_positions: dict[int, tuple[float, float, float]] | None = None,
    rule_check: Callable[[_repair.MappingState], list[str]] = default_rule_check,
    max_iters: int = 12,
    plateau_k: int = 4,
) -> tuple[_repair.MappingState, LoopResult]:
    """Drive ``policy`` over ``state``, returning the best valid mapping + a log.

    ``evaluate_fn(state)`` must return an object with ``.objective`` (float, lower =
    better) and ``.breakdown.terms`` (per-term R² for the observation). For real
    runs wrap ``agent.evaluate.evaluate_state``; tests inject a fast fake.
    """
    working = state
    working_ev = evaluate_fn(working)
    best_state, best_obj, best_ev = working, working_ev.objective, working_ev
    initial_obj = working_ev.objective
    traj: list[TrajectoryStep] = []
    last_result: str | None = None
    iters_since_improve = 0
    stop_reason = "budget"

    for it in range(1, max_iters + 1):
        obs = format_observation(
            working, working_ev, best_obj=best_obj, it=it, max_iters=max_iters,
            last_result=last_result,
        )
        action = policy.propose(obs)

        if action.submit:
            stop_reason = "submit"
            traj.append(TrajectoryStep(it, action.thought, [], False, False, None,
                                       best_obj, "submit"))
            break

        # apply (atomic) → rule-check → evaluate
        try:
            cand = apply_ops(working, action.ops, ref_positions=ref_positions)
        except Exception as e:  # invalid op / partition break
            last_result = f"REJECTED (invalid edit): {e}"
            traj.append(TrajectoryStep(it, action.thought, action.ops, False, False,
                                       None, best_obj, last_result))
            continue

        violations = rule_check(cand)
        if violations:
            last_result = "REJECTED (Martini rule): " + "; ".join(violations[:3])
            traj.append(TrajectoryStep(it, action.thought, action.ops, False, False,
                                       None, best_obj, last_result))
            continue

        cand_ev = evaluate_fn(cand)
        improved = cand_ev.objective < best_obj - 1e-12
        if improved:
            best_state, best_obj, best_ev = cand, cand_ev.objective, cand_ev
            working, working_ev = cand, cand_ev
            iters_since_improve = 0
            last_result = (f"ACCEPTED, new best objective {cand_ev.objective:.4f} "
                           f"(Δ {cand_ev.objective - initial_obj:+.4f} vs start)")
        else:
            # hill-climb: discard the non-improving edit; keep proposing from best
            working, working_ev = best_state, best_ev
            iters_since_improve += 1
            last_result = (f"DISCARDED: objective {cand_ev.objective:.4f} did not beat best "
                           f"{best_obj:.4f}; reverted to best — try a different edit")
        traj.append(TrajectoryStep(it, action.thought, action.ops, improved, improved,
                                   cand_ev.objective, best_obj, last_result))

        if iters_since_improve >= plateau_k:
            stop_reason = "plateau"
            break

    result = LoopResult(
        initial_objective=initial_obj, best_objective=best_obj,
        n_iterations=len(traj), stop_reason=stop_reason, trajectory=traj,
    )
    return best_state, result


def write_trajectory(result: LoopResult, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(asdict(result), indent=2))
    return path


# ---------- CLI ----------


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run the agentic CG-mapping repair loop")
    p.add_argument("--mapping", required=True)
    p.add_argument("--itp", required=True)
    p.add_argument("--aa-top", required=True)
    p.add_argument("--aa-traj", required=True, nargs="+")
    p.add_argument("--cg-struct", default=None, help="AA .gro for split ref positions (optional)")
    p.add_argument("--policy", choices=["llm", "scripted"], default="llm")
    p.add_argument("--scripted-ops", default=None,
                   help="JSON file: list of actions [{thought, ops:[...]}] for --policy scripted")
    p.add_argument("--model", default="anthropic/claude-sonnet-4.6",
                   help="OpenRouter/OpenAI-compatible model slug (default anthropic/claude-sonnet-4.6)")
    p.add_argument("--base-url", default=None,
                   help="OpenAI-compatible base URL (default OpenRouter, or $OPENAI_BASE_URL)")
    p.add_argument("--api-key-env", default="OPENROUTER_API_KEY",
                   help="env var holding the API key (default OPENROUTER_API_KEY)")
    p.add_argument("--frame-stride", type=int, default=25)
    p.add_argument("--max-iters", type=int, default=12)
    p.add_argument("--plateau-k", type=int, default=4)
    p.add_argument("--work-root", default=None)
    p.add_argument("--out", default=None, help="trajectory JSON output path")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    from agent.evaluate import evaluate_state

    args = _parse_args(argv)
    state = _repair.load_state(args.mapping, args.itp)
    mol = state.mapping.get("molecule", "MOL")
    work_root = Path(args.work_root or f"derived/{mol}/loop")
    ref = _repair.load_ref_positions(args.cg_struct) if args.cg_struct else None

    def evaluate_fn(st):
        return evaluate_state(
            st, aa_top=args.aa_top, aa_traj=args.aa_traj, work_root=work_root,
            molecule=mol, frame_stride=args.frame_stride,
        )

    if args.policy == "scripted":
        raw = json.loads(Path(args.scripted_ops).read_text())
        policy: Policy = ScriptedPolicy([Action(**a) for a in raw])
    else:
        policy = LLMPolicy(model=args.model, base_url=args.base_url,
                           api_key_env=args.api_key_env)

    best_state, result = run_loop(
        state, policy, evaluate_fn, ref_positions=ref,
        max_iters=args.max_iters, plateau_k=args.plateau_k,
    )

    print(f"loop finished ({result.stop_reason}) after {result.n_iterations} step(s)")
    print(f"  objective {result.initial_objective:.4f} -> {result.best_objective:.4f} "
          f"(Δ {result.best_objective - result.initial_objective:+.4f})")
    out_dir = work_root
    _repair.write_mapping(best_state, out_dir / f"{mol}_loop_best_mapping.json")
    _repair.write_itp(best_state, out_dir / f"{mol}_loop_best.itp", mol_name=mol)
    traj_path = Path(args.out) if args.out else out_dir / "trajectory.json"
    write_trajectory(result, traj_path)
    print(f"  best mapping: {out_dir / f'{mol}_loop_best_mapping.json'}")
    print(f"  trajectory  : {traj_path}")


if __name__ == "__main__":
    main()
