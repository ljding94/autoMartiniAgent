# program.md — AA → Martini 3 Mapping Protocol

> **Version**: 0.3 (draft, 2026-05-01) · uncalibrated thresholds, see §Calibration TODOs
> **Audience**: any MCP-aware agent (Claude Code, Codex CLI, OpenCode, Cursor, Continue, …) with the `auto-martini-agent` MCP server mounted and the skill (`SKILL.md` / `AGENTS.md`) loaded.
> **Pattern**: agent-runnable autoresearch *program* in Karpathy's `program.md` lineage — one file, fully self-contained, end-to-end. Agent-agnostic: tools resolve through MCP, no Claude-Code-specific syntax.

## Goal

Given an atomistic molecular structure and an AA simulation trajectory of that molecule, produce a Martini 3 coarse-grained **mapping** — atom-index → bead assignment, with bead type and bead size — that:

1. Satisfies Martini 3 hard rules (R/S/T sizing, functional-group integrity, symmetry, bead-count plausibility).
2. Yields bond and angle distributions, when computed on the AA trajectory projected through the mapping, that are well-described by single Gaussians (i.e., consistent with harmonic Martini bonded potentials).

The mapping is the deliverable. No CG simulation is run. No bonded force-field parameters are fit (that is a downstream post-process).

## Inputs (provided by the user)

- `structure`: SMILES string, or `.mol` / `.pdb` file with explicit atom indexing.
- `aa_traj`: AA trajectory file (`.xtc`, `.trr`, `.dcd`).
- `aa_top`: matching topology (`.tpr`, `.top`, `.psf`, or equivalent).
- `mol_name` (optional): a 3–4 char name for the residue. Validated against `martini_v3.0.0_solvents_v1.itp` reserved names.
- `budget` (optional): `{max_iters: 10, max_wall_minutes: 30}` defaults.

## Runtime requirements

Before this protocol runs, the agent's runtime must have:
- the `auto-martini-agent` MCP server registered and reachable (stdio transport),
- the `auto-martini-agent` skill loaded (`SKILL.md` for Claude Code, `AGENTS.md` for Codex/OpenCode, plain markdown otherwise),
- read access to `aa_traj` and `aa_top` paths (passed by path, not contents — too large for MCP payloads).

If any of these is missing the protocol must abort with a clear message — do not attempt manual workarounds.

## Tools available (MCP-registered, from `auto-martini-agent` MCP server)

| tool                      | input                                  | output                                                  |
|---------------------------|----------------------------------------|---------------------------------------------------------|
| `classify`                | `structure`                            | `{small, medium, oligomer, polymer, peptide}`           |
| `propose_mapping`         | `structure, category`                  | `mapping` (atom→bead, bead_type, bead_size)             |
| `project_trajectory`      | `aa_traj, aa_top, mapping`             | `cg_traj` (per-frame bead positions, mass-weighted COM) |
| `score_mapping`           | `cg_traj, mapping`                     | `{distributions, gaussianness, rule_violations, score}` |
| `repair_mapping`          | `mapping, score_report, structure`     | `mapping' + change_log`                                 |
| `martini_rules.lookup`    | `group_smiles`                         | recommended `{type, size}` for that fragment            |

Each tool is deterministic given its inputs (mod RNG seeds in the cold-start). Cache `project_trajectory` and `score_mapping` results keyed on `hash(mapping)` — the agent will revisit candidates.

## Acceptance criteria (all must hold to terminate as "success")

- **Distribution shape**, every bond and every angle:
  - `KL(empirical || best_fit_gaussian) < 0.05`
  - `|skew| < 0.5`
  - GMM-BIC selects 1 component (no bimodal signature)
- **Sizing rules**:
  - Regular (R) bead → contains exactly 4 heavy atoms
  - Small (S) bead → contains exactly 3 heavy atoms
  - Tiny (T) bead → contains exactly 2 heavy atoms
- **Functional-group integrity**: amide, ester, carboxylate, sulfonate, ammonium, sulfonamide groups are NOT split across beads.
- **Symmetry**: atoms equivalent under the molecule's symmetry group map to symmetry-equivalent beads.
- **Bead count plausibility**: total bead count is within ±25% of `heavy_atoms / 4`.

(Thresholds are first-cut; recalibrate once Phase 0 fixtures land — see §Calibration TODOs.)

## Protocol

```
1. classify
2. cold-start mapping (with stall recovery)
3. project AA trajectory
4. score
5. decide → repair and loop, or terminate
```

### Step 1 — Classify

Call `classify(structure)`. Record the category. If the user asserted `--repeat-unit`, force `category = polymer` regardless of heuristic output.

### Step 2 — Cold-start mapping

Branch on category:

- `small` or `medium`: call `propose_mapping(structure, category)`.
  - This delegates to AutoMARTINI3 with `--fpred` (always — ALOGPS lookups fail on charged/exotic groups otherwise) and a 5-minute subprocess timeout.
  - On timeout: fall back to **fragment-and-assemble**:
    1. RDKit BRICS fragmentation of the structure.
    2. Run AutoMARTINI3 per fragment (these tend to be small and converge).
    3. Reconnect fragments using bonds inferred from the AA topology.
- `polymer`: propose monomer-level mapping first. Defer polymer assembly (Polyply) to a later iteration. Note in the mapping spec that the bonded params at the inter-monomer boundary will need separate treatment.
- `peptide`: call Martinize2 stub. Protocol terminates here unless explicitly extended — peptide mappings are largely standardized and not the research focus.

Record the cold-start mapping as iteration 0.

### Step 3 — Project AA trajectory

Call `project_trajectory(aa_traj, aa_top, mapping)`. This computes per-frame mass-weighted center-of-mass for each CG bead.

Sanity check: if `n_frames < 1000` or the trajectory is shorter than 100 ns, log `confidence: low, reason: insufficient_sampling` and continue, but flag in the final report.

### Step 4 — Score

Call `score_mapping(cg_traj, mapping)`. Read out:

- per-bond Gaussianness (KL, skew, kurtosis, GMM components)
- per-angle Gaussianness
- rule violations (sizing, functional-group split, symmetry, bead count)
- scalar score (lower is better; sum of normalized term penalties)

### Step 5 — Decide

- **All acceptance criteria pass** → terminate. Return `{mapping, score_report, provenance_log, confidence: high}`.
- **Otherwise** → call `repair_mapping(mapping, score_report, structure)`, increment iteration, return to step 3.
- **Budget exhausted** (≥10 iterations or ≥30 wall-minutes) → terminate with `confidence: low`. Return the best-scored mapping seen, plus the list of remaining violations for human review.

## Decision rules used by `repair_mapping`

Consulted in this order — first match wins:

1. **Sizing violation, under-filled bead** (e.g., R-typed bead with only 2–3 heavy atoms):
   - If an adjacent bead exists such that combined heavy-atom count ≤ 4, **merge**.
   - Else **relabel** to S (3 heavy) or T (2 heavy) as appropriate.
2. **Sizing violation, over-filled bead** (≥ 5 heavy atoms):
   - **Split** along a rotatable bond. Prefer splits that preserve functional-group integrity.
3. **Bimodal angle distribution** (GMM-BIC ≥ 2 components):
   - A rotatable bond is being projected through a bead. Re-fragment so the rotatable bond becomes an inter-bead bond (i.e., a CG bond) rather than an internal degree of freedom. (See [Agentic Martini idea note] for the underlying argument.)
4. **Skewed bond distribution** (`|skew| > 0.5`):
   - Asymmetric atom assignment. Try shifting one heavy atom from the heavier bead to the lighter neighbor.
5. **Functional-group split**:
   - Recombine into the bead containing the functional group's anchor atom (carbonyl C for amide/ester, S for sulfonate, P for phosphate, N for ammonium).

If multiple rules apply to a single bead, apply them in the listed order, one repair per iteration. Log every change in `provenance_log` with rule index + rationale.

## Output

On success or budget exhaustion, write to `out/`:

- `mapping.json` — `{atoms: [{idx, name, element, bead}], beads: [{id, type, size, charge, smiles_fragment}]}`
- `mapping.itp` — GROMACS-format Martini 3 `.itp` skeleton (`[atoms]` filled, `[bonds]`/`[angles]` left as TODO for downstream)
- `score_report.json` — full distributions + Gaussianness scores + rule check
- `provenance.log` — chronological list of `{iter, action, rule, before, after, rationale}` entries
- `summary.md` — human-readable summary (one figure: AA distributions overlaid with best-fit Gaussian per bond/angle)

## Failure recovery

| symptom                                         | response                                                                  |
|-------------------------------------------------|---------------------------------------------------------------------------|
| AutoMARTINI3 stall past 5 min                   | kill subprocess; fall back to fragment-and-assemble                       |
| ALOGPS lookup failure                           | already mitigated by always passing `--fpred`; if it still fails, mark fragment polarity unknown and call `martini_rules.lookup` |
| AA traj too short / under-converged             | continue with `confidence: low`; flag in summary                          |
| `mol_name` collides with `martini_v3.0.0_solvents_v1.itp` | rename automatically (suffix `_X`); log the rename                |
| Repair budget exhausted with violations remaining | return best-scored mapping; populate `remaining_violations` in summary  |
| `score_mapping` returns NaN / inf for any term  | inspect distribution: usually a single-frame degenerate bead (overlapping atoms); flag and skip that term |

## Calibration TODOs

These thresholds and budgets are first-cut. Calibrate against Phase 0 fixtures and a held-out set of human-validated Martini 3 mappings (Martini library SI, M3 small-molecules paper SI):

- [ ] KL threshold (currently 0.05) — measure distribution on validated mappings.
- [ ] Skew bound (currently 0.5) — same.
- [ ] AutoMARTINI3 subprocess timeout (currently 5 min) — measure runtime distribution on neutral small molecules.
- [ ] Iteration budget (currently 10) — measure typical convergence iterations on Phase 0 fixtures.
- [ ] Wall-clock budget (currently 30 min) — same.
- [ ] AA-trajectory sampling-sufficiency thresholds (currently 1000 frames / 100 ns).

## Provenance

This protocol implements the loop described in the [Agentic Martini idea note](../Note-Work/Ideas/Agentic%20Martini.md): mapping is the optimization variable, AA-projected bonded-distribution Gaussianness is the objective, with hard Martini-3 rules as constraints. Tools and scope match the autoMartiniAgent project plan in `PROGRESS.md` (§2 scope, §3 tools, §4 phases).
