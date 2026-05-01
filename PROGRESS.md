# autoMartiniAgent — Progress Tracker

**Project goal**: build an agent-runnable research protocol (`program.md`) that, combined with the `auto-martini-agent` Claude Code skill, automatically converts AA simulation data + structure into a Martini 3 coarse-grained mapping (atom-to-bead assignment with bead types and sizes).

The mapping is the deliverable. We do not run CG simulations. Validation is in-the-loop: project the AA trajectory through the proposed mapping and inspect bead–bead distance and angle distributions.

Collaborators (ORNL): Lijie Ding (driver), Seonghan Kim (pipeline + reproduction notes, AA data), Chris Walker (chemistry validation, AA data), Jan Michael Carrillo.

---

## 1. What the email chain established

**Seonghan's reproducible pipeline (verified on 1-octanol):**
- env: conda `autom3`, Python 3.10, GROMACS ≥ 2021
- backend: [Auto-MartiniM3](https://github.com/Martini-Force-Field-Initiative/Automartini_M3) (`pip install -e .` from source) + RDKit
- force field files staged in `ff/`: `martini_v3.0.0.itp`, `martini_v3.0.0_solvents_v1.itp`, `water.gro`
- Pipeline: `auto_martiniM3 --smi … --mol … --canon -v` → `gmx insert-molecules` → `gmx solvate -radius 0.21` → `gmx grompp` → `gmx mdrun` (steep EM, emtol 100)
- Naming gotcha: `MOL=OCT` collides with octane in `martini_v3.0.0_solvents_v1.itp` → use `OCOL`. Always check solvents.itp before naming.
- (The CG-sim half of this pipeline is now out of scope for our deliverable, but the `.itp` generation half stays.)

**Seonghan's proposed three-stage workflow** (we adopt the routing logic, drop the CG-sim validation):
- *Initial CG generation* — classifier routes by molecule type:
  - small (≤25 heavy atoms) → AutoMARTINI3
  - medium / oligomer → fragment + AutoMARTINI3 + assembly
  - polymer / repeat unit → Polyply
  - protein / peptide → Martinize2
- *Refinement* — atomistic ref traj → bonded distribution targets → BI or Swarm-CG. **In our scope this is a post-process on the final mapping, not part of the inner loop.**
- *Validation via CG sim* — **dropped**. Replaced by AA-projected distribution scoring.

**Chris's stress-test (the real research problem):**
- PSBMA monomer (zwitterionic, ammonium + sulfonate, ~22 heavy atoms, SMILES `CC(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]`): AutoMARTINI3 hangs >30 min, no convergence.
- Without `--fpred`, ALOGPS fragment lookup fails → flag is mandatory in our wrapper.
- Smaller PMETAC (`CC(C)C(=O)OCC[N+](C)(C)`) succeeds in ~1 min but the mapping violates Martini 3 sizing rules: emits Regular (R) bead labels for fragments holding only 2–3 heavy atoms. Example output:

  | bead | type | atoms       | heavy | should be |
  |------|------|-------------|-------|-----------|
  | C01  | C3   | CCC         | 3     | S-bead    |
  | P01  | P1   | OC=O        | 3     | S-bead    |
  | C02  | C5   | CC          | 2     | T-bead    |
  | Q01  | Qd   | C[N+]C      | 3     | S-bead    |

  Martini 3 sizing: R = 4 heavy atoms, S = 3, T = 2.

- Required capability: **(a) change the number of CG beads, (b) re-pick bead size class, (c) recover when the underlying tool stalls.**

---

## 2. Scope (revised 2026-04-30)

**In scope**
- Input: AA trajectory + topology + chemical structure (SMILES / mol / pdb).
- Output: a *mapping* — atom-index → bead, bead type, bead size — plus a score report and provenance log.
- Inner-loop validator: AA→CG projection followed by bond/angle distribution scoring (Gaussianness + Martini-rule compliance).
- Repair loop on top of AutoMARTINI3 / Polyply: relabel size class, re-fragment, recover from stalls.
- A `program.md` protocol file that lets the agent run the full procedure end-to-end given inputs.
- Demonstration on Chris's PMETAC and PSBMA AA trajectories.

**Out of scope (for v1)**
- Running CG simulations (`gmx mdrun`, EM, NVT, NPT).
- Fitting bonded force-field parameters as part of the loop. (Optional post-process only.)
- Free-energy / partition-coefficient validation.
- Membranes, proteins beyond a Martinize2 stub.
- A GUI; the deliverable is CLI + `program.md`.

**Headline demonstration**: drop a new monomer + AA trajectory + `program.md` into any MCP-aware agent (Claude Code, Codex CLI, OpenCode, Cursor, Continue, …) with the skill loaded → agent autonomously produces a validated Martini 3 mapping, no human in the loop. This is the autoresearch-style proof point, in the Karpathy `program.md` lineage.

---

## 3. Deliverable architecture

The product is **agent-agnostic** — three portable layers:

```
autoMartiniAgent/
├── mcp_server/         # Python MCP server (stdio), runs inside autom3 conda env
│   ├── pyproject.toml
│   └── src/auto_martini_mcp/
├── skill/
│   ├── SKILL.md        # Claude Code skill format
│   ├── AGENTS.md       # cross-agent community standard (mirror of SKILL.md)
│   └── scripts/        # deterministic helpers callable without an LLM
├── program.md         # the autoresearch protocol
├── tests/fixtures/
└── PROGRESS.md
```

**Layer 1 — MCP server.** Single Python process exposing the tools below over stdio MCP. Runs inside the `autom3` conda env so AutoMARTINI3 + RDKit are available. Any MCP-aware agent mounts it. Large data (AA trajectories) passed by **path**, not contents.

**Layer 2 — Portable skill.** Markdown instructions on *when* and *how* to invoke the tools, plus reasoning prose (decision rules, recovery strategies). `SKILL.md` for Claude Code, `AGENTS.md` for Codex/OpenCode/etc. — same content, two filenames.

**Layer 3 — Research protocol.** `program.md`, plain markdown, agent-agnostic. References MCP tool names. No agent-specific syntax anywhere.

### MCP-registered tools

- `classify(structure)` → `category ∈ {small, medium, oligomer, polymer, peptide}`
- `propose_mapping(structure, category)` → initial mapping spec (cold-start via AutoMARTINI3 / Polyply / Martinize2 / fragment+assemble)
- `project_trajectory(aa_traj, aa_top, mapping)` → CG bead trajectory (per-frame mass-weighted COM)
- `score_mapping(cg_traj, mapping)` → `{bond_distributions, angle_distributions, gaussianness_per_term, rule_violations, scalar_score}`
- `repair_mapping(mapping, score_report, structure)` → revised mapping + change_log
- `martini_rules.lookup(group_smiles)` → recommended bead type/size for a chemical fragment

Portability invariant: nothing in the skill prose or `program.md` may reference Claude-Code-specific syntax, agent built-ins, or non-MCP tool names.

---

## 4. Phases

### Phase 0 — Environment + fixtures (week 1)
- Stand up `autom3` env per Seonghan's notes (Python 3.10, AutoMARTINI3 from source, RDKit). GROMACS not strictly required for the inner loop; needed only for replaying Seonghan's octanol reproduction.
- Reproduce 1-octanol `.itp` generation; capture as golden fixture.
- Reproduce Chris's PMETAC mis-sized output and PSBMA stall as regression fixtures.
- Ingest first AA trajectory from collaborators.
- **Deliverable**: `tests/fixtures/{octanol, pmetac, psbma}/` with `.itp`, AA `.xtc/.tpr` where available.

### Phase 1 — Read backends (week 1–2)
- AutoMARTINI3 source: fragmentation routine, bead-type assignment, the optimization loop that stalls on PSBMA, what `--fpred` swaps in, where size class is decided.
- Polyply: monomer `.itp` → polymer topology + coords; document input contract.
- Martinize2: identify the protein entry point (stub for v1).
- Martini 3 paper: bead types, R/S/T size rules, fragment-selection heuristics, ±1 polarity-shift soft rule.
- **Deliverable**: `notes/backends.md` — one-page contract per backend.

### Phase 2 — Classifier + dispatcher (week 2)
- Input parser: SMILES, `.mol`, `.pdb`, `--repeat-unit` flag.
- Heuristics: heavy-atom count via RDKit, peptide-bond detection, charge detection, ring/aromatic detection.
- Routing policy table (mirrors Seonghan's stage-1 logic).
- **Deliverable**: `agent/classify.py`, `agent/dispatch.py`, unit tests.

### Phase 3 — AA→CG projection + distribution scoring (week 3)
- AA→CG projector (mass-weighted center-of-mass per bead per frame), backed by MDAnalysis.
- Distribution scorer: bond/angle histograms; KL-divergence-vs-Gaussian, skew, kurtosis, GMM-component count (BIC).
- Rule checker: R/S/T sizing, functional-group integrity, symmetry, bead-count plausibility.
- **Deliverable**: `agent/project.py`, `agent/score.py`. Reproduces a sensible score on the AA-mapped octanol fixture.

### Phase 4 — QA + repair loop (week 4–5, the central research piece)
- Repair strategies, in priority order: relabel (R↔S↔T), re-fragment by merging adjacent under-filled beads, re-fragment by splitting over-loaded beads, LLM-in-the-loop for chemistry-judgment ties.
- Stall recovery: subprocess timeout on AutoMARTINI3; on timeout, fall back to fragment-and-assemble (BRICS fragmentation → AutoMARTINI3 per fragment → reconnect using bond inference from AA reference).
- **Deliverable**: `agent/qa.py`, `agent/repair.py`. Regression: PMETAC produces a Martini-3-rule-compliant mapping; PSBMA returns *something* within a 5-min budget.

### Phase 5a — MCP server + skill packaging (week 6)
- Wrap Phases 2–4 internals as an MCP server (stdio transport, Python MCP SDK), runs inside `autom3` conda env.
- Author `SKILL.md` (Claude Code format) + `AGENTS.md` (mirror) with the trigger surface, tool listing, and decision-rule prose.
- Scripts under `skill/scripts/` for deterministic steps.
- Calibrate `program.md` thresholds + decision rules from Phase 0–4 fixture data.
- **Deliverable**: `mcp_server/` installable + `skill/` directory + `program.md` running end-to-end on the octanol fixture inside Claude Code.

### Phase 5b — Portability verification (week 6)
- Mount the same MCP server in a second agent runtime (Codex CLI is the cheapest target; OpenCode / Cursor as stretch).
- Run `program.md` on the octanol fixture in that runtime; assert byte-identical mapping output and equivalent score report.
- Document any per-runtime mounting steps in `skill/AGENTS.md`.
- **Deliverable**: portability test report in `tests/portability/` showing identical output across two agents.

### Phase 6 — Demonstration (week 7)
- Run the full protocol on Chris's PMETAC and PSBMA AA trajectories with no human intervention.
- Generate a write-up showing: mapping diff vs. AutoMARTINI3 baseline, bond/angle distributions before/after repair, provenance log of repair decisions.
- **Deliverable**: `demos/{pmetac,psbma}/` with the agent's full session log + final mapping + score report. Headline figure for any paper / write-up.

### Open research questions (track, not blockers)
- How formally to encode Martini 3 bead-selection rules so an LLM can reason over them?
- Charged groups (Qd/Qa, ammonium/sulfonate) — single bead vs split? PSBMA is the canonical case.
- Polymers: monomer-only mapping vs requiring a dimer/trimer reference?
- Confidence scoring: how does the agent communicate "sketchy" vs "solid"?

### Risks / known traps
- **AutoMARTINI3 stalls** — hard timeout + fallback non-negotiable.
- **`--fpred` always-on** — ALOGPS lookups fail on charged/exotic; default it in our wrapper.
- **AA sampling adequacy** — under-converged AA trajectories look spuriously non-Gaussian. Need a sampling-sufficiency check before scoring.
- **Polyply input contract** — monomer `.itp` shape differs from free-molecule `.itp`; verify before pass-through.

---

## 5. Status log

Section numbers below refer to §3 architecture; phases to §4.

| date       | milestone                                                  | status |
|------------|------------------------------------------------------------|--------|
| 2026-04-30 | repo created, plan synthesized from email chain            | done   |
| 2026-04-30 | Obsidian project note created                              | done   |
| 2026-04-30 | scope revised (drop CG sim; AA-projection scoring instead) | done   |
| 2026-04-30 | first-cut `program.md` drafted                            | done   |
| 2026-05-01 | architecture committed: MCP server + portable skill + program.md | done   |
| 2026-05-01 | Phase 0: `autom3` env + AutoMARTINI3 M3 (commit 1fff05a) installed | done |
| 2026-05-01 | Phase 0: octanol golden fixture captured                   | done   |
| 2026-05-01 | Phase 0: PMETAC + PSBMA fixtures captured (M3 fork; behavior diverges from email predictions — see `tests/fixtures/README.md`) | done |
| —          | Phase 0: ingest first AA trajectory from collaborators     | blocked on data |
| —          | Phase 0: confirm PSBMA M3 mapping with Chris/Seonghan      | needs review |
| —          | Phase 1: backends.md                                       |        |
| —          | Phase 2: classifier + dispatcher                           |        |
| —          | Phase 3: AA→CG projection + scoring                        |        |
| —          | Phase 4: QA + repair loop                                  |        |
| —          | Phase 5a: MCP server + skill packaging                     |        |
| —          | Phase 5b: portability verification on second runtime       |        |
| —          | Phase 6: demonstration on PMETAC + PSBMA                   |        |
