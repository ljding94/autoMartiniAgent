# autoMartiniAgent — Project Outline

**Goal**: an agent-runnable workflow that takes AA simulation data + chemical structure (SMILES / mol / pdb) and produces a validated Martini 3 mapping — atom-index → bead, bead type, bead size — plus a score report and provenance log. The mapping is the deliverable. We do not run CG simulations. Validation is in-the-loop: project the AA trajectory through the proposed mapping and score bead-bead distance and angle distributions.

Collaborators (ORNL): Lijie Ding (driver, dingl1@ornl.gov), Seonghan Kim (kimsn@ornl.gov), Chris Walker (walkercc@ornl.gov), Jan Michael Carrillo (carrillojy@ornl.gov).

---

## The four ingredients

The project closes a loop: **stage → propose → evaluate → repair → re-evaluate**. Each ingredient must exist for the loop to close.

### 1. Stage — atomistic reference data

The substrate the agent reasons over.

**Need**: AA trajectory + topology + chemical structure for each test molecule. Trajectory must be long enough that bead-bead distance and angle distributions are converged.

**Have**:
- **AA PEO-20 (delivered by Seonghan 2026-05-15, inspected 2026-05-19)** at `reference/gromacs/`:
  - 20 ETHOX residues, 142 atoms total; **asymmetric end groups** — residue 1 is hydroxyl terminus (`OG311+HGP1`, 8 atoms), residue 20 is ethyl-ether terminus (`CG331+3×HGA3`, 8 atoms). Net chain is **HO-CH₂-CH₂-O-(CH₂-CH₂-O)₁₈-CH₂-CH₃** (i.e. HO-(CH₂CH₂O)₁₉-CH₂-CH₃; *not* `-CH₂-CH₂-CH₃` — only one methylene between the last ether O and the methyl C).
  - **CHARMM36 General** force field (via CHARMM-GUI converter), **TIP3P** water (4022 molecules), no ions.
  - Production: `step5_200.xtc` / `step5_200_center.xtc` — 200 ns, 2 fs timestep, PME, Force-switch vdW (1.0/1.2 nm), NPT 300 K / 1 bar, h-bond LINCS.
  - Reference structure: `step4.0_minimization.gro` or `equil.tpr`.
- **AA 20-mer charged-polymer archive on Kronos (Chris, available not yet pulled)**: PMETAC, PMPC, PSBMA, PNOMA, P2VPPS — all 300+ ns, 8 atactic sequences averaged. *All out of current v1 scope (charged)* — could broaden scope to include them.
- Seonghan's pipeline (`reference/email_chain.md`) is **CG-only** — never to be confused with AA.
- CG references: `tests/fixtures/octanol/OCOL.itp` (golden), `tests/fixtures/psbma/PSBA.itp` (out-of-scope zwitterion baseline), `reference/polyply_PEO20/` (Polyply CG reference for PEO-20 — pair with AA above; now also holds the AA→CG mapping artifact, see #2 status entry 2026-06-09).

**Open gap**: PDMAEMA monomer AA **is not** in Chris's archive (he has PMETAC, the charged-methylated version, not the neutral PDMAEMA). v1 small-molecule headline test currently has no AA data. Three options on the table — see "Scope decision (open)" below.

**Chris's process notes (worth honoring)**:
- Convergence check: bonded distributions in 50 ns chunks; discard unequilibrated head.
- One atactic sequence is sufficient for bonds/angles. Torsions need averaging across multiple sequences → defer torsions in v1.

**Sampling-sufficiency check** is part of this ingredient: trajectory must show converged second moments before the scorer (#3) is allowed to call rule violations. Under-sampled AA looks spuriously non-Gaussian.

**Fallback if no AA trajectories arrive**: build a thin AA-prep pipeline ourselves (SMILES → LigParGen → GROMACS solvate with TIP3P → NPT eq → NVT production) under `scripts/aa_prep/`. Adds scope but makes the project self-contained.

### 2. Process — initial AA→CG mapping generation

The cold start: chemical structure → first-cut mapping, before any AA-driven scoring.

**v1 scope** (per Seonghan 2026-05-10): **neutral small molecules + polymers in the Polyply built-in library** only. Charged / zwitterionic molecules drop from v1 — both small-molecule backends fail on them (AutoMARTINI3 rejects via ALOGPS; Martini Mapper silently mis-types because it has no Q-bead dictionary, which is dangerous in an automated pipeline).

**Small-molecule backends — two, run in parallel for v1**:
- **AutoMARTINI3** (M3 fork, `vendor/Automartini_M3` @ `1fff05a`). logP-based strategy. ≤25 heavy atoms. Fully CLI. **Always invoke with `--fpred`**.
- **Martini Mapper**. SMARTS rule-based strategy. No heavy-atom cap. Currently interactive — needs a `pexpect`-style wrapper.

When both produce identical bead types → high confidence cold start. When they disagree (e.g., PDMAEMA gives same 3-bead grouping but conflicting types) → the disagreement is a **research signal** the refiner agent adjudicates using AA distributions.

**Polymer backend**: **Polyply**, restricted in v1 to its built-in Martini 3 library: PEO, PS, PMMA, PE, PVA, PDMS, PSS. Out-of-library polymers defer to v2 (would need monomer parameterization first). Input contract from Seonghan's email:
```
polyply gen_params -lib martini3 -seq <NAME>:<N> -o <NAME><N>.itp -name <NAME><N>
polyply gen_coords -p topol.top -o <NAME><N>.gro -name <NAME><N> -box X Y Z
# random-walk coords need: gmx grompp/mdrun (EM) → gmx trjconv -pbc mol -center
```

**Cold-start fallback chain** (inside #2, so the loop always sees iter 0 with *some* mapping):
1. AutoMARTINI3 (small) or Polyply (polymer).
2. On crash / cap-exceeded: BRICS-fragment-and-assemble (fragment with RDKit BRICS → AutoMARTINI3 per fragment → reconnect using bond inference).
3. On all-backends-fail: naive 1-bead-per-heavy-atom.

**Dispatcher** (our code, not yet written):
- `agent/classify.py` — RDKit heuristics: heavy-atom count, peptide-bond detection, charge (refuses if charged in v1), aromaticity → category ∈ {small, polymer}. Out-of-scope categories return a clear error.
- `agent/dispatch.py` — route to backend; emit uniform mapping schema regardless of which backend ran; tag provenance (which backend, which fallback step).

**Naming traps**: never use `MOL=OCT` (collides with octane in `martini_v3.0.0_solvents_v1.itp`); always check the solvents `.itp` before naming.

**Known backend behavior** (from Phase 0 + Seonghan's update):
- **Octanol**: both backends succeed; golden reference.
- **PDMAEMA** (`CC(C(=O)OCCN(C)C)C`, neutral monomer): both backends produce 3-bead mappings with identical atom groupings but **conflicting bead types**. v1 headline test case.
- **PEO-20**: Polyply built-in `gen_params` produces a 20-bead `SN3r` chain. Polymer-route reference.
- **PMETAC + PSBMA** (charged): out of v1 scope. Kept under `tests/fixtures/known_failures/` as smoke tests for "tool refuses / silently mis-types charged groups."

### 3. Evaluation — distribution-based mapping scorer

The signal the agent optimizes against. Replaces the CG-simulation validation step from Seonghan's original three-stage workflow.

**Projector** (not yet written): `agent/project.py` — read AA trajectory + topology + proposed mapping; compute mass-weighted COM per bead per frame; emit a CG trajectory. Backed by MDAnalysis.

**Scorer** (not yet written): `agent/score.py` — over the CG trajectory:
- **Terms scored** (per Lijie 2026-05-10):
  - Bonded distance distributions (1-2 bead pairs along the chain).
  - 1-3 bead distance distributions (next-nearest along the chain — captures angle stiffness as an independent check).
  - Bonded angle distributions (1-2-3 triples).
- **Per-term metric**: fit a Gaussian to the histogram (unit-area normalized), report **RMSE between histogram and fit**. Lower = closer to the harmonic Martini ideal.
- Martini-rule checker (binary): R/S/T sizing (R = 4, S = 3, T = 2 heavy atoms), functional-group integrity, symmetry, bead-count plausibility, no Q-beads with neutral types in places that should be charged (defensive).
- Output: `{bond_distributions, angle_distributions, gaussian_fit_rmse_per_term, rule_violations, scalar_score, backend_disagreement}`. The last is a flag from #2 noting whether AutoMARTINI3 and Martini Mapper agreed on bead types.

**Acceptance criterion** (per Lijie 2026-05-10):
- Inner loop **minimizes scalar score within budget**, not against a fixed RMSE threshold (an outer autoresearch loop can re-engage with more budget if needed).
- Stop conditions: (a) zero rule violations AND no improvement in scalar score over last *K* iterations (plateau), or (b) iteration budget hit (default 10), or (c) wall-clock budget hit (default 10 min). Returns the best mapping seen.

**Bootstrap without AA data**: the scorer can be developed and unit-tested against a synthetic harmonic CG system (project a known-good model, assert score ≈ Gaussian and rules-clean). #3 does not block on #1.

### 4. Agent loop — iterative repair

Closure of the loop: score report → mapping revision → re-score, until acceptance criteria met or budget exhausted.

**Two modes, selected via CLI flag at runtime** (per Lijie 2026-05-10):
- `--mode tight` — deterministic repair only (verbs: `relabel_size_class`, `merge_beads`, `split_bead`). LLM acts as referee on ties / on parsing the score report into a repair choice. Reproducible.
- `--mode loose` — adds LLM-driven structural verbs (`change_bead_type`, `reassign_atom`). The agent reads the score report + AA-derived chemistry context + Martini-rule prose and proposes structured edits. Higher flexibility, lower reproducibility.

**Action vocabulary** (loose mode = full set, tight mode = first three only):
- `relabel_size_class(bead_id, R|S|T)` — fix R/S/T sizing without changing bead count.
- `merge_beads(bead_ids[])` — combine adjacent under-filled beads.
- `split_bead(bead_id, into=[{atoms[]}, ...])` — split an over-loaded bead.
- `change_bead_type(bead_id, new_type)` — adjudicate bead-type conflict (the dominant verb when the two backends disagreed in #2; e.g., PDMAEMA case).
- `reassign_atom(atom_id, from_bead, to_bead)` — move an atom between beads.

Each verb emits a structured JSON edit logged to a provenance trail.

**QA + repair** (not yet written):
- `agent/qa.py` — interpret the score report; decide accept / repair / escalate. Includes the plateau detector and budget tracker.
- `agent/repair.py` — chooses + applies action verbs based on score report; tight mode picks deterministically, loose mode delegates to LLM with structured-JSON output guard.
- Stall / crash recovery handled in #2 (cold-start fallback chain), not here.

**Packaging** (not yet written):
- `mcp_server/` — Python MCP server (stdio), runs in `autom3` env. Exposes: `classify`, `propose_mapping`, `project_trajectory`, `score_mapping`, `repair_mapping`, `martini_rules.lookup`. Large data passed by path.
- `skill/SKILL.md` (Claude Code) + `skill/AGENTS.md` (mirror for Codex / OpenCode / Cursor / Continue) — same content, two filenames. When and how to invoke the tools, decision rules, recovery strategies.
- `program.md` (already drafted) — agent-agnostic protocol that drives the loop end-to-end. References MCP tool names only; no agent-specific syntax.

**Headline demonstration**: drop a new monomer + AA trajectory + `program.md` into any MCP-aware agent with the skill loaded → agent autonomously produces a validated Martini 3 mapping, no human in the loop.

---

## Architecture

Three portable layers:

```
autoMartiniAgent/
├── mcp_server/         # Python MCP server (stdio), runs inside autom3 conda env
│   ├── pyproject.toml
│   └── src/auto_martini_mcp/
├── skill/
│   ├── SKILL.md        # Claude Code skill format
│   ├── AGENTS.md       # cross-agent mirror
│   └── scripts/        # deterministic helpers callable without an LLM
├── program.md          # the autoresearch protocol
├── tests/fixtures/
└── PROGRESS.md
```

Portability invariant: nothing in `skill/` or `program.md` may reference Claude-Code-specific syntax, agent built-ins, or non-MCP tool names.

---

## Alignment with Seonghan's three-stage workflow

We adopt his Stage 1 as-is, reinterpret Stage 2, and absorb Stage 3 into the inner loop:

| | Seonghan's framing | Ours (v1) |
|---|---|---|
| **Stage 1 — initial CG generation** | AutoMARTINI3 / Polyply / Martinize2 → first-cut `.itp` | Same. v1 small route runs **AutoMARTINI3 + Martini Mapper in parallel** so the agent has a backend-disagreement signal to adjudicate. |
| **Stage 2 — refinement** | AA reference → bonded distribution targets → **BI or Swarm-CG** → refined `.itp` (parameter fitting) | AA reference → **AA→CG projection + Gaussian-fit RMSE scoring** → **agent-driven mapping repair** (atom-bead grouping + bead types). *Not* parameter fitting. |
| **Stage 3 — validation** | Separate CG simulation → compare with AA-mapped distributions → report | Folded into Stage 2's inner loop. Optional v2: add Seonghan's CG-sim as a forward check on the final mapping. |

**Critical clarification (2026-05-10)**: Seonghan's reproduction script (`run.sh` for octanol) is a **Stage 1 + CG simulation** pipeline — its `.mdp` and force-field includes are Martini, and the trajectory it produces is CG, not AA. We **cannot** use it to generate the AA reference trajectories Stage 2 requires. AA data must come from a separate atomistic pipeline (LigParGen / OPLS-AA + TIP3P + GROMACS production), either provided by collaborators or built under `scripts/aa_prep/`.

---

## Status log

Section labels below map to the four ingredients above.

| date       | milestone                                                  | status |
|------------|------------------------------------------------------------|--------|
| 2026-04-30 | repo created, plan synthesized from email chain            | done   |
| 2026-04-30 | scope revised (drop CG sim; AA-projection scoring instead) | done   |
| 2026-04-30 | first-cut `program.md` drafted                             | done   |
| 2026-05-01 | `autom3` env + AutoMARTINI3 M3 fork installed              | done   |
| 2026-05-01 | octanol golden fixture captured                            | done   |
| 2026-05-01 | PMETAC + PSBMA fixtures captured (behavior diverges from email) | done |
| 2026-05-10 | outline restructured around the four ingredients           | done   |
| 2026-05-10 | scope tightened (neutral-only) per Seonghan's testing update | done |
| 2026-05-10 | acceptance criterion + modes (tight/loose) + action vocabulary locked | done |
| 2026-05-10 | confirmed Seonghan's pipeline is CG-only; AA data must come separately | done |
| 2026-05-10 | AA trajectories requested from Chris + Seonghan (octanol, PDMAEMA, PEO-20) | sent |
| 2026-05-10 | Seonghan offered AA PEO-20 (single chain, water, no ions); ETA ~2026-05-15 | replied |
| 2026-05-10 | Master TODO checklist added to Obsidian note               | done   |
| 2026-05-15 | **AA PEO-20 delivered** to `reference/gromacs/` (CHARMM36, TIP3P, 200 ns) | done |
| 2026-05-19 | AA PEO-20 inspected; end-group asymmetry (HO / CH₃) noted   | done |
| 2026-05-19 | Chris back from vacation — archive lists charged 20-mers (PMETAC, PMPC, PSBMA, PNOMA, P2VPPS); no PDMAEMA monomer | known |
| 2026-06-09 | **PEO-20 AA→CG mapping derived** — 20 ETHOX residues → 20 SN3r beads (1:1, mass-weighted, end-groups fold into bead 1 / bead 20). Artifacts: `reference/polyply_PEO20/PEO20_mapping.json` (canonical, atom-index keyed) + `PEO20.map` (Martini-style mirror). Hand-derived from CHARMM36 ETHOX + Polyply's `PEO20.itp` — *not* extracted from Polyply's built-in `.mapping` library (which targets OPLS atom names, not CHARMM). | done |
| 2026-06-09 | Helper scripts added: `scripts/build_peo20_mapping.py` (regenerate mapping deterministically from the two `.itp` files) + `scripts/check_peo20_mapping.py` (project `equil.gro` through mapping; assert per-bead mass sums; print 1-2 bond + 1-2-3 angle stats) | done |
| 2026-06-09 | Single-frame sanity check on `equil.gro`: mean 1-2 bond 0.310 nm (Polyply r₀ = 0.36), mean 1-2-3 angle 118° (Polyply θ₀ = 123°). All 20 bead-mass sums match expected AA sums to 0.01 g/mol → atom-index grouping verified. Whether the 0.05 nm bond offset is mapping-driven or parameter-driven is the scorer's call once the full xtc is projected. | done |
| 2026-06-09 | PROGRESS.md end-group description corrected — terminus is `-O-CH₂-CH₃` (ethyl-ether), not `-O-CH₂-CH₂-CH₃`. Off-by-one CH₂ in original note. | done |
| 2026-06-09 | **#3 Evaluation — AA→CG projector landed** at `agent/project.py`. Generic across any mapping JSON in our schema; CLI + library API. Reads AA top (`.tpr`/`.gro`) + AA traj (`.xtc`/...), writes CG `.xtc` + single-frame CG `.gro`. Mass-weighted COM via MDAnalysis AtomGroup; mass-sum validation guards against atom-index off-by-ones. Smoke tests at `tests/test_project.py` (5 pass) check shape + bit-exact agreement with from-scratch Python COM. | done |
| 2026-06-09 | Projector exercised on `step5_200_center.xtc` → `derived/PEO20/PEO20_cg.{xtc,gro}`. 1-2 bond mean **0.326 nm** (Polyply r₀=0.36, Δ=-0.034) ; 1-2-3 angle mean **131.3°** (Polyply θ₀=123°, Δ=+8°). These are real signals the scorer will need to act on — *not* projector bugs (positions verified bit-exact vs Python COM). | done |
| 2026-06-09 | **Sampling gap**: `step5_200_center.xtc` contains only **10 frames** (20 ns spacing for a 200 ns run). 190 bond + 180 angle observations is far below convergence for Gaussian-fit RMSE. Need a denser xtc — either re-`trjconv -dt` from Seonghan's original or ask him for a higher-frequency dump. Pre-requisite to the scorer (#3 evaluation half). | open |
| 2026-06-09 | autom3 env additions: `MDAnalysis 2.9.0`, `pytest 9.0.3` (+ transitive scipy/matplotlib). | done |
| —          | **Scope decision**: drop PDMAEMA / self-generate / broaden to Chris's charged 20-mers | open |
| —          | **#1 Stage**: sampling-sufficiency check                   | not started |
| —          | **#1 Stage**: build `scripts/aa_prep/` if no traj provided | contingency |
| —          | **#2 Process**: classifier + dispatcher (`classify.py`, `dispatch.py`) | not started |
| —          | **#2 Process**: Martini Mapper interactive-CLI wrapper     | not started |
| —          | **#2 Process**: cold-start fallback chain (BRICS + naive)  | not started |
| —          | **#2 Process**: Martini 3 rules table + lookup             | not started |
| —          | **#3 Evaluation**: AA→CG projector (`project.py`)          | not started |
| —          | **#3 Evaluation**: scorer (`score.py`) — bootstrap on synthetic data first | not started |
| —          | **#4 Loop**: QA + repair (`qa.py`, `repair.py`)            | not started |
| —          | **#4 Loop**: MCP server + skill packaging                  | not started |
| —          | **#4 Loop**: portability check on second runtime           | not started |
| —          | Demonstration on PDMAEMA + PEO-20 (depends on #1)          | not started |
| —          | Move PMETAC + PSBMA fixtures under `tests/fixtures/known_failures/` | not started |

---

## What we can do without AA data

#1 is gating, but #2 and #3 can advance in parallel without it:

- **#2** is structure-only — classifier, dispatcher, two-backend wrapper for AutoMARTINI3 + Martini Mapper, BRICS fallback. PDMAEMA cold-start divergence (the two backends disagree on bead types) can be captured as a fixture without any AA trajectory.
- **#3 scorer** can be built against synthetic harmonic CG systems (project a known-good model, assert Gaussian-fit RMSE ≈ 0, rules clean). The projector half *does* need AA data for any non-trivial test.
- **#4** glues #2 and #3 together; once both have minimal versions running, the QA-repair loop can be wired and exercised on synthetic-AA + the PDMAEMA disagreement fixture.

Net: AA data unblocks the projector and the demonstration, but does not block the rest of the architecture.

---

## Background — email chain and key findings

The April 2026 email chain (Walker / Kim / Ding / Carrillo) seeded the plan. Two findings from Phase 0 reproduction altered the original premise:

1. **Chris's failure cases were on `auto_martini` (M2), not `auto_martiniM3`.** His PSBMA stall and PMETAC mis-sizing came from the older Bereau & Kremer tool. On the M3 fork: PSBMA converges in ~5 s with a rule-compliant mapping; PMETAC fails differently (disconnected-fragment intermediates). Implication: the PSBMA-stall fallback is deprioritized; the PMETAC disconnected-fragment crash is the new headline repair target.
2. **`--fpred` is mandatory** in our wrapper. ALOGPS fragment lookup fails on charged / exotic fragments otherwise.

Reference paper from Carrillo (not yet read): https://doi.org/10.1021/acs.jcim.5c02903 (J. Chem. Inf. Model., 2026).

Reproduction notes: `tests/fixtures/README.md`. Email source: `reference/email_chain.md`. Polyply worked example: `reference/polyply_PEO20/`.

---

## Open questions

- How to encode Martini 3 bead-selection rules so an LLM can reason over them — lookup table, rule prose, or example library?
- Charged groups (Qd / Qa, ammonium / sulfonate): single bead vs split? PSBMA is the canonical case; chemistry sign-off from Chris / Seonghan still pending.
- Polymers: monomer-only mapping sufficient, or is a dimer / trimer reference required for the dispatcher to do its job?
- Confidence scoring: how does the agent communicate "sketchy" vs "solid" mappings to a downstream user?

## Risks

- **AA sampling adequacy** — under-converged trajectories look spuriously non-Gaussian. The sampling-sufficiency check in #1 is non-negotiable before #3 trusts its own output.
- **`--fpred` always-on** — default it in our wrapper; never expose as optional.
- **Polyply input contract drift** — monomer `.itp` shape differs from free-molecule `.itp`; verify on the PEO20 fixture before relying on pass-through.
- **Backend stalls / crashes** — hard subprocess timeout + fragment-and-assemble fallback non-negotiable in the dispatcher.
