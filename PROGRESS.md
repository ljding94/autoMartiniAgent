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
| 2026-06-09 | **Sampling gap**: `step5_200_center.xtc` contains only **10 frames** (20 ns spacing for a 200 ns run). 190 bond + 180 angle observations is far below convergence for Gaussian-fit RMSE. Need a denser xtc — either re-`trjconv -dt` from Seonghan's original or ask him for a higher-frequency dump. Pre-requisite to the scorer (#3 evaluation half). | ~~open~~ → resolved 2026-06-30 via NAMD production data |
| 2026-06-09 | autom3 env additions: `MDAnalysis 2.9.0`, `pytest 9.0.3` (+ transitive scipy/matplotlib). | done |
| 2026-06-30 | **AA PEO-20 production data delivered** — `reference/peo20_solu/namd/`: CHARMM-GUI Polymer Builder solvated box (12 208 atoms; same 142-atom S1P1/ETHOX polymer at indices 1-142 → `PEO20_mapping.json` works as-is), 200 NAMD production segments (`step5_{1..200}.dcd`, 10 frames @ 100 ps each), **2000 frames over 200 ns = 200× the previous GROMACS dump**. Closes the sampling gap from 2026-06-09. | done |
| 2026-06-30 | Projector CLI extended: `--aa-traj` is now `nargs='+'` and natural-sorts a single quoted glob (`step5_2.dcd` before `step5_10.dcd`). Library signature widened to `Sequence[str|Path]`. Projected the full 200-segment NAMD trajectory → `derived/PEO20_solu/PEO20_cg.{xtc,gro}` (2000 frames × 20 beads). | done |
| 2026-06-30 | Chris's `peo_validation` branch (PR #1, `8f21c0a`) lifted into production at **`agent/score.py`** — Polyply `.itp` parser (bonds + angles + targets), Gaussian fitter, JSON report, PDFs. Reads canonical Martini target μ from `.itp`, emits Δ vs target per term. Bonds + angles only (per scope decision — torsions deferred as non-Gaussian). 6 new tests at `tests/test_score.py` (all green); 11/11 total. | done |
| 2026-06-30 | **PEO-20 numbers on the full trajectory** (`agent.score` against Polyply `.itp`, end-exclude 2): bond μ_fit = **0.3301 nm** (target 0.36, **Δ = −0.033**), σ_fit = 0.0114, RMSE = 0.79; angle μ_fit = **129.86°** (target 123, **Δ = +7.8°**), σ_fit = 13.07, RMSE = 0.0012. 30 000 bond obs, 28 000 angle obs. The 10-frame baseline (μ_bond=0.326, μ_angle=131.3) and 2000-frame production agree on the means within ~0.4%, so the original Δ signals are reproducible and the underlying chemistry — not the projector or the sampling — drives them. | done |
| 2026-06-30 | autom3 env additions: `mdtraj 1.10.3` (+ netCDF4 / cftime). | done |
| 2026-07-07 | **Second molecule: PSBMA-20mer** — AA data delivered at `reference/PSBMA_20mer/` (99 MB xtc gitignored; the ~50 KB solute-only gro stays tracked). **25 001 frames over 500 ns @ 20 ps spacing** — 250× denser than the PEO-20 baseline. Zwitterionic methacrylate polymer with **6 beads per monomer** (SC1 backbone + TN5a/TP2a esters + Q1± ammonium/sulfonate + SC1 propyl linker), 782 AA atoms, ATRP terminator (`BR1` on residue 20). | done |
| 2026-07-07 | **Scorer refactor — multi-type support** at `agent/score.py`. Bonds now grouped by `(b0, kb)`, angles by `(θ0, kₐ)`, each group fit + scored independently. `.itp` parser also reads `[atoms]` so bead-type labels appear in each group's title. PDFs go multi-panel (one panel per group). ScoreReport schema: `bond` → `bond_terms: list[TermStats]`, same for angles. PEO regression unchanged (μ_bond=0.3301, μ_angle=129.86). 12/12 tests pass. | done |
| 2026-07-07 | **Mapping + topology builders for PSBMA-20mer**: `scripts/build_psbma20_mapping.py` derives the 120-bead AA→CG mapping (static heavy-atom table from PSBA.itp SMILES annotations + spatial-nearest H assignment on frame 0; end groups H6/BR1 fold into their monomer's bead 1). `scripts/build_psbma20_itp.py` tiles the single-monomer PSBA.itp × 20 + adds 19 inter-monomer SC1-SC1 backbone bonds (b0=0.27 nm, kb=7500 — Martini-3 methacrylate default). 782/782 AA atoms assigned; per-bead masses identical across mid-chain monomers. | done |
| 2026-07-07 | **PSBMA-20mer first numbers** (end_exclude=1, ~475k obs per bond group / ~450k per angle group). Bonds: all six groups fit clean Gaussians; Δ vs auto-martini target ranges from −0.036 nm (`TN5a-TP2a`) to +0.019 nm (`SC1-TN5a`). Backbone SC1-SC1 (our chosen default b0=0.27) came out at 0.255 nm, Δ=−0.001 — the default was well-picked. **Angles: systematically off**, Δ ranges from −32° (`SC1-Q1-SC1`) to +30° (`Q1-TP2a-SC1`), and several σ's are 25-45° (i.e. bimodal / non-Gaussian). This is exactly the signal the tests/fixtures/psbma README flagged as "surprise success (needs chemistry review)" — auto-martiniM3's angle parameterization for PSBMA doesn't survive AA validation. | done |
| 2026-07-07 | **Goodness-of-fit metric: RMSE → cross-entropy.** Amplitude-domain RMSE scales with peak density, so it's not comparable across terms (bond RMSE ~0.8 vs angle RMSE ~0.001 for the same PEO chemistry). Cross-entropy H(measured, fitted-Gaussian) in nats compares normalized probability distributions directly and puts a large penalty when measured has mass where the fitted Gaussian has near-zero density (bimodal / non-Gaussian). PEO baselines: H_bond=4.37, H_angle=4.24 (both near the minimum for our bin resolution). PSBMA bond H ranges 3.7–6.6; angle H ranges 4.8–6.1. Same numeric scale across all terms → per-group scoring is now cross-comparable. | superseded |
| 2026-07-07 | **Goodness-of-fit metric: cross-entropy → KL → R².** Cross-entropy carries the measured distribution's own entropy H(P), so it's neither comparable across groups nor binwidth-invariant; switched to forward KL (H − H(P)). But KL *inverted* vs visual fit quality on PSBMA: clean narrow peaks scored ~2 nats while broad bimodal blobs scored ~0.1. Root cause (verified): forward KL over the fixed [0,180°]/[0,0.8nm] window is dominated by eps-clipped tail bins — a narrow σ≈6° Gaussian leaves ~160 of 180 bins near-zero where the measured histogram still has trace mass, each contributing `p·log(p/eps)`. Every mass/overlap metric (KL, JS, Hellinger, Bhattacharyya) shares this flaw. **Settled on R² in density space** = "how well the fitted curve traces the histogram" — weights each bin's residual equally, ignores rare tails, matches eye/physics. PSBMA angle ranking under R²: clean SC1-TN5a-TP2a 0.990, bimodal SC1-Q1-TP2a 0.787, worst broad blob 0.464. `TermStats.fit_r2` replaces `fit_cross_entropy`; plots/CLI/JSON all say R². | done |
| 2026-07-08 | **`end_exclude` default 2 → 0 (was silently dropping one monomer).** Collaborator flagged PSBMA n_obs short by ~25000 (one instance/frame = one monomer). Cause: `_drop_end_excluded` trims by **bead index**, and `end_exclude=2` (calibrated for PEO where 1 bead = 1 monomer) excludes beads {0,1,118,119} — which for PSBMA's 6-bead monomers chops a *different single monomer* off each term (monomer 1 for early-bead bonds, monomer 20 for late-bead bonds, both for backbone), giving inconsistent 18/19/20-per-frame counts. Changed default to 0 (keep all monomers — a validation tool shouldn't silently drop ~5% of data); trimming stays opt-in via `--end-exclude`, with help noting it counts beads not monomers (use a multiple of beads-per-monomer to trim whole ends). Now intra-monomer terms = 20/frame (500020), backbone = 19/frame (475019, correct: 19 links between 20 monomers). PEO essentially unchanged by including its end monomers (bond μ 0.3301→0.3304, angle 129.86→129.72, R² flat), confirming the trim was costing data for no benefit. Test fixture pins `end_exclude=2` so PEO regression assertions are unaffected. 15/15 pass. | done |
| 2026-07-08 | **Projector PBC bug fixed — spurious ~10° angle peak.** All PSBMA angle distributions carried a small peak near 10° (three bonded beads can't fold that sharply). Root cause: `agent/project.py` computed each bead's `center_of_mass()` on **raw coordinates**, so a bead whose AA atoms straddle a periodic boundary collapsed to a garbage mid-box point — producing impossible ~2.5 nm "bonds" and fake small angles. It bit **52.5%** of frames (this 20-mer chain crosses the boundary constantly); the AA universe carries no bonds, so MDAnalysis `unwrap=True` was unavailable. Fix: `_bead_com()` shifts each atom to its minimum image relative to the bead's first atom (via `minimize_vectors`, triclinic-safe) before the mass-weighted average — no bond info needed, no-op for whole beads. Re-projected PSBMA (47 s): spurious <20° observations 1.9–4.6% → **0.000%**, every angle R² up (backbone 0.871→0.939, `Q1-SC1-Q1` 0.963→0.986), bonds cleaner too; the defined angle's Δ moved −27.3°→−24.6° once the garbage tail was gone. The `SC1-Q1-TP2a` bimodality **survived** → it's a real two-rotamer feature, not an artifact. PEO unaffected (its AA traj was pre-centred/whole). New `test_bead_com_unwraps_across_pbc`; 15/15 tests pass. | done |
| 2026-07-14 | **#4 Loop — first cut: agent-driven mapping repair.** Built the mapping-edit verbs (`agent/repair.py`: `reassign_atoms` / `merge_beads` / `split_bead` on a `MappingState`, each re-deriving masses + heavy-atom counts, maintaining the bead-bond graph, and emitting a parameter-free `.itp`) plus the loop evaluator (`agent/evaluate.py`: project → score → scalar **Gaussianity objective** = `mean(1−R²)` over all bonds + all bonded angles; content-hash cached; `frame_stride` for fast interactive eval — stride-25 ≈ 1000 frames in ~2 s vs ~47 s full). Driver = **LLM-in-the-loop** (agent reads the score report, proposes edits, re-evaluates, keeps improvements). Objective is deliberately **target-free** — pure distribution Gaussianity, matching "make the measured CG distributions as Gaussian as possible", not Δ-vs-force-field. **PSBMA-20 result** (full 25 001 frames): baseline **0.0400 → 0.0352 (−12 %)** via a rule-valid edit that shifts each sidechain's bead windows one carbon (ester → `-O-CH₂-CH₂-` = O2,C5,C6; ammonium → `N⁺(CH₃)₂-CH₂` = N,C7,C8,C9; propyl → C10,C11), re-centering the Q1 charge on N⁺: `Q1-TP2a` bond R² **0.892→0.996**, bimodal `SC1-Q1-TP2a` angle **0.840→0.894** (cost: `Q1-SC1 #2` bond 1.000→0.863, all beads stay ≤4 heavy). The *unconstrained* optimum (C9→ammonium only) hits 0.0337 (−19 %) but makes a 5-heavy Q1 bead → cleanly surfaces the **Gaussianity-vs-Martini-sizing tension** the rule checker (not yet built) will arbitrate. Deliverable at `derived/PSBMA20/repair/` (repaired mapping + itp, before/after PDFs, `repair_provenance.json`). New `tests/test_repair.py` (trajectory-free); **26/26 pass**. Still pending: `qa.py` plateau/budget + acceptance gate, force-field param regeneration for merged/split beads (downstream), MCP packaging. | done |
| 2026-07-07 | **Angle scoring restricted to bonded angles + coverage audit.** The auto-martiniM3 PSBMA `.itp` defines 10 angles/monomer but only `1-2-3` (SC1-TN5a-TP2a) is a classical bond-bond angle; the other 9 are non-adjacent "structural" restraints (vertex not bonded to both arms) that a single harmonic/Gaussian can't model. `_group_angles` now filters to vertex-bonded-both-arms only → PSBMA drops 10→1 scored angle group (PEO unaffected: its consecutive angles all survive). **Always-on `angle_coverage` audit** (topology-only) added to `ScoreReport`: enumerates every consecutive-bond angle bead-type type and flags which the `.itp` defines vs omits — PSBMA reports "6 types, itp defines 1, missing 5". New `--all-bonded-angles` mode measures the omitted ones from the AA-mapped reference (writes `*_all_angles.{json,pdf}`, canonical report untouched). Measured missing angles: most are good harmonic candidates (R² 0.87–0.96) except `SC1-Q1-TP2a`/3-4-5 (bimodal, 0.787); `Q1-TP2a-TN5a`/2-3-4 sits at 160.7° (near-linear, harmonic-risky). Even the one defined angle's θ₀ is 25° off the AA reference (97.5° vs itp 122.2°) → the itp angle block needs review; deferred as parameter-fitting (out of scope). 14/14 tests pass. | done |
| —          | **Scope decision**: drop PDMAEMA / self-generate / broaden to Chris's charged 20-mers | open |
| —          | **#1 Stage**: sampling-sufficiency check                   | not started |
| —          | **#1 Stage**: build `scripts/aa_prep/` if no traj provided | contingency |
| —          | **#2 Process**: classifier + dispatcher (`classify.py`, `dispatch.py`) | not started |
| —          | **#2 Process**: Martini Mapper interactive-CLI wrapper     | not started |
| —          | **#2 Process**: cold-start fallback chain (BRICS + naive)  | not started |
| —          | **#2 Process**: Martini 3 rules table + lookup             | not started |
| —          | **#3 Evaluation**: AA→CG projector (`project.py`)          | done (2026-06-09, extended 2026-06-30) |
| —          | **#3 Evaluation**: scorer (`score.py`) — bonds + angles via Polyply `.itp` target | done (2026-06-30) |
| —          | **#3 Evaluation**: Martini-rule checks (R/S/T sizing, Q-bead defaults) | not started |
| —          | **#4 Loop**: mapping-edit verbs (`repair.py`) + Gaussianity-objective evaluator (`evaluate.py`) | done (2026-07-14, LLM-in-the-loop) |
| —          | **#4 Loop**: QA acceptance gate + plateau/budget tracker (`qa.py`) | not started |
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
