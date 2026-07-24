# autoMartiniAgent

An agent-driven workflow that turns an **all-atom (AA) simulation** of a
molecule into a **Martini 3 coarse-grained (CG) mapping** — i.e. which
AA atoms get grouped into which CG bead, and what bead type each one is.

The deliverable is the mapping itself plus a score report — we do not
run CG simulations. Validation is in-the-loop: project the AA trajectory
through a proposed mapping, look at the resulting bead–bead distance
and angle distributions, and adjust the mapping until they look like
clean Martini ones.

## The pipeline (four ingredients)

```
  #1 Stage              #2 Process            #3 Evaluation         #4 Agent loop
  -----------           -----------           -----------           -----------
  AA trajectory   ───▶  initial CG     ───▶  project AA      ───▶  score report
  + structure          mapping guess         through mapping       drives mapping
                       (Polyply for           → CG trajectory       repair, until
                        polymers,             → bond/angle          a target score
                        AutoMARTINI3 +        distributions         or budget hit
                        Martini Mapper        → Gaussian-fit
                        for small mols)       RMSE + Martini
                                              rule check
```

See [`PROGRESS.md`](PROGRESS.md) for the full plan and status log.

## What works today

| ingredient | what's built | path |
|---|---|---|
| #1 Stage | AA PEO-20 reference data (CHARMM36, TIP3P, **NAMD 200 ns @ 100 ps → 2000 frames**) | `reference/peo20_solu/namd/` |
| #1 Stage | AA PSBMA-20mer reference data (**500 ns @ 20 ps → 25 001 frames**, zwitterionic methacrylate) | `reference/PSBMA_20mer/` |
| #2 Process | Polyply CG topology for PEO-20 (20 × SN3r) + canonical AA→CG mapping | `reference/polyply_PEO20/` |
| #2 Process | Auto-martiniM3 CG topology + AA→CG mapping tiled to PSBMA-20mer (120 beads: SC1 / TN5a / TP2a / Q1± / SC1 / Q1± per monomer + inter-monomer backbone bonds) | `reference/PSBMA_20mer/PSBMA20.itp`, `PSBMA20_mapping.json` |
| #2 Process | mapping + topology regenerators | `scripts/build_peo20_mapping.py`, `scripts/build_psbma20_mapping.py`, `scripts/build_psbma20_itp.py` |
| **#3 Evaluation** | **AA→CG trajectory projector** (library API + CLI; accepts a list/glob of trajectory files) | `agent/project.py` |
| **#3 Evaluation** | **CG mapping scorer** — bonds grouped by `(b0, kb)`, angles by `(θ0, kₐ)`, per-group Gaussian fit scored by **R²** (density-space, tracks visual/physical fit quality) + Δ vs target, multi-panel PDFs + JSON. Only **bonded** angles (vertex bonded to both arms) are scored; an always-on **angle-coverage audit** flags bonded angles the itp omits, and `--all-bonded-angles` measures every one of them | `agent/score.py` |
| **#4 Loop** | **Mapping-edit verbs** — `reassign_atoms` / `merge_beads` / `split_bead` on a `MappingState` (re-derive masses + heavy counts, maintain the bead-bond graph, emit a parameter-free `.itp`) | `agent/repair.py` |
| **#4 Loop** | **Repair evaluator** — project → score → scalar **Gaussianity objective** `mean(1−R²)` over all bonds + all bonded angles (target-free); content-hash cached, `frame_stride` for fast interactive eval. Driver mode is **LLM-in-the-loop**: the agent reads the report, proposes an edit, re-evaluates | `agent/evaluate.py` |
| tests | projector COM verification (incl. PBC-unwrap) + scorer end-to-end + repair verbs/objective (26 tests) | `tests/test_project.py`, `tests/test_score.py`, `tests/test_repair.py` |

What's **not** built yet: the molecule classifier and backend dispatcher
(#2), the Martini-rule checks (R/S/T sizing, Q-bead defaults) that would
let the loop *auto-reject* size-violating edits (#3), the autonomous loop
controller (`agent/loop.py`, #4), and the MCP + skill packaging
(portability layer). The repair loop's edit verbs and Gaussianity objective
(#4) now exist and are driven **LLM-in-the-loop**. Next is `agent/loop.py`:
a minimal, **SWE-agent-style** Python controller — the `repair.py` verbs are
the curated action space, the `evaluate.py` score report is the feedback
observation, and a swappable policy (LLM now, a deterministic baseline for
the agentic-vs-deterministic ablation) drives keep-best-valid iteration until
a budget/plateau/`submit` stop. See [`PROGRESS.md`](PROGRESS.md) §4.

## Try it (PEO-20, end-to-end)

Set up:

```sh
conda activate autom3      # env where autoMartini3, MDAnalysis, polyply live
```

Regenerate the AA→CG mapping from the two `.itp` files:

```sh
python scripts/build_peo20_mapping.py
# → reference/polyply_PEO20/PEO20_mapping.json   (canonical, atom-index keyed)
# → reference/polyply_PEO20/PEO20.map            (Martini-style mirror)
```

Sanity-check it on a single frame:

```sh
python scripts/check_peo20_mapping.py
# prints per-bead masses + 1-2 bond and 1-2-3 angle stats
```

Project the full AA trajectory (200 NAMD segments → 2000 CG frames) through the mapping:

```sh
python -m agent.project \
  --aa-top   reference/peo20_solu/namd/step3_input.psf \
  --aa-traj  'reference/peo20_solu/namd/step5_*.dcd' \
  --mapping  reference/polyply_PEO20/PEO20_mapping.json \
  --out-dir  derived/PEO20_solu
# → derived/PEO20_solu/PEO20_cg.xtc   (2000-frame CG trajectory)
# → derived/PEO20_solu/PEO20_cg.gro   (single-frame CG reference structure)
```

Score the mapping (Gaussian fit + Δ vs Polyply target):

```sh
python -m agent.score \
  --itp        reference/polyply_PEO20/PEO20.itp \
  --cg-struct  derived/PEO20_solu/PEO20_cg.gro \
  --cg-traj    derived/PEO20_solu/PEO20_cg.xtc \
  --out-dir    derived/PEO20_solu
# → derived/PEO20_solu/score_report.json   (includes an `angle_coverage` audit)
# → derived/PEO20_solu/bond_hists.pdf
# → derived/PEO20_solu/angle_hists.pdf
```

Every run prints an **angle-coverage** line — which consecutive-bond angle
types the topology has vs. which the `.itp` actually defines as harmonic terms —
so bonded angles missing from the force field surface automatically. To measure
the distributions of *all* bonded angles (including the ones the `.itp` omits),
add `--all-bonded-angles`; it writes `score_report_all_angles.json` and
`*_all_angles.pdf` alongside (leaving the canonical report untouched), reporting
each angle's fitted μ/σ/R² and a Δ-vs-target only where the `.itp` supplies one.

## Repair the mapping (agent loop, #4)

The loop's objective is a single **target-free Gaussianity error** — how far the
AA-projected distributions are from single Gaussians, averaged over every bond and
every bonded angle:

```
error = mean over (all bonds + all bonded angles) of (1 − R²_fit)      # 0 = ideal
```

Evaluate any mapping's error directly (use `--frame-stride` for a fast, lower-res
pass while iterating):

```sh
python -m agent.evaluate \
  --mapping reference/PSBMA_20mer/PSBMA20_mapping.json \
  --itp     reference/PSBMA_20mer/PSBMA20.itp \
  --aa-top  reference/PSBMA_20mer/PSBMA_20mer_no_water.gro \
  --aa-traj reference/PSBMA_20mer/PSBMA_20mer_no_water_skip10.xtc \
  --frame-stride 25
# objective (mean 1-R^2) = 0.0417 ... worst: SC1-Q1-TP2a, Q1-TP2a, ...
```

The agent then edits the mapping with the `agent.repair` verbs (`reassign_atoms`,
`merge_beads`, `split_bead`) and re-evaluates, keeping edits that lower the error.
Worked example on **PSBMA-20** (full 25 001 frames): shifting each sidechain's bead
windows one carbon — ester → `-O-CH₂-CH₂-`, ammonium → `N⁺(CH₃)₂-CH₂`, propyl → `CC`
— re-centers the Q1 charge on N⁺ and drops the error **0.0400 → 0.0352 (−12 %)**
(the broad `Q1-TP2a` bond R² 0.892→0.996, the bimodal `SC1-Q1-TP2a` angle
0.840→0.894), while keeping every bead ≤ 4 heavy atoms. The repaired mapping, its
itp, before/after PDFs and a provenance log land in `derived/PSBMA20/repair/`.
(The *unconstrained* optimum reaches −19 % but needs a 5-heavy bead — the
Gaussianity-vs-Martini-sizing tension a rule checker will arbitrate.)

Run the tests:

```sh
python -m pytest tests/ -v
```

## What the mapping artifact looks like

`PEO20_mapping.json` is the source of truth — atom-index keyed, easy for
the projector to consume:

```jsonc
{
  "molecule": "PEO20",
  "weighting": "mass",
  "beads": [
    {
      "bead_id": 1,
      "bead_name": "EC",
      "bead_type": "SN3r",
      "aa_residue": 1,
      "atom_indices": [1, 2, 3, 4, 5, 6, 7, 8],
      "atom_names": ["O1", "HO1", "C1", "H11", "H13", "C2", "H21", "H23"],
      "aa_mass_sum": 45.061,
      "comment": "HO terminus"
    },
    // ... 18 more mid-chain ether monomers ...
    { "bead_id": 20, "comment": "CH3 terminus", ... }
  ]
}
```

`PEO20.map` is the same mapping in Martini-style format (atom-name per
residue) for human inspection.

## First numbers

Goodness-of-fit is reported as **cross-entropy** `H(measured, fitted-Gaussian)`
in nats — H → self-entropy of the measured histogram when the fit is perfect
(a few nats for our bin resolution), and grows when the measured distribution
has mass where the fitted Gaussian is ≈ 0 (bimodal / non-Gaussian).

**PEO-20** (single bead type: SN3r; 2000 frames @ 100 ps; end_exclude=2):

| term | Gaussian fit μ | fit σ | H (nats) | target | **Δ vs target** |
|---|---:|---:|---:|---:|---:|
| bond `SN3r-SN3r` (nm) | 0.3301 | 0.0114 | 4.37 | 0.360 | **−0.033** |
| angle `SN3r-SN3r-SN3r` (°) | 129.86 | 13.07 | 4.24 | 123 | **+7.8** |

**PSBMA-20mer** (6 beads / monomer; 25 001 frames @ 20 ps; end_exclude=1):

Bond groups (auto-martiniM3 targets from `PSBA.itp`, inter-monomer backbone
tiled with Martini-3 methacrylate defaults b0=0.27 / kb=7500):

| bond group | μ_fit (nm) | σ | H (nats) | target | Δ | n_obs |
|---|---:|---:|---:|---:|---:|---:|
| `SC1-TN5a` (backbone→carbonyl) | 0.253 | 0.008 | 5.01 | 0.24 | **+0.019** | 475 019 |
| `TN5a-TP2a` (carbonyl→ester O) | 0.198 | 0.005 | 3.71 | 0.24 | **−0.036** | 500 020 |
| `Q1-SC1 #1` (propyl-sulfonate) | 0.270 | 0.007 | 5.18 | 0.26 | +0.016 | 475 019 |
| `SC1-SC1` (inter-monomer backbone) | 0.255 | 0.011 | **6.64** | 0.27 | −0.001 | 450 018 |
| `Q1-TP2a` (ester→ammonium) | 0.294 | 0.010 | **6.20** | 0.31 | −0.019 | 500 020 |
| `Q1-SC1 #2` (ammonium→propyl) | 0.305 | 0.005 | 4.78 | 0.33 | −0.024 | 500 020 |

Angle groups tell a much sharper story — many groups have Δ > 25° or σ > 25°,
so the single-Gaussian fit is really a goodness-of-fit metric that flags
"auto-martiniM3's angle parameterization for this zwitterion doesn't survive
AA validation". This is exactly the signal `tests/fixtures/psbma/README.md`
foreshadowed as "surprise success (needs chemistry review)".
See `derived/PSBMA20/{bond,angle}_hists.pdf`.

The projection itself is verified bit-exact against a hand-coded COM in
`tests/test_project.py`, so the Δ's above are real chemistry signals — not
projector or sampling artifacts.

## Repo layout

```
autoMartiniAgent/
├── agent/                      # core agent code
│   ├── project.py              # #3 — AA→CG trajectory projector
│   ├── score.py                # #3 — Gaussian-fit bond/angle scorer vs .itp target
│   ├── repair.py               # #4 — mapping-edit verbs (reassign/merge/split)
│   └── evaluate.py             # #4 — project→score→Gaussianity-error loop evaluator
├── scripts/                    # standalone helpers
│   ├── build_peo20_mapping.py  # build PEO20_mapping.json + PEO20.map
│   └── check_peo20_mapping.py  # single-frame mapping sanity check
├── reference/                  # inputs (committed)
│   ├── gromacs/                # AA PEO-20 (legacy 10-frame GROMACS dump)
│   ├── peo20_solu/             # AA PEO-20 production: NAMD, 2000 frames, solvated
│   ├── polyply_PEO20/          # PEO-20 CG topology + AA→CG mapping
│   ├── PSBMA_20mer/            # AA PSBMA-20mer 25k-frame trajectory (.xtc gitignored) + CG topology + mapping
│   └── email_chain.md          # April 2026 ORNL conversation seed
├── validation/                 # ad-hoc analysis notebooks/scripts (e.g. Chris's PR)
│   └── PEO20/                  # superseded by agent/score.py — kept as worked example
├── tests/                      # pytest suite
├── vendor/                     # third-party backends (gitignored)
│   └── Automartini_M3/         # logP-based small-molecule mapper
├── derived/                    # runtime outputs (committed as reference snapshots; regeneratable from inputs)
├── PROGRESS.md                 # full plan, scope, status log
├── program.md                  # agent-runnable protocol (draft)
└── README.md                   # this file
```

## Background

Triggered by an April 2026 ORNL email chain in which Seonghan Kim
circulated a reproducible SMILES → Martini 3 pipeline using
[Auto-MartiniM3](https://github.com/Martini-Force-Field-Initiative/Automartini_M3),
and Chris Walker stress-tested it on charged polymer monomers (PMETAC,
PSBMA). Two failure modes surfaced — subprocess stalls on zwitterionic
chemistry and Martini-3-rule-violating bead sizes — and motivated the
agent-driven QA + repair loop on top.

Collaborators (ORNL): Lijie Ding (driver), Seonghan Kim (Stage-1
pipeline, AA data), Chris Walker (chemistry validation, AA data), Jan
Michael Carrillo. Reference paper:
[Souza et al., *J. Chem. Inf. Model.* 2026](https://doi.org/10.1021/acs.jcim.5c02903).

## License

MIT.
