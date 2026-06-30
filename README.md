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
| #1 Stage | AA PEO-20 reference data (CHARMM36, TIP3P, **NAMD 200 ns @ 100 ps spacing → 2000 frames**) | `reference/peo20_solu/namd/` |
| #2 Process | Polyply CG topology for PEO-20 (20 × SN3r) + canonical AA→CG mapping | `reference/polyply_PEO20/` |
| #2 Process | mapping regenerator + sanity-check scripts | `scripts/build_peo20_mapping.py`, `scripts/check_peo20_mapping.py` |
| **#3 Evaluation** | **AA→CG trajectory projector** (library API + CLI; accepts a list/glob of trajectory files) | `agent/project.py` |
| **#3 Evaluation** | **CG mapping scorer** — Gaussian-fit bond + angle distributions, target μ pulled from Polyply ``.itp``, JSON report + PDFs | `agent/score.py` |
| tests | projector COM verification + scorer end-to-end | `tests/test_project.py`, `tests/test_score.py` |

What's **not** built yet: the molecule classifier and backend dispatcher
(#2), the Martini-rule checks (R/S/T sizing, Q-bead defaults) that the
scorer will eventually enforce (#3), the agent's repair loop (#4), the
MCP server, and the skill packaging.

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
# → derived/PEO20_solu/score_report.json
# → derived/PEO20_solu/bond_hists.pdf
# → derived/PEO20_solu/angle_hists.pdf
```

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

## PEO-20 first numbers

Projecting the 200 NAMD production segments (2000 frames @ 100 ps) through
the canonical PEO mapping and scoring against the Polyply ``.itp`` target:

| metric | Gaussian fit μ | fit σ | fit RMSE | Polyply M3 target | Δ vs target |
|---|---:|---:|---:|---:|---:|
| 1-2 bond (nm) | **0.3301** | 0.0114 | 0.79 | 0.360 | **−0.033** |
| 1-2-3 angle (°) | **129.86** | 13.07 | 0.0012 | 123 | **+7.8** |

(Both ends dropped by `--end-exclude 2` per Chris's convention; 30 000 bond
observations, 28 000 angle observations after exclusion.) The projection
itself is verified bit-exact against a hand-coded COM in
`tests/test_project.py`, so the Δ values are real chemistry signals — the
agent loop will use them to decide whether the mapping is acceptable.

## Repo layout

```
autoMartiniAgent/
├── agent/                      # core agent code
│   ├── project.py              # #3 — AA→CG trajectory projector
│   └── score.py                # #3 — Gaussian-fit bond/angle scorer vs .itp target
├── scripts/                    # standalone helpers
│   ├── build_peo20_mapping.py  # build PEO20_mapping.json + PEO20.map
│   └── check_peo20_mapping.py  # single-frame mapping sanity check
├── reference/                  # inputs (committed)
│   ├── gromacs/                # AA PEO-20 (legacy 10-frame GROMACS dump)
│   ├── peo20_solu/             # AA PEO-20 production: NAMD, 2000 frames, solvated
│   ├── polyply_PEO20/          # CG topology + AA→CG mapping artifact
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
