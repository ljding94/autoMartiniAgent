# Phase 0 fixtures — AutoMARTINI3 (M3) baseline behavior

**Backend under test**: [`Automartini_M3`](https://github.com/Martini-Force-Field-Initiative/Automartini_M3) at commit `1fff05a` (2026-03-23, "repare VS bonded beads bug").
**Environment**: `autom3` conda env, Python 3.10, RDKit 2026.3.1, AutoMARTINI3 installed editable from source.
**Date captured**: 2026-05-01.

## Reproduction

```bash
conda activate autom3
cd <fixture-dir>
python -m auto_martiniM3 --smi "<SMILES>" --mol <NAME> --canon -v --fpred
```

The exact SMILES per fixture is in each subdirectory's `README.md`.

## Findings vs. expectations from the email chain

The April 2026 email chain (Walker, Kim, Ding) made three predictions about AutoMARTINI3 behavior. Reproducing on the M3 fork shows two of them **do not hold**:

| molecule | predicted (from email)                                  | observed (M3 fork, 2026-05-01)                                       |
|----------|---------------------------------------------------------|----------------------------------------------------------------------|
| octanol  | works, EM converges                                     | **confirmed**. Converges in 12 iterations, ~1 s. Sizing rules satisfied (3× S-beads, 3 heavy atoms each). |
| PMETAC   | works in ~1 min, **mis-sizes beads** (2–3 heavy in R)  | **does not match**. M3 fork crashes with `gen_molecule_smi: Error. Only one molecule may be provided. C.C=O` — the optimizer produces disconnected SMILES intermediates. No `.itp` is generated. |
| PSBMA    | **stalls > 30 min**                                     | **does not match**. M3 fork converges in 5.4 s, produces a Martini-3-rule-compliant 6-bead mapping. |

Almost certainly because **Chris was running the original `auto_martini`** (Bereau & Kremer's M2-era tool), not the M3 fork:

- His PMETAC output uses M2-style bead names (`C3, P1, C5, Qd`), not M3's `SC1, TN5a, SP2d, Q1` style.
- The M3 fork has a restructured optimization that's faster and (anecdotally, on PSBMA) sizing-compliant.

## Implications for the project plan

1. The "QA + repair loop" central premise still holds, but the specific failure modes we'll repair are different from what the email chain documented. We need a broader test panel to characterize where M3 actually breaks.
2. **The PSBMA-stall-recovery feature (subprocess timeout + fragment-and-assemble fallback) is no longer on the critical path.** It may still be needed for some molecule classes, but PSBMA itself doesn't need it under M3.
3. **A new urgent failure mode**: the disconnected-fragment crash that bites PMETAC. This is what the agent's repair loop must actually address first.
4. We should validate the M3-generated PSBMA mapping with Chris and Seonghan — does the `SC1 / TN5a / TP2a / Q1(+) / SC1 / Q1(-)` partition match their chemistry intuition?
5. We should also test on the original `auto_martini` if we want to recover Chris's exact failure cases for parity (deferred).

## Per-fixture details

- `octanol/` — golden, M3 succeeds with rule-compliant output.
- `pmetac/` — failure case: M3 crash on disconnected-fragment intermediates.
- `psbma/` — surprise success: M3 produces a clean 6-bead mapping in 5 s; needs human review for chemistry correctness.
