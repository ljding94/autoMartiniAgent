# examples

Runnable inputs for the agent repair loop (`agent/loop.py`, `--policy scripted`).
A scripted-ops file is a JSON list of actions `[{ "thought": ..., "ops": [...] }]`
replayed in order (then the loop submits). Ops name **heavy atoms only** — each
heavy atom's hydrogens move with it (via the `--cg-struct` reference frame).

## `psbma_w2_ops.json`

The hand-found **W2** repair for PSBMA-20: shift each sidechain's bead windows one
carbon (C9 propyl→ammonium, C6 ammonium→ester), re-centering the Q1 charge on N⁺.
Reproduces the loop objective **0.0417 → 0.0370** (stride 25; 0.0400 → 0.0352 at
full resolution).

Run from the repo root with the `autom3` env active:

```sh
PYTHONPATH=. python -m agent.loop \
  --mapping   reference/PSBMA_20mer/PSBMA20_mapping.json \
  --itp       reference/PSBMA_20mer/PSBMA20.itp \
  --aa-top    reference/PSBMA_20mer/PSBMA_20mer_no_water.gro \
  --aa-traj   reference/PSBMA_20mer/PSBMA_20mer_no_water_skip10.xtc \
  --cg-struct reference/PSBMA_20mer/PSBMA_20mer_no_water.gro \
  --policy scripted --scripted-ops examples/psbma_w2_ops.json \
  --frame-stride 25
```

Outputs land in `derived/PSBMA20/loop/` (override with `--work-root`):
`PSBMA20_loop_best_mapping.json`, `PSBMA20_loop_best.itp`, and `trajectory.json`.
