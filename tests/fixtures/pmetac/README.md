# pmetac — failure case (M3 crash)

- **SMILES**: `CC(C)C(=O)OCC[N+](C)(C)` (PMETAC sidechain, 11 heavy atoms, +1 charge)
- **MOL name**: `PMTC`
- **Command**: `python -m auto_martiniM3 --smi "CC(C)C(=O)OCC[N+](C)(C)" --mol PMTC --canon -v --fpred`
- **Result**: **fails to produce output**. `gen_molecule_smi: Error. Only one molecule may be provided. C.C=O` — the optimizer's fragmentation produces disconnected SMILES intermediates that can't be re-parsed. The tool retries (3 occurrences in the 1274-line log) and exits without writing `.itp` / `.gro`.
- **Wall time before failure**: ~7 s.
- **Diagnosis (preliminary)**: bug in `auto_martiniM3.topology.gen_molecule_smi` rejecting disconnected intermediate fragments produced upstream by the optimizer. Charged species (the `[N+]`) seem to push the fragmentation into this state.

This is a **higher-priority failure mode** for our repair loop than the stall described in Chris's email — it's deterministic, reproducible, and bites a real polymer-relevant case.

## Note on email-chain divergence

Chris's original PMETAC test reportedly produced a *runnable but mis-sized* `.itp` in ~1 min. Bead names in his output (`C3, P1, C5, Qd`) are M2-style, indicating he ran the original `auto_martini` (not the M3 fork we're testing here). M3 fails earlier and harder.

## Files

- `auto_martiniM3.log` — full optimizer log (171 KB; see lines 326, 653, 1273 for the crash points)
- `run.stdout.log` — captured stdout/stderr
