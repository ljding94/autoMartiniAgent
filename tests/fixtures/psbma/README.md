# psbma — surprise success (needs chemistry review)

- **SMILES**: `CC(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]` (PSBMA monomer, zwitterionic, 18 heavy atoms)
- **MOL name**: `PSBA`
- **Command**: `python -m auto_martiniM3 --smi "..." --mol PSBA --canon -v --fpred`
- **Result**: **success in 5.4 s wall, 1 optimization iteration**. Contradicts Chris's email-chain prediction of a 30+ min stall on the original `auto_martini`.
- **Mapping**: 6 beads, sizing rules satisfied:

  | bead | type | size | atoms          | heavy | smiles fragment | charge |
  |------|------|------|----------------|-------|-----------------|--------|
  | C01  | SC1  | S    | C0, C1, C2     | 3     | CCC             |  0     |
  | N01  | TN5a | T    | C3, O4         | 2     | C=O             |  0     |
  | P01  | TP2a | T    | O5, C6         | 2     | CO              |  0     |
  | 101  | Q1   | R    | C7, N8, C9, C10| 4     | C[N+](C)C       | +1     |
  | C02  | SC1  | S    | C11, C12, C13  | 3     | CCC             |  0     |
  | 102  | Q1   | R    | S14, O15, O16, O17 | 4 | O=[SH](=O)[O-]  | −1     |

- **Bonds + angles**: present, 5 bonds, 10 angles.
- **Mapping additivity ratio**: 3.65 (whole vs sum: −0.27 vs 0.72) — flagged by the tool itself as a warning condition; meaning the assumption that bead free energies are additive doesn't hold well for this zwitterion. Worth investigating.

## Open question for collaborators

This is the canonical Phase 0 ask: **does this 6-bead partition match Chris and Seonghan's chemistry intuition for PSBMA?** Specifically:

- Is splitting the ester into `TN5a (C=O)` + `TP2a (CO)` the right move, or should they merge into a single S-bead?
- Q1 with explicit ±1 charges for the ammonium and sulfonate — correct in principle, but is `Q1` the right subtype for both, or should the sulfonate be `Q5n` / similar?
- The ammonium-sulfonate distance (4 bonds) — does the zwitterionic interaction need any special bonded treatment, or do the standard non-bonded LJ + Coulomb suffice?

## Files

- `PSBA.itp`, `PSBA.gro` — generated outputs
- `auto_martiniM3.log` — full optimizer log (530 KB)
- `run.stdout.log` — captured stdout
