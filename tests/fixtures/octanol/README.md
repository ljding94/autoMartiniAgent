# octanol — golden fixture

- **SMILES**: `CCCCCCCCO` (1-octanol, 9 heavy atoms)
- **MOL name**: `OCOL` (avoids `OCT`/octane collision in `martini_v3.0.0_solvents_v1.itp`)
- **Command**: `python -m auto_martiniM3 --smi "CCCCCCCCO" --mol OCOL --canon -v --fpred`
- **Result**: success, 12 optimization iterations, ~1 s wall.
- **Mapping**: 3 beads — `SC1` (CCC), `SC1` (CCC), `SP2d` (CCO). All S-typed (small), 3 heavy atoms each. Sizing rules satisfied.
- **Bonds + angles**: present.
