# Polyply PEO-20 — CG topology + AA→CG mapping

CG cold-start artifacts for the PEO-20 chain whose AA reference lives at
`reference/gromacs/` (CHARMM36, TIP3P, 200 ns, single chain).

## Files

| file | source | role |
|---|---|---|
| `PEO20.itp` | `polyply gen_params -lib martini3 -seq PEO:20 -o PEO20.itp -name PEO20` | CG topology — 20 SN3r beads, central-bond r₀=0.36 nm k=7000, central-angle θ₀=123° k=80 |
| `PEO20_raw.gro` | `polyply gen_coords` random walk | initial CG coords pre-EM |
| `PEO20_clean.gro` | `gmx trjconv -pbc mol -center` after `em.tpr` minimisation | post-EM CG coords |
| `topol.top`, `em.*` | GROMACS EM artifacts | provenance only |
| `PEO20_mapping.json` | `scripts/build_peo20_mapping.py` | **canonical AA→CG mapping** — atom-index keyed, mass-weighted, loaded by the projector |
| `PEO20.map` | same script | Martini-style human-readable mirror of the same mapping |

## Mapping rule

One SN3r bead per ETHOX residue: bead *i* ← residue *i*. All 7 atoms of a
mid-chain ETHOX (`O1, C1, H11, H13, C2, H21, H23`) project into one bead by
mass-weighted COM. End groups fold in:

- **Bead 1 (HO terminus)** — residue 1 contributes 8 atoms: standard 7 +
  the hydroxyl proton `HO1`. AA mass 45.061 g/mol.
- **Bead 20 (CH₃ terminus)** — residue 20 contributes 8 atoms: standard 7
  + the third methyl-H `H22` (since `C2` is `CG331` not `CG321`). AA mass
  45.061 g/mol.

Mid-chain beads (2–19): 7 atoms, AA mass 44.053 g/mol. The Polyply `.itp`
uses 45.0 for every bead — the ~1 g/mol overshoot on mid-chain beads is
the standard Martini PEO rounding and is not a mapping error.

## Sanity check on `equil.gro` (single frame)

| metric | value | Polyply equilibrium | note |
|---|---|---|---|
| 1-2 bond, mean | 0.310 nm | 0.360 nm | single equilibrated frame, distribution expected to broaden over the 200 ns trajectory |
| 1-2 bond, range | 0.250 – 0.365 nm | — | |
| 1-2-3 angle, mean | 118.2° | 123° | |
| 1-2-3 angle, range | 100.3° – 155.0° | — | |
| bead mass sums | match expected AA sums to 0.01 g/mol | — | rules out atom-index errors |

The 0.05 nm offset on the mean 1-2 bond is the kind of signal the scorer
(#3 in `PROGRESS.md`) is built to detect — CHARMM PEO with this mapping
may want a tighter equilibrium bond than the Polyply default. Whether
that's a mapping problem or a parameter-fitting one is the scorer's job to
adjudicate.

## Re-generating

```sh
conda activate autom3
python scripts/build_peo20_mapping.py    # rebuilds .json + .map from the two .itp files
python scripts/check_peo20_mapping.py    # asserts masses + prints 1-2 and 1-2-3 stats
```

## AA system topology — confirmed structure

From `reference/gromacs/toppar/S1P1.itp` (142 atoms, 20 ETHOX residues):

- Residue 1: `O1(OG311)-HO1-C1-2H-C2-2H` — hydroxyl terminus
- Residues 2–19: `O1(OG301)-C1-2H-C2-2H` — ether monomer
- Residue 20: `O1(OG301)-C1-2H-C2(CG331)-3H` — methyl-terminated

Net chain: **HO-CH₂-CH₂-O-(CH₂-CH₂-O)₁₈-CH₂-CH₃**, i.e. an ethyl-ether
terminated 19-mer of ethylene oxide. The `PROGRESS.md` note that called
the terminus `-CH₂-CH₂-CH₃` was off by one CH₂ — the topology has only one
methylene between the last ether O and the methyl C.
