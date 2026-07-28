# All-atom reference structures — PP, PLA, PS, PVP, PMMA

Water-free single-chain structures for the five polymers marked "newly added" in
[`docs/TEST_SET.md`](../docs/TEST_SET.md). Each is the starting point for AA → Martini 3
mapping and for extracting target bonded distributions.

| File | Atoms | Monomers | Residue names | Equilibrated box (nm) |
|---|---|---|---|---|
| `pp/pp_no_water.gro` | 182 | 20 | `PROR` / `PROS` | 7.197 |
| `pla/pla_no_water.gro` | 182 | 20 | `LACTR` / `LACTS` | 7.019 |
| `ps/ps_no_water.gro` | 322 | 20 | `STYRR` / `STYRS` | 6.858 |
| `pvp/pvp_no_water.gro` | 342 | 20 | `VIPR` / `VIPS` | 6.832 |
| `pmma/pmma_no_water.gro` | 302 | 20 | `MMAR` / `MMAS` | 6.955 |

All five are N = 20 as specified in `TEST_SET.md`. Atom counts equal
20 × (monomer formula) + 2 terminal hydrogens — e.g. PVP 20 × C6H9NO (17) + 2 = 342.
The `R`/`S` residue-name pairs are the two stereo-configurations of the same monomer;
the chains are atactic, so both appear along one backbone.

## Provenance

Source: CHARMM-GUI Polymer Builder → NAMD, single chain solvated in a cubic water box.
Simulations live on external drive T9 at
`backup_research_onrl/systems/<sys>/namd/` (not in git):

- production `step5_1.dcd` … `step5_200.dcd` = 200 × 1 ns = **200 ns**
- timestep 2 fs, `dcdfreq` 50000 → 1 frame / 100 ps → 2000 frames total
- `wrapAll on`; cubic box; CHARMM36 (`toppar_all36_synthetic_polymer`)

Each `.gro` here is the **first production frame** (t = 100 ps). Water (`SOLV`
segment, TIP3) was stripped; the systems contain no ions (`npos = nneg = 0`), so
what remains is only the polymer chain (segment `S1P1`). The chain was made whole
across periodic boundaries via a bond-graph unwrap using the `!NBOND` list from
`step3_input.psf`, then shifted so its centre of geometry sits inside the primary
cell. Bond lengths were verified to stay within 1.1–1.7 Å in every frame.

## Trajectories

The full 200 ns water-free trajectories are committed alongside each structure as
multi-frame `.gro`, since the single frames above are not sufficient for fitting
bonded distributions:

| File | Frames | Spacing | Size |
|---|---|---|---|
| `pp/pp_no_water_1-200ns.gro` | 2000 | 100 ps | 16 MB |
| `pla/pla_no_water_1-200ns.gro` | 2000 | 100 ps | 16 MB |
| `ps/ps_no_water_1-200ns.gro` | 2000 | 100 ps | 28 MB |
| `pvp/pvp_no_water_1-200ns.gro` | 2000 | 100 ps | 29 MB |
| `pmma/pmma_no_water_1-200ns.gro` | 2000 | 100 ps | 26 MB |

Frame *i* is at *i* x 100 ps, so the span is 100 ps -> 200000 ps. Coordinate
precision is 0.001 nm.

> **Note on format choice.** These total ~115 MB; the identical data as `.xtc` is
> ~15 MB, an 8x saving at the same 0.001 nm precision, and the repo's `.gitignore`
> otherwise keeps bulk MD trajectories out of git (see its "Bulk MD trajectories"
> section and `PSBMA_20mer_no_water_skip10.xtc`). The `.gro` form was committed by
> explicit request. If repo size becomes a problem the `.xtc` equivalents are on T9
> at `backup_research_onrl/systems/combined_nowat/` and are drop-in replacements --
> `gmx` and MDAnalysis read either.

The converter (`dcd2gmx.py`) and its own README are on T9 in that same directory.

## Notes for the mapping pipeline

- **No velocities.** DCD stores coordinates only, so these files have 3 columns per
  atom, unlike `PSBMA_20mer/PSBMA_20mer_no_water.gro` which carries velocities.
- **Real residues, not a single `LIG`.** Residues are numbered 1…20 with the
  chemistry-specific names above, whereas `PSBMA_20mer_no_water.gro` puts every atom
  in one residue `LIG` resid 1. Code that assumes a single-residue ligand may need a
  flattening pass.
- **Box size.** `TEST_SET.md` describes the setup as a cubic 8 nm box; the actual
  CHARMM-GUI cells were ~7.0–7.2 nm before NPT and 6.83–7.20 nm after. Worth
  reconciling in that doc. The stated rationale (extended N = 20 contour ≈ 5–6 nm plus
  ~1.2 nm cutoff on each side) is close to marginal at 6.8 nm for the longest chains.
- **Per-frame re-centring** means absolute positions are not continuous in time, so
  the trajectories on T9 cannot be used for centre-of-mass diffusion. Internal
  geometry — bonds, angles, dihedrals, Rg, end-to-end — is unaffected, which is what
  bonded-parameter fitting needs.
