# CG Mapping Test Set

Benchmark polymers for the autoMartiniAgent AA → Martini 3 mapping pipeline.
Each system is a deliberately chosen stress point, not an arbitrary example: the
set is laid out so that difficulty and available scaffolding fall off together,
turning the table into a graded stress ladder for the agent's QA/repair loop.

## Two axes

**Test role**

- **Validation** — a Martini 3 reference model exists, so the pipeline's output
  can be scored against an accepted target.
- **Novel** — no reference exists; this tests *de novo* mapping generation.

**Difficulty** — defined mechanically by the structural feature that stresses the
automatic pipeline, *not* by subjective judgment:

- **Easy** — linear backbone with small apolar side groups (bead typing trivial;
  no automatic step fails).
- **Medium** — heteroatom-bearing side chain or backbone stereocenter, no ring or
  charge (stresses bead-type assignment).
- **Hard** — ring (stresses rigid-body construction and symmetry handling).
- **Extreme** — formal charge/zwitterion with no reference model and no builder
  support, built from in-house all-atom models — the regime that originally broke
  the automatic mapper.

Difficulty increases monotonically with the loss of external scaffolding
(Martini 3 target, then AA builder), so the table doubles as a stress ladder.

## The set

| Polymer | Difficulty | Role | Martini 3 model | Rationale — what it adds to the set |
|---|---|---|---|---|
| **PEO**† <br>poly(ethylene oxide) | Easy | Validation | Available | Positive control on the simplest ether backbone; confirms the pipeline reproduces a well-validated, accepted mapping. |
| **PP** <br>polypropylene | Easy | Novel | Unavailable | Novel-side floor: polyethylene backbone plus one methyl; lowest-complexity *de novo* mapping with no reference to lean on. |
| **PMMA** <br>poly(methyl methacrylate) | Medium | Validation | Available | Validation at intermediate complexity — pendant ester, no ring or charge; stresses bead-type assignment against an existing target. |
| **PLA** <br>poly(lactic acid) | Medium | Novel | Unavailable | Only polyester backbone in the set: ester linkages and a chiral center alter backbone torsions and flexibility, probing bonded-parameter derivation beyond vinyl chains; high biomedical relevance. |
| **PS**\* <br>polystyrene | Hard | Validation | Available (provisional) | Aromatic ring with a reference model still under revision; tests whether the QA loop flags a questionable published target — the flagship "catch the error" case. |
| **PVP** <br>poly(vinylpyrrolidone) | Hard | Novel | Unavailable | Polar lactam ring (tertiary amide, H-bond acceptor only): a ring failure mode of different polarity than PS, with no reference. |
| **PSBMA**†‡ <br>poly(sulfobetaine methacrylate) | Extreme | Novel | Unavailable | Zwitterion — two opposite charges per monomer; in-house antifouling system that originally broke the automatic mapper (subprocess stalls, Martini-rule-violating bead sizes). |
| **PNOMA**†‡ | Extreme | Novel | Unavailable | Fully in-house antifouling monomer with no external scaffolding at any stage; the true end-to-end stress test on the target application domain. |

† Already part of the current dataset; the remaining five are newly added.
‡ In-house all-atom model built by Seonghan Kim (not available in Polymer Builder).
\* Martini 3 polystyrene parameters are provisional — bonded terms were carried
over from the Martini 2 model and are under active revision
([polyply discussion #379](https://github.com/marrink-lab/polyply_1.0/discussions/379)).

## Balance check

- **Difficulty**: Easy 2 / Medium 2 / Hard 2 / Extreme 2 — fully balanced.
- **Role**: Validation 3 (PEO, PMMA, PS) / Novel 5.
- **Backbone chemistry**: ether (PEO), ester (PLA), carbon-vinyl (the rest) — 3 types.
- **Rings**: PS (apolar aromatic) vs PVP (polar lactam) — two ring failure modes of
  different polarity.
- **Charge**: PSBMA, PNOMA — the target application domain fills this tier.

## All-atom reference generation

All-atom reference structures were generated with
[CHARMM-GUI Polymer Builder](https://www.charmm-gui.org/input/polymer) for the six
commercially standard polymers (PEO, PP, PMMA, PLA, PS, PVP) at a degree of
polymerization of N = 20, each as a single chain solvated in a cubic 8 nm box of
water; PSBMA and PNOMA, which are not available in Polymer Builder, were built and
simulated in-house by Seonghan Kim at the same chain length.

> **Box size rationale.** The 8 nm cube keeps a fully extended N = 20 chain (max
> contour length ≈ 5–6 nm) from contacting its own periodic image once the
> non-bonded cutoff (~1.2 nm) is added on each side. Use a single box size across
> all systems so box dimension is not a variable in the comparison. This applies to
> the single-chain-in-solution setup used for extracting bonded distributions; a
> melt setup would be sized by target density instead.
