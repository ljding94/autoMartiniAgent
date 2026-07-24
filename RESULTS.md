# Results notebook — evidence for the paper

A running log of paper-worthy findings as we build. This is the source of truth
for [`result_note.html`](result_note.html) (generated via
`python scripts/build_result_note.py`). Add an entry whenever a result lands.

Each entry: the **claim** (one line, paper-ready), the **figure(s)**, the key
**numbers**, the **finding**, and **for the paper** (where/how it's used).

**Scope so far**: PSBMA-20 (zwitterionic methacrylate; 6 beads/monomer; 25 001 AA
frames). Everything is **simulation-free** — we project the *existing* AA trajectory
through the mapping; no AA or CG simulation is run.

## R1 · The agent autonomously repairs a CG mapping
- **Claim**: Given only a simulation-free Gaussianity objective and a curated action space, an LLM agent autonomously discovers a chemically-motivated mapping repair.
- **Figure**: derived/PSBMA20/repair/bond_hists_compare.pdf
- **Figure**: derived/PSBMA20/repair/angle_hists_compare.pdf
- **Numbers**: objective `0.0400 → 0.0352` (−12 %, full 25 001 frames). `Q1-TP2a` bond R² 0.892→0.996; the bimodal `SC1-Q1-TP2a` angle 0.840→0.894. Every bead stays ≤4 heavy (rule-valid).
- **Finding**: opus-4.8 (via OpenRouter), reasoning explicitly about re-centering the ammonium charge on N⁺, proposes the coupled edit `C6→ester + C9→ammonium` — the same repair we had found by hand. Found reliably with best-of-3.
- **For the paper**: the headline "the method works" result. The model's reasoning trace (re-centring the charge) is worth quoting verbatim.

## R2 · The objective orders mappings by quality (rank)
- **Claim**: AA-projected bonded-distribution Gaussianity is a valid mapping-quality signal — it degrades monotonically as a good mapping is chemically degraded.
- **Figure**: derived/PSBMA20/rank_validation.pdf
- **Numbers**: good mapping `0.0417`; random-scramble severity ×{1,2,4,8} means `0.058 / 0.067 / 0.118 / 0.194` (monotonic); every functional-group split scores worse (split_ammonium +0.067, split_sulfonate +0.048); 34/36 perturbations worse than the good mapping.
- **Finding**: perturbing the good mapping (named functional-group splits + a random-scramble ladder, fixed bead count) monotonically raises the objective.
- **For the paper**: validation figure 1 — answers the reviewer's core objection, "is the objective a valid target?", using the good mapping as a relative ground truth (no simulation).

## R3 · The agent recovers damaged mappings (recover)
- **Claim**: The objective does not just *rank* mappings — it *guides* the agent back to good ones from deliberate damage.
- **Figure**: derived/PSBMA20/recover_validation.pdf
- **Numbers**: mean **107 %** of the degradation gap closed; every functional-group split recovered; 3/4 over-recover to the W2 optimum; grouping match 77–85 %.
- **Finding**: degrade the good mapping (functional-group split), run the loop, and it climbs back — often past the good baseline to the W2 optimum.
- **For the paper**: validation figure 2 — "the objective's landscape leads to good mappings." Pairs with R2 to complete the sim-free validation.

## R4 · Chemical reasoning keeps the search honest (ablation)
- **Claim**: The Gaussianity objective is *exploitable*; the LLM's chemical reasoning supplies validity that a blind optimizer needs explicit rules to match.
- **Figure**: derived/PSBMA20/ablation.pdf
- **Numbers**: greedy (objective-only) `0.0217` but **INVALID** (splits the sulfonate); greedy (chem-constrained) `0.0417`, **stuck**; LLM (no FG rule) `0.0370`, **valid**, **8 evals** (vs greedy's 79).
- **Finding**: steepest-ascent single-move greedy games the objective by splitting the sulfonate (S in one bead, its three O's in another) — chemically nonsensical. Constrained by functional-group integrity it can't improve at all. The LLM finds the valid W2 repair *without being given the FG rule*, and far more evaluation-efficiently.
- **For the paper**: the "why an LLM, not a hill-climber" figure — three wins at once (validity for free, capability, efficiency) plus a cautionary result that pure objective optimization is exploitable (this motivated the rule checker).

## Methods notes (rigor — for the Methods section)
- **PBC projection fix**: computing each bead's COM on raw coordinates collapses beads whose atoms straddle the periodic boundary → spurious ~2.5 nm "bonds" and a fake ~10° angle peak (hit 52.5 % of PSBMA frames). Fixed with a minimum-image unwrap (`minimize_vectors`) before the mass-weighted average. Essential caveat for anyone projecting AA→CG.
- **Goodness-of-fit = R² in density space**: RMSE is not comparable across terms (scales with peak density); KL / JS / cross-entropy *invert* vs. visual fit quality (dominated by ε-clipped tail bins). R² weights every bin's residual equally, so it tracks the eye and the physics.
- **Simulation-free**: the method never runs AA or CG simulations — it only projects existing AA trajectories through the candidate mapping. This is both a scope boundary and a selling point.
