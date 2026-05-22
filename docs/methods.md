# Methods

OpenFibrilFold (OFF) is a fine-tune of [OpenFold3](https://github.com/aqlaboratory/openfold-3) (an open-source reimplementation of AlphaFold 3) targeting amyloid fibril structure prediction. The model retains OF3's full architecture and is adapted through (1) targeted loss-module changes, (2) a ribbon-symmetric cropping strategy, (3) polymorph-aware training and evaluation, and (4) a fibril-specific dataset.

## Data

- **Source.** Publicly deposited amyloid-fibril structures from the PDB.
- **Filtering.** Resolution ≤ 9 Å (cryo-EM), minimum chain count to ensure ribbon-axis assignment, deduplicated by entity sequence and binding-site clustering.
- **Splits.** Train / validation split is held fixed across all experiments; validation PDB identifiers are not released to keep the holdout clean for future benchmarks.
- **MSAs and templates.** OpenFold3-style alignments computed against the standard sequence databases. Templates include both PDB hits and a synthetic "ribbon template" that injects a stacked-β prior into structures lacking close PDB precedents.

## Architecture changes

OFF preserves the OpenFold3 trunk (Pairformer + diffusion module + confidence heads) verbatim — no new modules. The fine-tuning logic is concentrated in the data, loss, and config layers.

### Ribbon-symmetric cropping
- Standard OF3 spatial / interface crops can clip a fibril asymmetrically — e.g., 6 chains on one face and 1 on the other — yielding training samples whose β-sheet stacking is geometrically inconsistent.
- Custom post-processor symmetrizes crops about the ribbon axis: after the spatial crop is selected, chains on the under-represented side are re-included up to a ratio cap (`symmetric_max_ratio = 1.3`).
- Token budget held at 484; symmetric cap holds the worst-case crop ≈ 629 tokens (~25% memory headroom vs an uncapped symmetrize).

### Loss-module adaptations
- **Cross-ribbon distogram weighting.** Adds an explicit weight for distogram pairs that span ribbon halves, emphasizing the fold-defining cross-β contacts that vanilla distogram loss treats as "just another long-range pair."
- **Inter-chain weighting.** Diffusion loss inter-chain weights set to `intra=1.0`, `adjacent=100.0`, `other=5.0` — the adjacent-chain term dominates so the gradient is shaped by the geometry of stacked, in-register β-strands.
- **Bond loss disabled.** The AF3 bond term is built from `token_bonds × (is_polymer × is_ligand)` — i.e. covalent polymer–ligand linkages. Amyloid + small-molecule binders in this dataset are noncovalent, so the term is identically zero on every batch and is left at weight 0 to avoid logging confusion.
- **Polymorph diversity loss.** Diffusion samples are soft-assigned to the per-structure polymorph set via `softmax(lDDT / T)`; the entropy of the mean assignment is pushed toward `log(J)` where `J` is the number of polymorphs. Penalizes mode-collapse onto a single conformer.

### Polymorph-aware training and evaluation
- Each PDB carries a precomputed polymorph group (different fibril packings of the same protein). During training, every sample is paired with the per-sample best-matching polymorph as ground truth; the matching is run on each diffusion sample and on the lDDT scoring.
- During validation, the same matching is applied (`polymorph_match_val_gt = true`). Without it, a sample matching polymorph 3 would be scored against polymorph 0 and produce a spurious negative pLDDT/lDDT correlation.
- Ligand binding-site polymorphs: each distinct ligand binding site contributes a polymorph variant where non-site ligand atoms are masked out of the loss, teaching the confidence head "if you commit to this site, the others should disappear."

## Training

- **Stages.** AF3-style staged fine-tune. Stage 1 (structural) trains the trunk on `mse + smooth_lddt + distogram` with confidence heads disabled; Stage 2 freezes the trunk and trains only the confidence heads (pLDDT / PAE / PDE / experimentally_resolved) under AF3's published Stage-2 weights.
- **Hardware.** 3× NVIDIA RTX A6000 Ada (DDP), bf16-mixed precision, accumulation factor 8. Effective batch size 24 structures per optimizer step. End-to-end training (data prep, both stages, eval) completes in ~1 day on this machine.
- **Schedule.** Linear warmup 63 steps → peak LR 1e-3 → step decay every 125 steps with factor 0.9 starting at step 375. EMA decay 0.999.
- **Diffusion.** 32 diffusion samples per structure during training and validation.

## Evaluation

- All numbers reported in this repository are from a single full-sequence validation pass (`ValidationPDBDataset`, no cropping at val time) on the held-out fibril set.
- Comparison to OpenFold3 base is run end-to-end using the same configuration (same val pool, same diffusion sampling, same polymorph matching) — the only difference is the model weights.
- Metrics reported: lDDT (intra/inter, protein/ligand/complex), distogram loss, GDT-HA / GDT-TS, dRMSD, RMSD, Pearson(lDDT, pLDDT), and clash rates.

## What is not released

- Validation PDB identifiers and their MSAs/templates, to preserve the held-out set for future benchmarks.
- The raw curated training dataset and dataset-construction scripts.
- Training run logs containing validation-set fingerprints.

The model weights, model and loss configuration, and model-architecture diff against upstream OpenFold3 are released.
