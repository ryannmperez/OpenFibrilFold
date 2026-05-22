<p align="center">
  <img src="assets/logo.svg" alt="OpenFibrilFold logo" width="240" />
</p>

<h1 align="center">OpenFibrilFold</h1>

<p align="center">
  <em>An OpenFold3 fine-tune for amyloid fibril structure prediction.</em>
</p>

<p align="center">
  <a href="https://github.com/aqlaboratory/openfold-3"><img alt="upstream" src="https://img.shields.io/badge/built_on-OpenFold3-2b6cb0.svg"></a>
  <a href="LICENSE"><img alt="license" src="https://img.shields.io/badge/license-Apache--2.0-1f6feb.svg"></a>
  <img alt="status" src="https://img.shields.io/badge/status-research_preview-orange.svg">
</p>

---

## Overview

OpenFibrilFold adapts [OpenFold3](https://github.com/aqlaboratory/openfold-3) (an open reimplementation of AlphaFold 3) to the structural-prediction problem of **amyloid fibrils**: ribbon-like protein assemblies built from cross-β–stacked monomers, often in many coexisting *polymorphs*. Off-the-shelf folding models tend to either collapse onto a single fibril packing or produce geometrically-inconsistent ribbons because their cropping, distogram, and confidence machinery weren't designed for periodic ribbon assemblies. OpenFibrilFold introduces:

- a **ribbon-symmetric crop** that preserves both faces of the cross-β stack during training,
- a **cross-ribbon distogram weight** that emphasizes the fold-defining inter-strand contacts,
- **polymorph-aware training and evaluation**, including per-sample best-matching ground-truth selection and a sample-diversity loss against fibril polymorphs,
- a **fibril-specific dataset** of publicly deposited PDB structures, curated for split cleanliness and binding-site diversity.

We are also interested in the model's capability to **predict ligand poses bound to fibrils**, since many recent cryo-EM amyloid depositions include co-bound small molecules at potential drug-discovery sites. Ligand-specific lDDT metrics are reported alongside the protein metrics below.

## Highlights

> [!NOTE]
> Validation is a single full-sequence pass (`trainer.validate`) over the held-out fibril set, run end-to-end against unmodified OpenFold3 base weights using the *same* config (val pool, diffusion sampling, polymorph matching, MSAs, templates). The only difference between the two columns is the model weights. Two dataloaders, 13 PDBs total: `val_unique_seq` (n=6, rare-fold apo + ligand fibrils on unseen sequences) and `val_ligand` (n=7, ligand-bound fibrils). Numbers below are aggregates across both; per-dataloader and per-PDB breakdowns are in [`docs/val_compare.md`](docs/val_compare.md).

<p align="center">
  <img src="figures/val_ligand_tm.png" alt="Per-PDB fibril-assembly TM-score on val_ligand: base OF3 vs OpenFibrilFold (paired bars)" width="780" />
</p>

<sub>Per-PDB fibril-assembly TM-score on `val_ligand` (USalign multi-chain complex mode). OpenFibrilFold scores higher than base OF3 on every PDB in the set; per-PDB margins range from +0.15 to +0.32. The `val_unique_seq` breakdown is in [`docs/val_compare.md`](docs/val_compare.md).</sub>

<!-- METRICS_TABLE_START -->
| Metric | OpenFold3 base | OpenFibrilFold |
| --- | --- | --- |
| **Structure (lDDT / TM, ↑)** |  |  |
| Intra-protein lDDT | 0.477 | **0.614** |
| Inter-protein lDDT | 0.193 | **0.482** |
| Intra-complex lDDT | 0.479 | **0.615** |
| TM-score (USalign, complex) | 0.333 | **0.539** |
| **Ligand (lDDT, ↑)** |  |  |
| Intra-ligand lDDT | 0.876 | **0.893** |
| Intra-ligand lDDT (uha) | 0.729 | **0.754** |
| Inter-ligand lDDT | 0.232 | **0.313** |
| Protein–ligand lDDT | 0.069 | **0.086** |
| **Geometric (↓ lower better; GDT ↑)** |  |  |
| Distogram loss | 1.375 | **1.180** |
| Scaled distogram loss | 0.0412 | **0.0354** |
| Intra-protein dRMSD (Å) | 15.31 | **11.97** |
| Intra-ligand dRMSD (Å) | 0.800 | **0.698** |
| Complex RMSD (Å) | 30.53 | **19.22** |
| GDT-TS | 0.024 | **0.129** |
| GDT-HA | 0.004 | **0.057** |
| **Confidence (pLDDT magnitude, ↑)** |  |  |
| pLDDT (protein) | 0.254 | **0.440** |
| pLDDT (complex) | 0.253 | **0.434** |
| pLDDT (ligand) | 0.175 | **0.230** |

<sub>Both columns: same `trainer.validate` pass on the same held-out fibril set (13 PDBs), identical config — only model weights differ. Aggregates are weighted across the two dataloaders (n=6 + n=7). Ligand rows are reported from `val_ligand` only. TM-score is USalign multi-chain complex alignment on the val pipeline's pred + GT CIFs; it scores the full fibril assembly (monomer fold *and* stacking geometry) in one number. Pearson(lDDT, pLDDT) is intentionally not reported in this headline table — at the current pLDDT magnitudes (still in the 0.2–0.5 band) the per-PDB Pearson swings on tiny shifts and is not a reliable calibration signal yet. See [`docs/val_compare.md`](docs/val_compare.md) for the per-dataloader breakdown and the Pearson numbers with caveats.</sub>
<!-- METRICS_TABLE_END -->

## Head-to-head predictions

Three held-out fibrils, viewed end-on (looking straight down the protofilament stacking axis) so the cross-section that defines the fold is visible. **Top row:** OpenFold3 base. **Middle:** experimental cryo-EM ground truth. **Bottom:** OpenFibrilFold. Columns were picked by per-PDB ΔTM (OpenFibrilFold − base OF3) within each dataloader, so the figure spans the full range of behaviour on the held-out set: **9qlu** — largest ΔTM on `val_unique_seq` (+0.253); **9ug1** — largest ΔTM on `val_ligand` (+0.324); **9ljb** — smallest ΔTM in either set (+0.066), a 1010-residue assembly where both models still struggle. Each column is rendered at the same physical scale; the prediction panels are coloured by per-residue lDDT against the experimental reference — saturated colour = aligned, faded toward white = mis-aligned. The number under each prediction is the per-PDB intra-complex lDDT computed with the *same* function the val pipeline uses (`openfold3.core.metrics.validation_all_atom.lddt`) on heavy atoms after Hungarian chain-permutation matching, so values are on the same scale as the metric table above.

<!-- HEAD_TO_HEAD_START -->
<p align="center">
  <img src="figures/h2h_grid.png" alt="3×3 head-to-head: OF3 base vs ground truth vs OpenFibrilFold on 9qlu (largest val_unique_seq ΔTM), 9ug1 (largest val_ligand ΔTM), 9ljb (smallest ΔTM in either val set)" width="780" />
</p>
<!-- HEAD_TO_HEAD_END -->

## Methods (in brief)

- **Backbone:** OpenFold3 trunk (Pairformer + diffusion module + confidence heads), unmodified.
- **Training set:** publicly deposited cryo-EM amyloid-fibril structures from the PDB, filtered by resolution (≤ 9 Å) and minimum chain count for ribbon-axis assignment, then deduplicated by entity sequence and ligand binding-site clustering. Train / validation split is held fixed across all experiments.
- **Polymorph determination:** structures are grouped by protein identity (UniProt accession), then clustered within each group by mean-absolute-difference of CA–CA distance matrices (≈ 4 Å threshold) — distinct fibril packings of the same protein become distinct polymorphs. Ligand-bound structures additionally contribute one binding-site polymorph per resolved site.
- **Crop:** ribbon-symmetric post-processor — ratio-capped at 1.3, token budget 484.
- **Distogram:** explicit cross-ribbon weight on β-stacking pairs.
- **Diffusion loss:** inter-chain weights tuned for stacked, in-register β-strands (`intra=1`, `adjacent=100`, `other=5`).
- **Polymorph awareness:** per-sample best-matching polymorph as ground truth (training and val), plus a softmax-entropy diversity loss across diffusion samples.
- **Stages:** AF3-style two-stage fine-tune (structural trunk first, confidence-only second).
- **Hardware:** 3× NVIDIA RTX A6000 Ada (DDP), bf16-mixed precision, 32 diffusion samples per structure at train and val. End-to-end training in ~1 day.

For the full description, including the rationale for each loss-module change, see [`docs/methods.md`](docs/methods.md).

## Upcoming

The following will be released alongside the public flip of this repository:

- **Curated fibril dataset** — preprocessed structures, splits, MSAs, ribbon assignments, and polymorph groups.
- **Trained model weights** — final OpenFibrilFold checkpoint and the inference config to use it.
- **Training scripts** — full reproducible training pipeline (data prep → stage 1 → stage 2 → eval), including the OpenFold3 patch set.
- **Manuscript / preprint** — describing the architectural changes, training recipe, and benchmarking results in detail.

## How it was built

OpenFibrilFold was developed entirely as a **"vibe-coded"** project using [Claude Code](https://www.anthropic.com/claude-code) (Opus 4.6 and 4.7) with agent teams parallelizing data curation, patch development, debugging, ablation analysis, and figure generation. The OpenFold3 patch set, training infrastructure, dataset pipeline, polymorph machinery, and this repository were all written in collaboration with the model. End-to-end training fit on **3× NVIDIA RTX A6000 Ada** over the course of a single day.

## Acknowledgements

OpenFibrilFold is a fine-tune of [OpenFold3](https://github.com/aqlaboratory/openfold-3) (AQ Laboratory). Both OpenFold3 and the reference [AlphaFold 3](https://www.nature.com/articles/s41586-024-07487-w) algorithm are foundational to this work. Cross-β fibril structures used for training and validation come from depositions in the [Protein Data Bank](https://www.rcsb.org/).

## Citation

If this work is useful to you, please cite the repository:

```bibtex
@software{perez_openfibrilfold_2026,
  author       = {Perez, Ryann},
  title        = {OpenFibrilFold: an OpenFold3 fine-tune for amyloid fibril structure prediction},
  year         = {2026},
  publisher    = {GitHub},
  url          = {https://github.com/ryannmperez/OpenFibrilFold}
}
```

## License

[Apache 2.0](LICENSE).

## Release notes

- **2026-05-07** — initial private snapshot. Renamed from FibrilFold → OpenFibrilFold.
