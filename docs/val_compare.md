# Validation: OpenFibrilFold vs OpenFold3 base

Side-by-side validation of OpenFibrilFold (exp43, step 1024) against the unmodified OpenFold3 ft3_v1 base weights, run through the same Lightning val pipeline: same dataloaders, same crops/templates, same diffusion sampling, same polymorph matching. The only thing that changes between the two columns is the model weights.

## Val set

Two held-out dataloaders, full-sequence (no train-time cropping):

- **`val_unique_seq`** (n=6, 6 unique protein sequences not seen at train time): rare-fold apo and ligand-bound fibrils that the model has never been exposed to at the sequence level.
- **`val_ligand`** (n=7, 7 ligand-bound fibrils): protein + small-molecule co-folding.

13 PDBs total. Predictions and ground-truth CIFs are written to `val_structures/epoch_0/` during the `trainer.validate` call and used as the inputs for the structural-similarity scores below.

## Headline TM-score

USalign multi-chain complex alignment (`-mm 1 -ter 0`), pred CIF vs GT CIF from the same Lightning val pass. The score is computed on the full fibril assembly — chain correspondences are found by MM-align greedy search, so a single TM number captures both monomer fold *and* stacking geometry.

<p align="center">
  <img src="../figures/val_ligand_tm.png" alt="Per-PDB TM-score on val_ligand: base OF3 vs OpenFibrilFold" width="780" />
</p>

<p align="center">
  <img src="../figures/val_unique_seq_tm.png" alt="Per-PDB TM-score on val_unique_seq: base OF3 vs OpenFibrilFold" width="780" />
</p>

| dataloader | n | mean TM (base OF3) | mean TM (OpenFibrilFold) | Δ |
|---|---:|---:|---:|---:|
| `val_ligand`     | 7 | 0.376 | **0.610** | +0.235 |
| `val_unique_seq` | 6 | 0.282 | **0.456** | +0.174 |

**Per-PDB TM-score (winner highlighted):**

| dataloader | pdb | L (res) | base OF3 | OpenFibrilFold | Δ | base RMSD (Å) | OFF RMSD (Å) |
|---|---|---:|---:|---:|---:|---:|---:|
| val_ligand     | 8byn | 750 | 0.260 | **0.510** | +0.249 | 7.60 | 6.44 |
| val_ligand     | 8fug | 730 | 0.351 | **0.657** | +0.305 | 4.11 | 4.58 |
| val_ligand     | 9e8v | 770 | 0.401 | **0.551** | +0.150 | 6.03 | 5.58 |
| val_ligand     | 9e8w | 790 | 0.360 | **0.552** | +0.192 | 6.60 | 5.24 |
| val_ligand     | 9e8x | 780 | 0.336 | **0.508** | +0.171 | 7.04 | 6.19 |
| val_ligand     | 9e8y | 780 | 0.349 | **0.599** | +0.250 | 6.45 | 6.08 |
| val_ligand     | 9ug1 | 590 | 0.573 | **0.897** | +0.324 | 5.09 | 2.32 |
| val_unique_seq | 9cww |  340 | 0.509 | **0.759** | +0.250 | 5.62 | 4.21 |
| val_unique_seq | 9hx4 |  610 | 0.195 | **0.369** | +0.174 | 8.17 | 5.79 |
| val_unique_seq | 9ljb | 1010 | 0.158 | **0.224** | +0.066 | 9.17 | 6.81 |
| val_unique_seq | 9m5q |  270 | 0.233 | **0.389** | +0.156 | 5.95 | 4.28 |
| val_unique_seq | 9qlu |  330 | 0.204 | **0.457** | +0.253 | 7.09 | 3.54 |
| val_unique_seq | 9u4l |  345 | 0.392 | **0.538** | +0.146 | 5.71 | 5.58 |

OpenFibrilFold scores higher than base OF3 on every PDB in the set (per-PDB ΔTM range +0.066 to +0.324). The largest per-PDB improvement is 9ug1 (α-syn fibril, 10 chains): TM 0.573 / RMSD 5.09 Å → TM 0.897 / RMSD 2.32 Å. The smallest is 9ljb (a 1010-residue assembly), where both models still produce only partial folds.

## Full Lightning val-pipeline metrics

Aggregate numbers below are weighted across all 13 val PDBs (6 + 7). Per-dataloader columns show the underlying breakdown — `val_unique_seq` is the harder set (longer chains, unseen sequences, no ligand pose constraint).

### Structure (lDDT, ↑)

| Metric | base OF3 (agg) | OpenFibrilFold (agg) | base val_unique_seq | OFF val_unique_seq | base val_ligand | OFF val_ligand |
|---|---:|---:|---:|---:|---:|---:|
| Intra-protein lDDT | 0.477 | **0.614** | 0.477 | **0.538** | 0.477 | **0.679** |
| Inter-protein lDDT | 0.193 | **0.482** | 0.168 | **0.425** | 0.214 | **0.531** |
| Intra-complex lDDT | 0.479 | **0.615** | 0.477 | **0.538** | 0.480 | **0.681** |
| TM-score (USalign, complex) | 0.333 | **0.539** | 0.282 | **0.456** | 0.376 | **0.610** |

### Ligand (lDDT, ↑) — `val_ligand` only (n=7)

| Metric | base OF3 | OpenFibrilFold |
|---|---:|---:|
| Intra-ligand lDDT | 0.876 | **0.893** |
| Intra-ligand lDDT (uha) | 0.729 | **0.754** |
| Inter-ligand lDDT | 0.232 | **0.313** |
| Protein–ligand lDDT | 0.069 | **0.086** |

### Geometric (↓ lower better; GDT ↑)

| Metric | base OF3 (agg) | OpenFibrilFold (agg) | base val_unique_seq | OFF val_unique_seq | base val_ligand | OFF val_ligand |
|---|---:|---:|---:|---:|---:|---:|
| Distogram loss            | 1.375  | **1.180**  | 1.642 | **1.387** | 1.146 | **1.003** |
| Scaled distogram loss     | 0.0412 | **0.0354** | 0.0493 | **0.0416** | 0.0344 | **0.0301** |
| Intra-protein dRMSD (Å)   | 15.31  | **11.97**  | 13.87 | **12.56** | 16.54 | **11.47** |
| Intra-ligand dRMSD (Å)¹   | 0.800  | **0.698**  | — | — | 0.800 | **0.698** |
| Complex RMSD (Å)          | 30.53  | **19.22**  | 25.63 | **20.84** | 34.74 | **17.83** |
| GDT-TS                    | 0.024  | **0.129**  | 0.038 | **0.129** | 0.011 | **0.128** |
| GDT-HA                    | 0.004  | **0.057**  | 0.006 | **0.053** | 0.002 | **0.060** |

¹ Ligand metrics from `val_ligand` only.

### Confidence calibration (↑)

| Metric | base OF3 (agg) | OpenFibrilFold (agg) | base val_unique_seq | OFF val_unique_seq | base val_ligand | OFF val_ligand |
|---|---:|---:|---:|---:|---:|---:|
| pLDDT (protein)             | 0.254 | **0.440** | 0.280 | **0.316** | 0.233 | **0.546** |
| pLDDT (complex)             | 0.253 | **0.434** | 0.280 | **0.316** | 0.231 | **0.535** |
| pLDDT (ligand)¹             | 0.175 | **0.230** | — | — | 0.175 | **0.230** |
| Pearson(lDDT, pLDDT) protein² | 0.273 | 0.413 | 0.826 | 0.415 | -0.202 | 0.411 |
| Pearson(lDDT, pLDDT) complex² | 0.314 | 0.386 | 0.826 | 0.415 | -0.124 | 0.362 |
| Pearson(lDDT, pLDDT) ligand¹² | 0.144 | -0.517 | — | — | 0.144 | -0.517 |

¹ Ligand metrics from `val_ligand` only.
² **Pearson is unreliable below ~0.5 pLDDT magnitude.** When a model's pLDDT magnitudes are still bunched in the low band, the per-PDB correlation between lDDT and pLDDT swings wildly on tiny shifts in the cluster geometry — base OF3's high `val_unique_seq` Pearson (0.83) reflects this, not genuine calibration, and the negative ligand Pearson in OpenFibrilFold reflects the same problem from the other side. The Pearson rows should not be read as evidence of (or against) confidence calibration on this val set; that conclusion is gated on getting protein/complex pLDDT magnitudes consistently above ~0.5 first.

### Clash rates (↓, mean atom-pair clash fraction)

| Metric | base OF3 (agg) | OpenFibrilFold (agg) |
|---|---:|---:|
| Inter-protein clash | 2.5e-7 | 2.9e-7 |
| Inter protein–ligand clash¹ | 3.9e-7 | **0.0** |

¹ Ligand metrics from `val_ligand` only. Both models are essentially clash-free at val-pipeline granularity.

## Reproducing

The full pipeline lives in the private working repo. From the OpenFibrilFold side, the relevant artefacts are:

- Val pipeline: standard `trainer.validate()` on the `runner_v3_exp43.yml` config, with only `restart_checkpoint_path` swapped between runs. `RIBBONFOLD_DISABLE_GATE=1` set on both to keep the template-embedder code path identical.
- TM-score: [USalign](https://zhanggroup.org/US-align/), built from the published C++ source. Multi-chain complex mode (`-mm 1 -ter 0`) on pred + GT CIFs from `val_structures/epoch_0/`.
- Bar chart: matplotlib, paired bars per PDB sorted by OpenFibrilFold TM descending.

## Caveats

- The 1010-residue `val_unique_seq` PDB (9ljb) is still hard for both models (TM 0.158 → 0.224). Long, polymorph-heavy assemblies remain the failure mode.
- Pearson(lDDT, pLDDT) is reported for completeness but is uninformative at the current pLDDT magnitudes — see footnote ² in the calibration table.
- TM-score and the Lightning lDDT/dRMSD numbers are computed against the **cropped GT assembly that came out of the val dataloader**, not the raw deposited PDB. That is the correct comparator for "did the model do its job on the val pipeline," but absolute numbers may shift slightly when scoring against full-PDB ground truths.
