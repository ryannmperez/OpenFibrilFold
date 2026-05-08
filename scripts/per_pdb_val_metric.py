"""Compute per-PDB lDDT using OpenFold3's lddt() function (same formula as
the val pipeline) on saved val CIFs, plus a per-residue breakdown for
visualization.

Usage:
    pymol -cq scripts/per_pdb_val_metric.py -- PRED_CIF GT_CIF OUT_JSON

What it does
------------
1. Load pred and gt CIFs in PyMOL.
2. cmd.super for global structural alignment (handles full-seq pred vs
   resolved-only gt).
3. Hungarian assignment on per-chain post-super RMSD to resolve the
   chain-permutation ambiguity for homo-multimer fibrils.
4. Build matched atom tensors (heavy atoms only, matching val pipeline).
5. Call openfold3.core.metrics.validation_all_atom.lddt() with
   intra-chain only filter (matches val/lddt_intra_protein).
6. Extract per-atom lDDT, aggregate to per-residue (mean over heavy atoms
   per residue) for use as a coloring channel.
7. Write a JSON: {scalar, per_residue: {chain: {resi: lddt}}}.
"""

import json
import sys
from pathlib import Path

import numpy as np
import torch
from pymol import cmd

# Use the exact lddt function the val pipeline uses.
from openfold3.core.metrics.validation_all_atom import lddt as of3_lddt

if len(sys.argv) < 4:
    print("usage: pymol -cq per_pdb_val_metric.py -- PRED_CIF GT_CIF OUT_JSON")
    sys.exit(2)

pred_cif, gt_cif, out_json = sys.argv[1], sys.argv[2], sys.argv[3]

cmd.reinitialize()
cmd.load(pred_cif, "pred")
cmd.load(gt_cif,   "gt")

# Global super: aligns pred to gt by structure (necessary because pred
# is full-seq and gt is resolved-only).
cmd.super("pred and polymer.protein", "gt and polymer.protein")


def collect_atoms(obj):
    """Return {chain: {(resi, atom_name): xyz}} for every protein heavy
    atom (skip H), plus a chain → list of (resi, atom_name) order map.

    Iterates a single get_model and groups in Python — sidesteps PyMOL's
    chain-selection-language bug with multi-digit chain IDs.
    """
    out = {}
    m = cmd.get_model(f"{obj} and polymer.protein and not elem H")
    for at in m.atom:
        try:
            resi = int(at.resi)
        except ValueError:
            continue
        out.setdefault(at.chain, {})[(resi, at.name)] = np.array(at.coord)
    return out


pred_atoms = collect_atoms("pred")
gt_atoms   = collect_atoms("gt")

# ── chain permutation: Hungarian on per-chain RMSD ───────────────────────
gt_keys   = list(gt_atoms.keys())
pred_keys = list(pred_atoms.keys())
n_gt, n_pred = len(gt_keys), len(pred_keys)

cost = np.full((n_gt, n_pred), 1e6, dtype=float)
for i, gch in enumerate(gt_keys):
    gd = gt_atoms[gch]
    for j, pch in enumerate(pred_keys):
        pd = pred_atoms[pch]
        # Use Cα-only for the cost matrix (stable, fast).
        common = [k for k in gd if k in pd and k[1] == "CA"]
        if not common:
            continue
        a = np.stack([gd[k] for k in common])
        b = np.stack([pd[k] for k in common])
        cost[i, j] = float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))

from scipy.optimize import linear_sum_assignment
ri, ci = linear_sum_assignment(cost)
chain_map = {gt_keys[i]: pred_keys[j] for i, j in zip(ri, ci)}

# ── build aligned atom tensors ───────────────────────────────────────────
# For each (gt_chain, mapped_pred_chain), iterate gt atoms and pull the
# matching pred atom by (resi, atom_name).
gt_pos, pred_pos, asym_id, atom_resi, atom_chain = [], [], [], [], []
asym_counter = 0
for gch, pch in chain_map.items():
    gd = gt_atoms[gch]
    pd = pred_atoms[pch]
    for (resi, name), gxyz in gd.items():
        if (resi, name) in pd:
            gt_pos.append(gxyz)
            pred_pos.append(pd[(resi, name)])
            asym_id.append(asym_counter)
            atom_resi.append(resi)
            atom_chain.append(gch)
    asym_counter += 1

if not gt_pos:
    Path(out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(out_json).write_text(json.dumps({
        "lddt_intra_complex": float("nan"),
        "n_atoms_aligned": 0,
        "per_residue": {},
    }))
    sys.exit(0)

gt_pos   = torch.tensor(np.stack(gt_pos),   dtype=torch.float32)
pred_pos = torch.tensor(np.stack(pred_pos), dtype=torch.float32)
asym_id_t = torch.tensor(asym_id, dtype=torch.long)
N = gt_pos.shape[0]
atom_mask  = torch.ones(N, dtype=torch.float32)
intra_filt = torch.ones(N, dtype=torch.float32)
inter_filt = torch.zeros(N, N, dtype=torch.float32)  # intra-only metric

# Pairwise distances.
gt_pair   = torch.cdist(gt_pos.unsqueeze(0),   gt_pos.unsqueeze(0)  ).squeeze(0)
pred_pair = torch.cdist(pred_pos.unsqueeze(0), pred_pos.unsqueeze(0)).squeeze(0)

# ── call OF3's lddt() — same fn the val pipeline uses ────────────────────
intra_score, _ = of3_lddt(
    pair_dist_pred_pos=pred_pair,
    pair_dist_gt_pos=gt_pair,
    all_atom_mask=atom_mask,
    intra_mask_filter=intra_filt,
    inter_mask_filter=inter_filt,
    asym_id=asym_id_t,
)
lddt_intra_complex = float(intra_score.item()) if intra_score is not None else float("nan")

# ── per-atom (and per-residue) lDDT for visualization ────────────────────
# Replicate lddt's inner computation but keep per-atom granularity.
cutoff = 15.0
thresholds = (0.5, 1.0, 2.0, 4.0)
intra_mask = (asym_id_t.unsqueeze(-1) == asym_id_t.unsqueeze(-2))
dists_to_score = (gt_pair < cutoff) & (~torch.eye(N, dtype=torch.bool)) \
                 & intra_mask & (atom_mask.unsqueeze(-1).bool()) \
                 & (atom_mask.unsqueeze(-2).bool())
err = torch.abs(gt_pair - pred_pair)
score = torch.zeros_like(err)
for t in thresholds:
    score += (err < t).float()
score /= len(thresholds)
denom = dists_to_score.float().sum(dim=-1).clamp(min=1.0)
per_atom_lddt = (score * dists_to_score.float()).sum(dim=-1) / denom
# Atoms with no in-cutoff partners get nan.
per_atom_lddt[dists_to_score.sum(dim=-1) == 0] = float("nan")

# Aggregate to per-residue (mean over the residue's heavy atoms).
per_res_sum, per_res_cnt = {}, {}
for k, v in zip(zip(atom_chain, atom_resi), per_atom_lddt.tolist()):
    if v == v:  # skip nan
        per_res_sum[k] = per_res_sum.get(k, 0.0) + v
        per_res_cnt[k] = per_res_cnt.get(k, 0)   + 1
per_residue = {}
for (ch, resi), s in per_res_sum.items():
    per_residue.setdefault(ch, {})[str(resi)] = s / per_res_cnt[(ch, resi)]

out = {
    "lddt_intra_complex": lddt_intra_complex,
    "n_atoms_aligned": int(N),
    "n_chains_matched": len(chain_map),
    "chain_map": chain_map,
    "per_residue": per_residue,
}
Path(out_json).parent.mkdir(parents=True, exist_ok=True)
Path(out_json).write_text(json.dumps(out, indent=2))
print(f"lddt_intra_complex={lddt_intra_complex:.4f}  "
      f"n_atoms={N}  n_chains_matched={len(chain_map)}  → {out_json}")
