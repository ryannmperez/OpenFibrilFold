"""Per-PDB TM-score for a saved val pair (pred + gt CIF).

Method: US-align with multi-chain protein flags, the same invocation
RibbonFold uses (Mingchenchen/RibbonFold, src/metrics.py:404):

    USalign {pred} {gt} -mol prot -mm 1 -ter 1

US-align finds the optimal chain mapping (handles homo-oligomeric fibril
chain-permutation ambiguity), then iteratively maximises TM-score over
the alignment. Output reports two TM-scores; we take the one normalised
by Structure_2 (the GT), per US-align's own recommendation
("You should use TM-score normalized by length of the reference structure").

Usage:
    python scripts/per_pdb_tm_score.py PRED_CIF GT_CIF [OUT_JSON]

USalign discovery: $USALIGN env var, then PATH.
Source: https://zhanggroup.org/US-align/
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


def find_usalign():
    env = os.environ.get("USALIGN")
    if env and os.path.exists(env):
        return env
    p = shutil.which("USalign")
    if p:
        return p
    sys.exit("USalign not found. Set $USALIGN or install on PATH "
             "(https://zhanggroup.org/US-align/).")


def parse_usalign(text):
    out = {}
    for line in text.splitlines():
        if line.startswith("Length of Structure_1:"):
            out["L_pred"] = int(line.split(":")[1].split()[0])
        elif line.startswith("Length of Structure_2:"):
            out["L_gt"] = int(line.split(":")[1].split()[0])
        elif line.startswith("Aligned length="):
            parts = line.split(",")
            out["L_aligned"] = int(parts[0].split("=")[1].strip())
            out["rmsd"] = float(parts[1].split("=")[1].strip())
        elif "normalized by length of Structure_1" in line:
            out["tm_norm_pred"] = float(line.split("=")[1].split("(")[0].strip())
        elif "normalized by length of Structure_2" in line:
            out["tm_norm_gt"] = float(line.split("=")[1].split("(")[0].strip())
    return out


def main():
    if len(sys.argv) < 3:
        print("usage: per_pdb_tm_score.py PRED_CIF GT_CIF [OUT_JSON]")
        sys.exit(2)
    pred, gt = sys.argv[1], sys.argv[2]
    out_json = sys.argv[3] if len(sys.argv) > 3 else None

    cmd = [find_usalign(), pred, gt, "-mol", "prot", "-mm", "1", "-ter", "1"]
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    parsed = parse_usalign(res.stdout)
    parsed["pred_cif"] = pred
    parsed["gt_cif"] = gt
    parsed["cmd"] = " ".join(cmd)

    print(f"tm_score (norm GT) = {parsed.get('tm_norm_gt'):.4f}  "
          f"RMSD = {parsed.get('rmsd'):.2f} Å  "
          f"L_aligned = {parsed.get('L_aligned')}  "
          f"L_pred = {parsed.get('L_pred')}  "
          f"L_gt = {parsed.get('L_gt')}")
    if out_json:
        Path(out_json).parent.mkdir(parents=True, exist_ok=True)
        Path(out_json).write_text(json.dumps(parsed, indent=2))


main()
