#!/usr/bin/env python3
"""Build the README metrics table from two metrics.csv files.

Usage:
    build_metrics_table.py --base BASE_METRICS_CSV --tuned TUNED_METRICS_CSV
        [--base-step N] [--tuned-step N]

For each csv, picks the latest row with a non-empty `val/lddt_intra_protein`
unless --*-step is given. Writes a markdown table to stdout.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from typing import Dict, List, Optional, Tuple


GROUPS: List[Tuple[str, List[Tuple[str, str, str]]]] = [
    ("Structure (lDDT, ↑)", [
        ("val/lddt_intra_protein",          "Intra-protein lDDT",         "ratio"),
        ("val/lddt_inter_protein_protein",  "Inter-protein lDDT",         "ratio"),
        ("val/lddt_intra_complex",          "Intra-complex lDDT",         "ratio"),
        ("val/lddt_intra_modified_residues","Modified residues lDDT",     "ratio"),
    ]),
    ("Ligand (lDDT, ↑)", [
        ("val/lddt_intra_ligand",           "Intra-ligand lDDT",          "ratio"),
        ("val/lddt_intra_ligand_uha",       "Intra-ligand lDDT (uha)",    "ratio"),
        ("val/lddt_inter_ligand_ligand",    "Inter-ligand lDDT",          "ratio"),
        ("val/lddt_inter_protein_ligand",   "Protein–ligand lDDT",        "ratio"),
    ]),
    ("Geometric (↓ lower better; GDT ↑)", [
        ("val/distogram_loss",              "Distogram loss",             "loss"),
        ("val/scaled_distogram_loss",       "Scaled distogram loss",      "loss"),
        ("val/drmsd_intra_protein",         "Intra-protein dRMSD (Å)",    "dist"),
        ("val/drmsd_intra_ligand",          "Intra-ligand dRMSD (Å)",     "dist"),
        ("val/rmsd",                        "Complex RMSD (Å)",           "dist"),
        ("val/gdt_ts",                      "GDT-TS",                     "ratio"),
        ("val/gdt_ha",                      "GDT-HA",                     "ratio"),
    ]),
    ("Confidence calibration (↑)", [
        ("val/pearson_correlation_lddt_plddt_protein", "Pearson(lDDT, pLDDT) protein", "ratio"),
        ("val/pearson_correlation_lddt_plddt_ligand",  "Pearson(lDDT, pLDDT) ligand",  "ratio"),
        ("val/pearson_correlation_lddt_plddt_complex", "Pearson(lDDT, pLDDT) complex", "ratio"),
    ]),
]


def latest_val_row(path: str, target_step: Optional[int]) -> Dict[str, str]:
    with open(path) as f:
        rows = list(csv.DictReader(f))
    val_rows = [r for r in rows if r.get("val/lddt_intra_protein", "")]
    if not val_rows:
        sys.exit(f"no val rows in {path}")
    if target_step is None:
        return val_rows[-1]
    for r in val_rows:
        try:
            if int(r["step"]) == target_step:
                return r
        except (KeyError, ValueError):
            continue
    sys.exit(f"no val row at step {target_step} in {path}")


def fmt(value: str, kind: str) -> str:
    if not value:
        return "—"
    try:
        f = float(value)
    except ValueError:
        return "—"
    if math.isnan(f):
        return "—"
    if kind == "ratio":
        return f"{f:.4f}"
    if kind == "loss":
        return f"{f:.4f}"
    if kind == "dist":
        return f"{f:.3f}"
    return str(f)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="OF3 base metrics.csv")
    ap.add_argument("--tuned", required=True, help="FibrilFold metrics.csv")
    ap.add_argument("--base-step", type=int, default=None)
    ap.add_argument("--tuned-step", type=int, default=None)
    ap.add_argument("--base-label", default="OpenFold3 base")
    ap.add_argument("--tuned-label", default="FibrilFold")
    args = ap.parse_args()

    base = latest_val_row(args.base, args.base_step)
    tuned = latest_val_row(args.tuned, args.tuned_step)

    out: List[str] = []
    out.append(f"| Metric | {args.base_label} | {args.tuned_label} |")
    out.append("| --- | --- | --- |")
    for header, items in GROUPS:
        out.append(f"| **{header}** |  |  |")
        for col, label, kind in items:
            b = fmt(base.get(col, ""), kind)
            t = fmt(tuned.get(col, ""), kind)
            if b == "—" and t == "—":
                continue
            out.append(f"| {label} | {b} | {t} |")

    print("\n".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
