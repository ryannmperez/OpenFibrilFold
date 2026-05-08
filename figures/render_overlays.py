# PyMOL rendering script for OpenFibrilFold head-to-head overlays.
#
# For each PDB, loads the experimental ground-truth CIF, the OpenFold3 base
# prediction CIF, and the OpenFibrilFold prediction CIF. Aligns both
# predictions onto the experimental backbone using cmd.super (structure-based;
# robust to (a) the GT containing only resolved residues while predictions are
# full-sequence, and (b) homo-multimer fibrils where sequence-based align
# matches a degenerate subset of chains). Renders a 1600×1200 PNG with all
# three structures overlaid:
#   ground truth      grey   (translucent cartoon)
#   OpenFold3 base    orange (cartoon, ligands as sticks)
#   OpenFibrilFold    blue   (cartoon, ligands as sticks)
#
# Usage:
#     pymol -cq render_overlays.py -- \
#         GT_PATH BASE_PRED_PATH FF_PRED_PATH OUT_PNG
#
# Tested with PyMOL 2.5+.

from pymol import cmd
import os
import sys

if len(sys.argv) < 5:
    print("usage: pymol -cq render_overlays.py -- GT BASE_PRED FF_PRED OUT_PNG")
    sys.exit(1)

gt_path, base_path, ff_path, out_png = sys.argv[1:5]

cmd.reinitialize()
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("ray_shadows", 0)
cmd.set("ambient", 0.20)
cmd.set("specular", 0.30)
cmd.set("antialias", 2)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("stick_radius", 0.18)

# ── load ──────────────────────────────────────────────────────────────────
cmd.load(gt_path,   "gt")
cmd.load(base_path, "base")
cmd.load(ff_path,   "ff")

# Predicted/ground-truth CIFs lack HELIX/SHEET records; assign secondary
# structure now so cartoons render as colored arrows/sheets, not pale tubes.
for obj in ("gt", "base", "ff"):
    cmd.dss(obj)

# ── align both predictions onto the experimental backbone ────────────────
# cmd.super does structure-based superposition (no sequence requirement) and
# handles homo-multimer fibrils + GT/pred length mismatch correctly. If super
# fails to converge (n_cycles == 1, no outliers rejected), fall back to
# cealign on a single best-matching chain pair so the picture still renders
# something interpretable rather than the unconverged identity transform.
def superpose(mover, anchor):
    sel_m = f"{mover} and polymer.protein and name CA"
    sel_a = f"{anchor} and polymer.protein and name CA"
    if cmd.count_atoms(sel_m) == 0 or cmd.count_atoms(sel_a) == 0:
        print(f"  [warn] no Cα in {mover} or {anchor}; skipping")
        return
    res = cmd.super(sel_m, sel_a)
    rmsd, n_aln, n_cycles = res[0], res[1], res[2]
    print(f"  super({mover}->{anchor}): rmsd={rmsd:.2f}  n_aln={n_aln}  n_cycles={n_cycles}")
    if n_cycles <= 1:
        # Unconverged. Pick the best chain-pair via cealign.
        best = None
        for cm in cmd.get_chains(mover):
            for ca in cmd.get_chains(anchor):
                try:
                    r = cmd.cealign(f"{anchor} and chain {ca} and polymer.protein",
                                    f"{mover} and chain {cm} and polymer.protein")
                    if best is None or r["RMSD"] < best[0]:
                        best = (r["RMSD"], cm, ca, r["alignment_length"])
                except Exception:
                    pass
        if best:
            print(f"  fallback cealign best: chain {best[1]}->{best[2]}  "
                  f"rmsd={best[0]:.2f}  len={best[3]}")
            cmd.cealign(f"{anchor} and chain {best[2]} and polymer.protein",
                        f"{mover} and chain {best[1]} and polymer.protein")

superpose("base", "gt")
superpose("ff",   "gt")

# ── style ─────────────────────────────────────────────────────────────────
cmd.hide("everything")
cmd.show("cartoon", "all")

# Ligands and any non-polymer hetero atoms (e.g. metal ions) as sticks/spheres
cmd.show("sticks",  "all and not polymer and not resn HOH")
cmd.show("spheres", "all and symbol Mg+Ca+Zn+Fe+Mn+Na+K")
cmd.set("sphere_scale", 0.4)

# Colors — built-in named colors are most reliable across PyMOL builds
cmd.color("grey70",   "gt")
cmd.color("orange",   "base")
cmd.color("marine",   "ff")

# Carbon-only color for ligand sticks so heteroatoms keep CPK
cmd.color("grey70", "gt   and not polymer and elem C")
cmd.color("orange", "base and not polymer and elem C")
cmd.color("marine", "ff   and not polymer and elem C")

# GT is the reference structure. Render as a thick, dark, opaque cartoon-putty
# tube so it pops against the (larger, full-sequence) predictions, which we
# render semi-transparent so all three layers read where they overlap.
cmd.color("grey20", "gt")
cmd.color("grey20", "gt and not polymer and elem C")
cmd.cartoon("putty", "gt and polymer.protein")
cmd.set("cartoon_putty_radius",       0.6, "gt")
cmd.set("cartoon_putty_quality",      20,  "gt")
cmd.set("cartoon_putty_scale_min",    1.0, "gt")
cmd.set("cartoon_putty_scale_max",    1.0, "gt")
cmd.set("cartoon_putty_scale_power",  0.0, "gt")
cmd.set("cartoon_transparency", 0.0, "gt")
cmd.set("cartoon_transparency", 0.5, "base")
cmd.set("cartoon_transparency", 0.5, "ff")

# ── view ──────────────────────────────────────────────────────────────────
# Orient on GT so the experimental fibril sets the camera, but zoom on all
# objects so misaligned predictions don't clip the frame edges.
cmd.orient("gt and polymer.protein")
cmd.zoom("all", 4)

os.makedirs(os.path.dirname(out_png), exist_ok=True)
cmd.ray(1600, 1200)
cmd.png(out_png, dpi=200)
print(f"[wrote] {out_png}")
