"""Assemble the 3×3 head-to-head grid from the per-panel PNGs.

Layout:
    cols (left→right): 9qlu (best val_unique_seq), 9ug1 (best val_ligand),
                       9ljb (hardest case in either val set)
    rows (top→bottom): OF3 base, Ground truth, OpenFibrilFold

Each cell shows the panel image and lDDT + TM-score labels below.
Row labels appear as rotated text on the left margin; column labels
across the top.

Column selection: per-PDB FT-vs-OF3 TM-score delta. 9qlu has the largest
ΔTM (+0.253) within val_unique_seq (rare-fold apo / ligand fibrils on
unseen sequences). 9ug1 has the largest ΔTM (+0.324) within val_ligand.
9ljb has the smallest ΔTM (+0.066) and the lowest absolute TM in either
dataloader — a 1010-residue assembly where both models still struggle.
"""

from PIL import Image, ImageDraw, ImageFont, ImageChops
import os

PANELS_DIR = os.path.join(os.path.dirname(__file__), "panels")
OUT = os.path.join(os.path.dirname(__file__), "h2h_grid.png")

PDBS = ["9qlu", "9ug1", "9ljb"]
ROWS = [
    ("OF3 base",        "base", (227, 121, 21)),
    ("Ground truth",    "gt",   (50, 50, 50)),
    ("OpenFibrilFold",  "off",  (28, 99, 173)),
]

# Per-PDB lDDT computed via openfold3.core.metrics.validation_all_atom.lddt
# (the function the val pipeline uses) on heavy atoms after Hungarian
# chain-permutation matching. Same formula → same numerical scale as the
# README metric table. Computed by scripts/per_pdb_val_metric.py on the
# val_compare pred / gt CIFs (epoch 0 of trainer.validate on exp43
# step 1024 vs base OF3 ft3_v1).
LDDT = {
    "9qlu": {"base": 0.406, "off": 0.492},
    "9ug1": {"base": 0.436, "off": 0.774},
    "9ljb": {"base": 0.395, "off": 0.400},
}

# Per-PDB TM-score from US-align in multi-chain complex mode
# (USalign … -mm 1 -ter 0), so the score covers the full fibril assembly
# with chain correspondences found by MM-align greedy search.
# Computed by scripts/compute_val_compare_tm.py on the same val_compare
# pred / gt CIFs used for LDDT above.
TM = {
    "9qlu": {"base": 0.204, "off": 0.457},
    "9ug1": {"base": 0.573, "off": 0.897},
    "9ljb": {"base": 0.158, "off": 0.224},
}

CELL = 600                # input panel size (matches render_panel.py)
CAPTION_H = 96            # space below each panel for the two-line caption
COL_HEADER_H = 70         # space above for PDB IDs
ROW_LABEL_W = 90          # left margin for row labels
PAD = 12                  # inner padding around each cell

def find_font(size):
    for path in [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
        "/usr/share/fonts/TTF/DejaVuSans-Bold.ttf",
    ]:
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()

font_label   = find_font(38)
font_header  = find_font(40)
font_caption = find_font(28)

cell_w = CELL + 2 * PAD
cell_h = CELL + CAPTION_H + 2 * PAD
W = ROW_LABEL_W + 3 * cell_w
H = COL_HEADER_H + 3 * cell_h

img = Image.new("RGB", (W, H), "white")
draw = ImageDraw.Draw(img)

# Column headers
for ci, pdb in enumerate(PDBS):
    x = ROW_LABEL_W + ci * cell_w + cell_w // 2
    bbox = draw.textbbox((0, 0), pdb, font=font_header)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    draw.text((x - tw // 2, (COL_HEADER_H - th) // 2),
              pdb, fill="black", font=font_header)

# Row labels (rotated 90° CCW so text reads bottom-up on the left margin)
for ri, (label, _, color) in enumerate(ROWS):
    txt_img = Image.new("RGBA", (cell_h, ROW_LABEL_W), (0, 0, 0, 0))
    td = ImageDraw.Draw(txt_img)
    bbox = td.textbbox((0, 0), label, font=font_label)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    td.text(((cell_h - tw) // 2, (ROW_LABEL_W - th) // 2 - 4),
            label, fill=color, font=font_label)
    rot = txt_img.rotate(90, expand=True)
    img.paste(rot, (0, COL_HEADER_H + ri * cell_h), rot)

def trim_to_content(im, margin=20):
    """Crop white margins, then re-pad to a square so the structure fills."""
    bg = Image.new(im.mode, im.size, (255, 255, 255))
    diff = ImageChops.difference(im, bg)
    bbox = diff.getbbox()
    if bbox is None:
        return im
    cropped = im.crop(bbox)
    side = max(cropped.size) + 2 * margin
    sq = Image.new("RGB", (side, side), "white")
    sq.paste(cropped, ((side - cropped.size[0]) // 2,
                       (side - cropped.size[1]) // 2))
    return sq

# Pre-load panels and normalise scale per column: each column's cells are
# resized to the largest content extent in that column, so OF3 / GT / OFF
# panels for the same PDB share a visual scale.
loaded = {}
content_side = {}
for ci, pdb in enumerate(PDBS):
    cells = {}
    max_side = 0
    for _, src, _ in ROWS:
        path = os.path.join(PANELS_DIR, f"{pdb}_{src}.png")
        cropped = trim_to_content(Image.open(path).convert("RGB"))
        cells[src] = cropped
        max_side = max(max_side, cropped.size[0])
    loaded[pdb] = cells
    content_side[pdb] = max_side

# Panels
for ri, (label, src, color) in enumerate(ROWS):
    for ci, pdb in enumerate(PDBS):
        cropped = loaded[pdb][src]
        # Pad cropped to the column's max side, then resize to CELL.
        s = content_side[pdb]
        sq = Image.new("RGB", (s, s), "white")
        sq.paste(cropped, ((s - cropped.size[0]) // 2,
                           (s - cropped.size[1]) // 2))
        panel = sq.resize((CELL, CELL), Image.LANCZOS)
        x = ROW_LABEL_W + ci * cell_w + PAD
        y = COL_HEADER_H + ri * cell_h + PAD
        img.paste(panel, (x, y))

        # Caption: GT row is the reference (single line); prediction rows
        # get two stacked lines — intra-complex lDDT and TM-score.
        if src == "gt":
            cap_lines = ["experimental reference"]
        else:
            cap_lines = [
                f"intra-complex lDDT = {LDDT[pdb][src]:.3f}",
                f"TM-score = {TM[pdb][src]:.3f}",
            ]
        cx = ROW_LABEL_W + ci * cell_w + cell_w // 2
        # Stack lines vertically inside the caption area.
        line_h = font_caption.getbbox("Ag")[3] - font_caption.getbbox("Ag")[1]
        block_h = line_h * len(cap_lines) + 6 * (len(cap_lines) - 1)
        y0 = COL_HEADER_H + ri * cell_h + PAD + CELL + (CAPTION_H - block_h) // 2
        for i, line in enumerate(cap_lines):
            bbox = draw.textbbox((0, 0), line, font=font_caption)
            tw = bbox[2] - bbox[0]
            draw.text((cx - tw // 2, y0 + i * (line_h + 6)),
                      line, fill=color, font=font_caption)

img.save(OUT, optimize=True)
print(f"[wrote] {OUT}  ({img.size[0]}×{img.size[1]})")
