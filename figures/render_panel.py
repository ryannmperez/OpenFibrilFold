# Render a single structure top-down (looking down the protofilament
# stacking axis). Used to build the 3×3 head-to-head figure.
#
# Approach: pick the principal direction along which chain centroids
# project most evenly (the stacking axis — a real fibril has one chain
# layer per period). Rotate the structure so that axis lies along world
# +z, then PyMOL's default camera gives the end-on view.
#
# Optional residue highlighting: pass a per-residue lDDT JSON
# (chain → resi → score) as a 4th argument. Residues are coloured by
# blending the base colour with white based on score (high = saturated,
# low = faded), which makes "the part that aligned to GT" pop visually.
#
# Optional reference alignment: pass --ref REF_CIF to super the input
# onto the reference and use the reference's stacking axis for the
# rotation. Without this, every panel picks its own PCA-derived axis;
# PCA's sign ambiguity then lets pred and GT panels land on opposite
# faces of the cross-section. With --ref REF (= GT cif), pred panels
# share the GT's frame and end-on orientation.
#
# Usage:
#     pymol -cq render_panel.py -- INPUT_CIF OUT_PNG COLOR \
#         [PER_RES_JSON] [--ref REF_CIF]
#
# COLOR is a PyMOL named color (e.g. "orange", "marine", "grey50").
# Renders 600×600 PNG, no transparency, white background.

from pymol import cmd
import argparse
import os, sys, json
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("input")
ap.add_argument("out_png")
ap.add_argument("color")
ap.add_argument("per_res_json", nargs="?", default=None)
ap.add_argument("--ref", default=None,
                help="reference CIF; obj is super'd to it and the rotation "
                     "is computed from ref's stacking axis so panels share "
                     "an end-on orientation")
args = ap.parse_args()

inp, out_png, color = args.input, args.out_png, args.color
per_res_json = args.per_res_json
ref_cif = args.ref

cmd.reinitialize()
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("ray_shadows", 0)
cmd.set("ambient", 0.25)
cmd.set("specular", 0.30)
cmd.set("antialias", 2)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("stick_radius", 0.18)

cmd.load(inp, "obj")
cmd.dss("obj")

if ref_cif:
    cmd.load(ref_cif, "ref")
    cmd.dss("ref")
    # Manual chain-matched Kabsch alignment.
    # cmd.super gets confused on homo-oligomer fibrils (every chain has
    # the same sequence, so the seq-then-struct heuristic aliases chains
    # and converges to a 20+ Å fit). We instead alternate between Kabsch
    # (given a chain mapping, solve for rotation) and Hungarian (given
    # a rotation, solve for chain mapping). Two iterations are enough
    # to recover the right protofilament-to-protofilament correspondence
    # for multi-protofilament fibrils where pred chain "1" might map to
    # gt chain "6" instead of gt chain "1".
    def _ca_by_chain(obj):
        m = cmd.get_model(f"{obj} and polymer.protein and name CA")
        out = {}
        for at in m.atom:
            try:
                resi = int(at.resi)
            except ValueError:
                continue
            out.setdefault(at.chain, {})[resi] = np.array(at.coord)
        return out

    def _kabsch(P, Q):
        cP, cQ = P.mean(axis=0), Q.mean(axis=0)
        H = (P - cP).T @ (Q - cQ)
        U, _, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
        t = cQ - R @ cP
        return R, t

    def _matched_atoms(obj_ca, ref_ca, mapping):
        P, Q = [], []
        for och, rch in mapping.items():
            ofd, rfd = obj_ca[och], ref_ca[rch]
            for resi, oxyz in ofd.items():
                if resi in rfd:
                    P.append(oxyz)
                    Q.append(rfd[resi])
        return np.array(P), np.array(Q)

    obj_ca = _ca_by_chain("obj")
    ref_ca = _ca_by_chain("ref")
    obj_chains = list(obj_ca.keys())
    ref_chains = list(ref_ca.keys())

    # Initial mapping: by chain ID (works for 1-protofilament fibrils;
    # may be wrong for multi-protofilament cases — Hungarian below fixes it).
    mapping = {c: c for c in obj_chains if c in ref_ca}

    if len(mapping) >= 1:
        from scipy.optimize import linear_sum_assignment
        for _ in range(3):
            P, Q = _matched_atoms(obj_ca, ref_ca, mapping)
            if len(P) < 3:
                break
            R, t = _kabsch(P, Q)

            # Re-evaluate chain mapping on transformed obj centroids.
            obj_cent = {ch: R @ np.mean(list(d.values()), axis=0) + t
                        for ch, d in obj_ca.items()}
            ref_cent = {ch: np.mean(list(d.values()), axis=0)
                        for ch, d in ref_ca.items()}
            cost = np.array([[np.linalg.norm(obj_cent[oc] - ref_cent[rc])
                              for rc in ref_chains] for oc in obj_chains])
            ri, ci = linear_sum_assignment(cost)
            new_mapping = {obj_chains[i]: ref_chains[j] for i, j in zip(ri, ci)}
            if new_mapping == mapping:
                break
            mapping = new_mapping

        # Final Kabsch with the converged mapping.
        P, Q = _matched_atoms(obj_ca, ref_ca, mapping)
        R, t = _kabsch(P, Q)
        cmd.transform_object("obj", [
            float(R[0, 0]), float(R[0, 1]), float(R[0, 2]), float(t[0]),
            float(R[1, 0]), float(R[1, 1]), float(R[1, 2]), float(t[1]),
            float(R[2, 0]), float(R[2, 1]), float(R[2, 2]), float(t[2]),
            0.0, 0.0, 0.0, 1.0,
        ])

# ── style ─────────────────────────────────────────────────────────────────
# Style the rendered object only — ref (if loaded) stays hidden, so it
# only contributes its centroids/transform and never paints into the PNG.
cmd.hide("everything")
cmd.show("cartoon", "obj")
cmd.show("sticks",  "obj and not polymer and not resn HOH")
cmd.show("spheres", "obj and symbol Mg+Ca+Zn+Fe+Mn+Na+K")
cmd.set("sphere_scale", 0.4)
cmd.color(color, "obj")
cmd.color(color, "obj and not polymer and elem C")

# ── per-residue lDDT highlighting ────────────────────────────────────────
# Residues that aligned to GT are bold, mis-aligned residues fade. The
# mapping is a sharp sigmoid centred at lDDT = 0.5: above ~0.6 → almost
# full base colour; below ~0.4 → almost invisible (95% white). This
# makes the locally-correct vs locally-wrong regions read clearly.
# Cartoons of low-lDDT residues are also dropped to a thin loop so the
# saturated regions pop visually.
if per_res_json and os.path.exists(per_res_json):
    with open(per_res_json) as f:
        per_res = json.load(f).get("per_residue", {})

    base_rgb = np.array(cmd.get_color_tuple(cmd.get_color_index(color)))
    used_colors = set()
    bad_sel_parts = []   # residues that get a thinner cartoon
    for ch, residues in per_res.items():
        for resi_str, l in residues.items():
            t = 1.0 / (1.0 + np.exp(-12.0 * (l - 0.5)))
            blend = base_rgb * t + np.array([1.0, 1.0, 1.0]) * (1.0 - t)
            bucket = int(round(t * 100))
            cname = f"_l_{bucket}"
            if cname not in used_colors:
                cmd.set_color(cname, [float(blend[0]),
                                      float(blend[1]),
                                      float(blend[2])])
                used_colors.add(cname)
            cmd.color(cname, f"obj and chain {ch} and resi {resi_str}")
            if t < 0.3:
                bad_sel_parts.append(f"(chain {ch} and resi {resi_str})")

    if bad_sel_parts:
        bad_sel = "obj and (" + " or ".join(bad_sel_parts) + ")"
        cmd.set("cartoon_transparency", 0.55, bad_sel)

# ── compute stacking axis from chain centroids ───────────────────────────
# Always use obj's own centroids — for predictions whose actual stacking
# axis differs structurally from gt's (e.g. a 1-protofilament prediction
# of a 2-protofilament gt), forcing gt's axis would render obj side-on.
# When --ref is supplied, ref's stacking axis is computed too and used
# only to break PCA's sign ambiguity (so pred and gt agree on which
# face of the stack is "front" once Kabsch has put them in a common frame).
def _centroids(obj_name):
    mm = cmd.get_model(f"{obj_name} and polymer.protein and name CA")
    bc = {}
    for at in mm.atom:
        bc.setdefault(at.chain, []).append(at.coord)
    return np.array([np.mean(np.array(pts), axis=0) for pts in bc.values()])

centroids = _centroids("obj")
ref_centroids = _centroids("ref") if ref_cif else None

# ── physically rotate the object so its stacking axis aligns with +z ────
# Then the default PyMOL camera (looking down +z) gives the end-on view.
#
# Picking the stacking axis: PCA's PC1 of chain centroids isn't reliable
# for multi-protofilament fibrils — for a 2-protofilament structure with
# few chains per filament the cross-protofilament axis can outsize the
# stacking axis. So evaluate all three PCs and pick the one whose
# average per-chain rise is closest to the cross-β period (~4.7 Å).
if len(centroids) >= 2:
    centred = centroids - centroids.mean(axis=0)
    U, S, Vt = np.linalg.svd(centred, full_matrices=False)

    # Score each PC by how evenly the centroids project onto it.
    # The stacking axis projects centroids at uniform layer spacing
    # (every chain = one layer), so its spacing-CV is small.
    # Cross-protofilament axes give bimodal projections (a large gap
    # between protofilament clusters dominates the diff distribution)
    # → high CV. PCs with nearly-zero extent give meaningless tiny
    # spacings → also rejected.
    def score(pts_centred, direction):
        proj = sorted(pts_centred @ direction)
        if proj[-1] - proj[0] < 3.0:
            return 1e9   # axis with no real extent
        diffs = np.array([proj[i+1] - proj[i] for i in range(len(proj)-1)])
        m = float(diffs.mean())
        if m < 1e-3:
            return 1e9
        return float(diffs.std()) / m   # CV — lower is more evenly spaced

    best = min(range(3), key=lambda k: score(centred, Vt[k]))
    stacking = Vt[best] / np.linalg.norm(Vt[best])

    # Sign-match to ref's stacking axis when --ref is supplied — this is
    # what stops 7nck-style PCA-sign-ambiguity from showing pred and gt
    # on opposite faces. Magnitude doesn't matter, only sign of the dot.
    if ref_centroids is not None and len(ref_centroids) >= 2:
        ref_centred = ref_centroids - ref_centroids.mean(axis=0)
        _, _, Vt_r = np.linalg.svd(ref_centred, full_matrices=False)
        ref_best = min(range(3), key=lambda k: score(ref_centred, Vt_r[k]))
        ref_stacking = Vt_r[ref_best] / np.linalg.norm(Vt_r[ref_best])
        if float(np.dot(stacking, ref_stacking)) < 0:
            stacking = -stacking

    target = np.array([0., 0., 1.])
    axis   = np.cross(stacking, target)
    n      = np.linalg.norm(axis)
    if n > 1e-6:
        angle_deg = float(np.degrees(np.arctan2(
            n, float(np.dot(stacking, target))
        )))
        u = axis / n
        com = centroids.mean(axis=0)
        cmd.rotate([float(u[0]), float(u[1]), float(u[2])],
                   angle_deg, "obj", camera=0,
                   origin=[float(com[0]), float(com[1]), float(com[2])])

if ref_cif:
    cmd.delete("ref")  # remove from the rendered scene

cmd.zoom("obj", 2)

os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
cmd.ray(600, 600)
cmd.png(out_png, dpi=200)
print(f"[wrote] {out_png}")
