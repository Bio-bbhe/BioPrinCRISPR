# 06_rna_ss_vis.py with argparse added

from __future__ import division, print_function
from builtins import zip
import forgi.threedee.utilities.vector as ftuv
import math, logging, itertools, colorsys, forgi
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def _clashfree_annot_pos(pos, coords):
    for c in coords:
        dist = ftuv.vec_distance(c, pos)
        if dist < 14:
            return False
    return True


def _annotate_rna_plot(ax, cg, coords, annotations, text_kwargs):
    annot_dict = {elem: elem for elem in cg.defines}
    if annotations is None:
        annot_dict = {elem: "" for elem in cg.defines}
    else:
        annot_dict.update(annotations)
    stem_coords = {}
    for stem in cg.stem_iterator():
        stem_start = np.mean([coords[cg.defines[stem][0] - 1], coords[cg.defines[stem][3] - 1]], axis=0)
        stem_end = np.mean([coords[cg.defines[stem][1] - 1], coords[cg.defines[stem][2] - 1]], axis=0)
        stem_center = np.mean([stem_start, stem_end], axis=0)
        stem_coords[stem] = (stem_start, stem_center, stem_end)
        if annot_dict[stem]:
            stem_vec = stem_end - stem_start
            norm_vec = (stem_vec[1], -stem_vec[0])
            norm_vec /= ftuv.magnitude(norm_vec)
            annot_pos = np.array(stem_center) + 23 * norm_vec
            if not _clashfree_annot_pos(annot_pos, coords):
                annot_pos = np.array(stem_center) - 23 * norm_vec
                if not _clashfree_annot_pos(annot_pos, coords):
                    log.info(f"Cannot annotate {stem} as '{annot_dict[stem]}', insufficient space.")
                    annot_pos = None
            if annot_pos is not None:
                ax.annotate(annot_dict[stem], xy=annot_pos, ha="center", va="center", **text_kwargs)
    for hloop in cg.hloop_iterator():
        hc = [coords[nt - 1] for nt in cg.define_residue_num_iterator(hloop, adjacent=True)]
        annot_pos = np.mean(hc, axis=0)
        if _clashfree_annot_pos(annot_pos, coords):
            ax.annotate(annot_dict[hloop], xy=annot_pos, ha="center", va="center", **text_kwargs)
        else:
            log.info(f"Cannot annotate {hloop} as '{annot_dict[hloop]}' INSIDE, trying outside...")
            nt1, nt2 = cg.define_a(hloop)
            start = np.mean([coords[nt1 - 1], coords[nt2 - 1]], axis=0)
            vec = annot_pos - start
            annot_pos2 = annot_pos + vec * 3
            if _clashfree_annot_pos(annot_pos2, coords):
                ax.annotate(annot_dict[hloop], xy=annot_pos2, ha="center", va="center", **text_kwargs)
            else:
                log.info(f"Cannot annotate {hloop} as '{annot_dict[hloop]}', insufficient space.")
    for iloop in cg.iloop_iterator():
        s1, s2 = cg.connections(iloop)
        annot_pos = np.mean([stem_coords[s1][2], stem_coords[s2][0]], axis=0)
        if _clashfree_annot_pos(annot_pos, coords):
            ax.annotate(annot_dict[iloop], xy=annot_pos, ha="center", va="center", **text_kwargs)
        else:
            log.info(f"Cannot annotate {iloop} as '{annot_dict[iloop]}', trying outside...")
            loop_vec = stem_coords[s2][0] - stem_coords[s1][2]
            norm_vec = (loop_vec[1], -loop_vec[0])
            norm_vec /= ftuv.magnitude(norm_vec)
            annot_pos_p = annot_pos + 25 * norm_vec
            annot_pos_m = annot_pos - 25 * norm_vec
            plus = sum(ftuv.vec_distance(annot_pos_p, coords[nt - 1]) < ftuv.vec_distance(annot_pos_m, coords[nt - 1]) for nt in cg.define_residue_num_iterator(iloop))
            minus = len(list(cg.define_residue_num_iterator(iloop))) - plus
            preferred_pos = annot_pos_p if plus > minus else annot_pos_m
            if _clashfree_annot_pos(preferred_pos, coords):
                ax.annotate(annot_dict[iloop], xy=preferred_pos, ha="center", va="center", **text_kwargs)
            else:
                log.info(f"Cannot annotate {iloop} as '{annot_dict[iloop]}', insufficient space.")
    for mloop in itertools.chain(cg.floop_iterator(), cg.tloop_iterator(), cg.mloop_iterator()):
        nt1, nt2 = cg.define_a(mloop)
        res = list(cg.define_residue_num_iterator(mloop))
        if len(res) == 0:
            anchor = np.mean([coords[nt1 - 1], coords[nt2 - 1]], axis=0)
        elif len(res) % 2 == 1:
            anchor = coords[res[len(res) // 2] -1]
        else:
            anchor = np.mean([coords[res[len(res)//2 - 1] -1], coords[res[len(res)//2] - 1]], axis=0)
        loop_vec = coords[nt1 - 1] - coords[nt2 - 1]
        norm_vec = (loop_vec[1], -loop_vec[0])
        norm_vec /= ftuv.magnitude(norm_vec)
        annot_pos = anchor - norm_vec * 18
        if _clashfree_annot_pos(annot_pos, coords):
            ax.annotate(annot_dict[mloop], xy=annot_pos, ha="center", va="center", **text_kwargs)
        else:
            log.info(f"Cannot annotate {mloop} as '{annot_dict[mloop]}', insufficient space.")


def plot_rna(cg, ax=None, offset=(0, 0), color=True, text_kwargs=None, backbone_kwargs=None,
             basepair_kwargs=None, mismatch_kwargs=None, lighten=0.0, annotations=None):
    if text_kwargs is None:
        text_kwargs = {}
    if backbone_kwargs is None:
        backbone_kwargs = {}
    if annotations is None:
        annotations = {}
    if mismatch_kwargs is None:
        mismatch_kwargs = {}
    if basepair_kwargs is None:
        basepair_kwargs = {}
    import RNA
    import matplotlib.colors as mc
    RNA.cvar.rna_plot_type = 1

    if ax is None:
        ax = plt.gca()
    if offset is None:
        offset = (0, 0)
    elif offset is True:
        offset = (ax.get_xlim()[1], ax.get_ylim()[1])

    bp_string = cg.to_dotbracket_string()
    el_string = cg.to_element_string()
    el_to_color = {'f': 'orange', 't': 'orange', 's': 'green', 'h': 'blue', 'i': 'yellow', 'm': 'red'}

    vrna_coords = RNA.get_xy_coordinates(bp_string)
    coords = np.array([(offset[0] + vrna_coords.get(i).X, offset[1] + vrna_coords.get(i).Y) for i in range(len(bp_string))])

    # Backbone
    bkwargs = {"color": "black", "zorder": 0}
    bkwargs.update(backbone_kwargs)
    ax.plot(coords[:, 0], coords[:, 1], **bkwargs)

    # Basepairs and mismatches
    basepairs, basepairs_bases, mismatched_basepairs = [], [], []
    matched_pairs = [{'A', 'u'}, {'G', 'C'}, {'U', 'A'}, {'C', 'G'}]

    for stem in cg.stem_iterator():
        for p1, p2 in cg.stem_bp_iterator(stem):
            nt1, nt2 = cg.seq[p1], cg.seq[p2]
            pair = {nt1, nt2}
            is_matched = any(pair == m for m in matched_pairs)
            if is_matched:
                basepairs.append([coords[p1 -1], coords[p2 -1]])
                basepairs_bases.append(pair)
            else:
                mismatched_basepairs.append([coords[p1 -1], coords[p2 -1]])

    if basepairs:
        basepairs = np.array(basepairs)
        c = "red" if color else "black"
        bpkwargs = {"color": c, "zorder": 0, "linewidth": 3}
        bpkwargs.update(basepair_kwargs)
        ax.plot(basepairs[:, :, 0].T, basepairs[:, :, 1].T, **bpkwargs)

    if mismatched_basepairs:
        mismatched_basepairs = np.array(mismatched_basepairs)
        mismatch_kwargs.update({"color":"red", "zorder": 0, "linewidth": 1, "linestyle": "--"})
        ax.plot(mismatched_basepairs[:, :, 0].T, mismatched_basepairs[:, :, 1].T, **mismatch_kwargs)

    # Nucleotides circles and labels
    for i, coord in enumerate(coords):
        if color:
            c = el_to_color.get(el_string[i], 'black')
            h, l, s = colorsys.rgb_to_hls(*mc.to_rgb(c))
            if lighten > 0.0:
                l += (1 - l) * min(1, lighten)
            else:
                l += l * max(-1, lighten)
            c = colorsys.hls_to_rgb(h, l, s)
            circle = plt.Circle(coord, color=c)
        else:
            circle = plt.Circle(coord, edgecolor="black", facecolor="white")
        ax.add_artist(circle)
        if cg.seq and i+1 in cg.seq:
            if "fontweight" not in text_kwargs:
                text_kwargs["fontweight"] = "bold"
            ax.annotate(cg.seq[i + 1], xy=coord, ha="center", va="center", **text_kwargs)

    all_coords = list(coords)
    ntnum_kwargs = {"color": "gray"}
    ntnum_kwargs.update(text_kwargs)
    for nt in range(10, cg.seq_length, 10):
        annot_pos = _find_annot_pos_on_circle(nt, all_coords, cg)
        if annot_pos is not None:
            ax.annotate(str(nt), xy=coords[nt - 1], xytext=annot_pos,
                        arrowprops={"width": 1, "headwidth": 1, "color": "gray"},
                        ha="center", va="center", zorder=0, **ntnum_kwargs)
            all_coords.append(annot_pos)

    _annotate_rna_plot(ax, cg, all_coords, annotations, text_kwargs)
    xmin, xmax = min(coords[:, 0]), max(coords[:, 0])
    ymin, ymax = min(coords[:, 1]), max(coords[:, 1])
    datalim = ((xmin, ymin), (xmax, ymax))
    width, height = xmax - xmin, ymax - ymin
    ax.set_aspect('equal', 'datalim')
    ax.update_datalim(datalim)
    ax.autoscale_view()
    ax.set_axis_off()
    return (ax, coords)


def _find_annot_pos_on_circle(nt, coords, cg):
    for i in range(5):
        for sign in [-1, 1]:
            a = np.pi / 4 * i * sign
            if cg.get_elem(nt)[0] == "s":
                bp = cg.pairing_partner(nt)
                anchor = coords[bp - 1]
            else:
                anchor = np.mean([coords[nt - 2], coords[nt]], axis=0)
            vec = coords[nt - 1] - anchor
            vec = vec / ftuv.magnitude(vec)
            rotated_vec = np.array([
                vec[0] * math.cos(a) - vec[1] * math.sin(a),
                vec[0] * math.sin(a) + vec[1] * math.cos(a)
            ])
            annot_pos = coords[nt - 1] + rotated_vec * 18
            if _clashfree_annot_pos(annot_pos, coords):
                return annot_pos
    return None


def save_rna_svg(rna_structure_file, sequence, rna_id, outdir):
    fig, ax = plt.subplots(figsize=(20, 20))
    cg = forgi.load_rna(rna_structure_file, allow_many=False)
    ax, coords = plot_rna(
        cg,
        ax=ax,
        text_kwargs={"fontsize": 12, "fontweight": "black"},
        lighten=0.5,
        backbone_kwargs={"color": 'blue', "linewidth": 1},
        basepair_kwargs={"linewidth": 3, "color": 'blue'},
        mismatch_kwargs={"color": "red", "zorder": 0, "linewidth": 1, "linestyle": "--"},
        color=False,
        annotations={}
    )
    for i, (base, coord) in enumerate(zip(sequence, coords)):
        x, y = coord
        # Optional: color bases by position; here simplified as black
        ax.text(x, y, base, color='black', fontsize=12, ha='center', va='center', fontweight='bold')

    plt.title(f"RNA Secondary Structure of {rna_id}")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    svg_file = os.path.join(outdir, f"{rna_id}.svg")
    plt.savefig(svg_file, format="svg")
    plt.close(fig)
    log.info(f"SVG saved to {svg_file}")


def main():
    parser = argparse.ArgumentParser(description="Visualize RNA secondary structure as SVG.")
    parser.add_argument('-s', '--structure_file', required=True,
                        help='RNA secondary structure file (.fx format)')
    parser.add_argument('-q', '--sequence', required=True,
                        help='RNA sequence corresponding to the structure')
    parser.add_argument('-i', '--id', required=True,
                        help='Identifier for RNA (used for output file naming)')
    parser.add_argument('-o', '--outdir', required=True,
                        help='Output directory to save SVG file')
    args = parser.parse_args()

    save_rna_svg(args.structure_file, args.sequence, args.id, args.outdir)


if __name__ == "__main__":
    main()
