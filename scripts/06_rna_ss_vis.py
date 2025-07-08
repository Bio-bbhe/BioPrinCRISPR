# modified from https://github.com/ViennaRNA/forgi/blob/master/forgi/visual/mplotlib.py
from __future__ import division
from __future__ import print_function
from builtins import zip
import forgi.threedee.utilities.vector as ftuv
import math, logging, itertools, colorsys, forgi
import matplotlib.pyplot as plt
import numpy as np

log = logging.getLogger(__name__)


def _clashfree_annot_pos(pos, coords):
    for c in coords:
        dist = ftuv.vec_distance(c, pos)
        # log.debug("vec_dist=%s", dist)
        if dist < 14:
            return False
    return True


def _annotate_rna_plot(ax, cg, coords, annotations, text_kwargs):
    # Plot annotations
    annot_dict = {elem: elem for elem in cg.defines}
    if annotations is None:
        annot_dict = {elem: "" for elem in cg.defines}
    else:
        annot_dict.update(annotations)
    stem_coords = {}
    for stem in cg.stem_iterator():
        stem_start = np.mean([coords[cg.defines[stem][0] - 1],
                              coords[cg.defines[stem][3] - 1]],
                             axis=0)
        stem_end = np.mean([coords[cg.defines[stem][1] - 1],
                            coords[cg.defines[stem][2] - 1]],
                           axis=0)
        stem_center = np.mean([stem_start, stem_end], axis=0)
        stem_coords[stem] = (stem_start, stem_center, stem_end)
        if annot_dict[stem]:
            stem_vec = stem_end - stem_start
            norm_vec = (stem_vec[1], -stem_vec[0])
            norm_vec /= ftuv.magnitude(norm_vec)
            annot_pos = np.array(stem_center) + 23 * norm_vec
            # log.debug("Checking clashfree for %s, %s", stem, annot_pos)
            if not _clashfree_annot_pos(annot_pos, coords):
                log.debug(
                    "Cannot annotate %s as %s ON THE RIGHT HAND SIDE, because of insufficient space. Trying left side...",
                    stem, annot_dict[stem])
                annot_pos = np.array(stem_center) - 23 * norm_vec
                # log.debug("Checking clashfree OTHER SIDE for %s, %s", stem, annot_pos)
                if not _clashfree_annot_pos(annot_pos, coords):
                    log.info("Cannot annotate %s as '%s', because of insufficient space.", stem, annot_dict[stem])
                    annot_pos = None
            # log.debug("%s", annot_pos)
            if annot_pos is not None:
                ax.annotate(annot_dict[stem], xy=annot_pos,
                            ha="center", va="center", **text_kwargs)
    for hloop in cg.hloop_iterator():
        hc = []
        for nt in cg.define_residue_num_iterator(hloop, adjacent=True):
            hc.append(coords[nt - 1])
        annot_pos = np.mean(hc, axis=0)
        if _clashfree_annot_pos(annot_pos, coords):
            ax.annotate(annot_dict[hloop], xy=annot_pos,
                        ha="center", va="center", **text_kwargs)
        else:
            log.info("Cannot annotate %s as '%s' ON THE INSIDE, because of insufficient space. Trying outside...",
                     hloop, annot_dict[hloop])
            nt1, nt2 = cg.define_a(hloop)
            start = np.mean([coords[nt1 - 1], coords[nt2 - 1]], axis=0)
            vec = annot_pos - start
            annot_pos = annot_pos + vec * 3
            if _clashfree_annot_pos(annot_pos, coords):
                ax.annotate(annot_dict[hloop], xy=annot_pos,
                            ha="center", va="center", **text_kwargs)
            else:
                log.info("Cannot annotate %s as '%s', because of insufficient space.", hloop, annot_dict[hloop])
    for iloop in cg.iloop_iterator():
        s1, s2 = cg.connections(iloop)
        annot_pos = np.mean([stem_coords[s1][2], stem_coords[s2][0]], axis=0)
        if _clashfree_annot_pos(annot_pos, coords):
            ax.annotate(annot_dict[iloop], xy=annot_pos,
                        ha="center", va="center", **text_kwargs)
        else:
            log.debug("Cannot annotate %s as '%s' ON THE INSIDE, because of insufficient space. Trying outside...",
                      iloop, annot_dict[iloop])
            loop_vec = stem_coords[s2][0] - stem_coords[s1][2]
            norm_vec = (loop_vec[1], -loop_vec[0])
            norm_vec /= ftuv.magnitude(norm_vec)
            annot_pos_p = np.array(annot_pos) + 25 * norm_vec
            annot_pos_m = np.array(annot_pos) - 25 * norm_vec
            # iloops can be asymmetric (more nts on one strand.)
            # plot the label on the strand with more nts.
            plus = 0
            minus = 0
            for nt in cg.define_residue_num_iterator(iloop):
                if ftuv.vec_distance(annot_pos_p, coords[nt - 1]) < ftuv.vec_distance(annot_pos_m, coords[nt - 1]):
                    plus += 1
                else:
                    minus += 1
            if plus > minus:
                if _clashfree_annot_pos(annot_pos_p, coords):
                    ax.annotate(annot_dict[iloop], xy=annot_pos_p,
                                ha="center", va="center", **text_kwargs)
                else:
                    log.info(
                        "Cannot annotate %s as '%s' (only trying inside and right side), because of insufficient space.",
                        iloop, annot_dict[iloop])

            else:
                if _clashfree_annot_pos(annot_pos_m, coords):
                    ax.annotate(annot_dict[iloop], xy=annot_pos_m,
                                ha="center", va="center", **text_kwargs)
                else:
                    log.info(
                        "Cannot annotate %s as '%s' (only trying inside and left side), because of insufficient space.",
                        iloop, annot_dict[iloop])
    for mloop in itertools.chain(cg.floop_iterator(), cg.tloop_iterator(), cg.mloop_iterator()):
        nt1, nt2 = cg.define_a(mloop)
        res = list(cg.define_residue_num_iterator(mloop))
        if len(res) == 0:
            anchor = np.mean([coords[nt1 - 1], coords[nt2 - 1]], axis=0)
        elif len(res) % 2 == 1:
            anchor = coords[res[int(len(res) // 2)] - 1]
        else:
            anchor = np.mean([coords[res[int(len(res) // 2) - 1] - 1],
                              coords[res[int(len(res) // 2)] - 1]],
                             axis=0)
        loop_vec = coords[nt1 - 1] - coords[nt2 - 1]
        norm_vec = (loop_vec[1], -loop_vec[0])
        norm_vec /= ftuv.magnitude(norm_vec)
        annot_pos = anchor - norm_vec * 18
        if _clashfree_annot_pos(annot_pos, coords):
            ax.annotate(annot_dict[mloop], xy=annot_pos,
                        ha="center", va="center", **text_kwargs)
        else:
            log.info("Cannot annotate %s as '%s' , because of insufficient space.", mloop, annot_dict[mloop])


def plot_rna(cg, ax=None, offset=(0, 0), color=True, text_kwargs=None, backbone_kwargs=None,
             basepair_kwargs=None, mismatch_kwargs=None, lighten=0.0, annotations=None):
    '''
    Plot an RNA structure given a set of nucleotide coordinates

    .. note::

        This function calls set_axis_off on the axis. You can revert this by
        using ax.set_axis_on() if you like to see the axis.

    :param color: color to draw paired bases
    :param turning_points: List of indices in the sequence where the structure should "turn" to form a branched shape.
                           If None, turning points are determined automatically based on RNA structure features.
    :param cg: A forgi.threedee.model.coarse_grain.CoarseGrainRNA structure
    :param ax: A matplotlib plotting area
    :param offset: Offset the plot by these coordinates. If a simple True is passed in, then
                   offset by the current width of the plot
    :param text_kwargs: keyword arguments passed to matplotlib.pyplot.annotate
                        for plotting of the sequence
    :param backbone_kwargs: keyword arguments passed to matplotlib.pyplot.plot
                        for plotting of the backbone links
    :param basepair_kwargs: keyword arguments passed to matplotlib.pyplot.plot
                        for plotting of the basepair links
    :param mismatch_kwargs: for mismatched base pairs  # hbb added  2024-09-10
    :param lighten: Make circles lighter. A percent value where 1 makes
                    everything white and 0 leaves the colors unchanged
    :param annotations: A dictionary {elem_name: string} or None.
                        By default, the element names (e.g. "s0") are plotted
                        next to the element. This dictionary can be used to
                        override the default element names by costum strings.
                        To remove individual annotations, assign an empty string to the key.
                        To remove all annotations, set this to None.

                        .. warning::

                            Annotations are not shown, if there is not enough space.
                            Annotations not shown are logged with level INFO
    :return: (ax, coords) The axes and the coordinates for each nucleotide
    '''
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
    log.info("Starting to plot RNA...")
    import RNA
    import matplotlib.colors as mc
    RNA.cvar.rna_plot_type = 1

    coords = []
    # colors = []
    # circles = []

    bp_string = cg.to_dotbracket_string()
    # get the type of element of each nucleotide
    el_string = cg.to_element_string()
    # i.e. eeesssshhhhsssseeee
    el_to_color = {'f': 'orange',
                   't': 'orange',
                   's': 'green',
                   'h': 'blue',
                   'i': 'yellow',
                   'm': 'red'}

    if ax is None:
        ax = plt.gca()

    if offset is None:
        offset = (0, 0)
    elif offset is True:
        offset = (ax.get_xlim()[1], ax.get_ylim()[1])
    else:
        pass

    vrna_coords = RNA.get_xy_coordinates(bp_string)
    # TODO Add option to rotate the plot
    for i, _ in enumerate(bp_string):
        coord = (offset[0] + vrna_coords.get(i).X,
                 offset[1] + vrna_coords.get(i).Y)
        coords.append(coord)
    coords = np.array(coords)

    # First plot backbone
    bkwargs = {"color": "black", "zorder": 0}
    bkwargs.update(backbone_kwargs)
    ax.plot(coords[:, 0], coords[:, 1], **bkwargs)
    # Now plot basepairs
    basepairs = []
    basepairs_bases = []  # hbb added,2024-09-10, for mismatched base pairs
    mismatched_basepairs = []  # hbb added,2024-09-10, for mismatched base pairs
    for s in cg.stem_iterator():
        for p1, p2 in cg.stem_bp_iterator(s):
            # --------------- hbb added,2024-09-10, for mismatched base pairs: --------------- #
            matched_pair = [{'A', 'u'}, {'G', 'C'}, {'U', 'A'}, {'C', 'G'}]
            # print(p1, p2)
            nt1 = cg.seq[p1]  # Access nucleotide at position p1
            nt2 = cg.seq[p2]  # Access nucleotide at position p2
            pair = {nt1, nt2}
            matched = [pair.intersection(p) for p in matched_pair if len(pair.intersection(p)) == 2]
            if not matched:
                # Append mismatched base pair coordinates to the list
                mismatched_basepairs.append([coords[p1 - 1], coords[p2 - 1]])
            else:
                basepairs.append([coords[p1 - 1], coords[p2 - 1]])
                basepairs_bases.append(pair)  # hbb added,2024-09-10, for mismatched base pairs

    if basepairs:
        basepairs = np.array(basepairs)
        # print(f'basepairs is {basepairs}')
        # print(f'basepairs_bases i {basepairs_bases}')
        if color:
            c = "red"
        else:
            c = "black"
            bpkwargs = {"color": c, "zorder": 0, "linewidth": 3}
            bpkwargs.update(basepair_kwargs)
            ax.plot(basepairs[:, :, 0].T, basepairs[:, :, 1].T, **bpkwargs)

    # --------------- hbb added,2024-09-10, for mismatched base pairs: --------------- #
    if mismatched_basepairs:
        mismatched_basepairs = np.array(mismatched_basepairs)
        mismatch_dict = {"color": "red", "zorder": 0, "linewidth": 1, "linestyle": "--"}
        mismatch_dict.update(mismatch_kwargs)
        ax.plot(mismatched_basepairs[:, :, 0].T, mismatched_basepairs[:, :, 1].T, **mismatch_kwargs)

    # --------------- hbb added,2024-09-10, end of  editing --------------- #

    # Now plot circles
    for i, coord in enumerate(coords):
        if color:
            c = el_to_color[el_string[i]]
            h, l, s = colorsys.rgb_to_hls(*mc.to_rgb(c))
            if lighten > 0.0:
                l += (1 - l) * min(1, lighten)
            else:
                l += l * max(-1, lighten)
            if l > 1.0 or l < 0.0:
                print(l)
            c = colorsys.hls_to_rgb(h, l, s)
            circle = plt.Circle((coord[0], coord[1]),
                                color=c)
        else:
            circle = plt.Circle((coord[0], coord[1]),
                                edgecolor="black", facecolor="white")

        ax.add_artist(circle)
        if cg.seq:
            if "fontweight" not in text_kwargs:
                text_kwargs["fontweight"] = "bold"
            ax.annotate(cg.seq[i + 1], xy=coord, ha="center", va="center", **text_kwargs)

    all_coords = list(coords)
    ntnum_kwargs = {"color": "gray"}
    ntnum_kwargs.update(text_kwargs)
    for nt in range(10, cg.seq_length, 10):
        # We try different angles
        annot_pos = _find_annot_pos_on_circle(nt, all_coords, cg)
        if annot_pos is not None:
            ax.annotate(str(nt), xy=coords[nt - 1], xytext=annot_pos,
                        arrowprops={"width": 1, "headwidth": 1, "color": "gray"},
                        ha="center", va="center", zorder=0, **ntnum_kwargs)
            all_coords.append(annot_pos)

    _annotate_rna_plot(ax, cg, all_coords, annotations, text_kwargs)
    datalim = ((min(list(coords[:, 0]) + [ax.get_xlim()[0]]),
                min(list(coords[:, 1]) + [ax.get_ylim()[0]])),
               (max(list(coords[:, 0]) + [ax.get_xlim()[1]]),
                max(list(coords[:, 1]) + [ax.get_ylim()[1]])))

    '''
    min_coord = min(datalim[0][0], datalim[0][1])
    max_coord = max(datalim[1][0], datalim[1][1])
    datalim = ((min_coord, min_coord), (max_coord, max_coord))

    print "min_coord:", min_coord
    print "max_coord:", max_coord
    print "datalime:", datalim
    '''

    width = datalim[1][0] - datalim[0][0]
    height = datalim[1][1] - datalim[0][1]

    # ax.set_aspect(width / height)
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
            rotated_vec = np.array([vec[0] * math.cos(a) - vec[1] * math.sin(a),
                                    vec[0] * math.sin(a) + vec[1] * math.cos(a)])
            annot_pos = coords[nt - 1] + rotated_vec * 18
            if _clashfree_annot_pos(annot_pos, coords):
                log.debug("Annot pos on c is %s", annot_pos)
                return annot_pos
    return None


def get_base_color(position):
    if position < 10:
        return 'blue'  # Bases 1-9
    elif 10 <= position < 20:
        return 'green'  # Bases 10-19
    elif 20 <= position < 30:
        return 'orange'  # Bases 20-29
    else:
        return 'red'  # Bases 30 and above


def save_rna_svg(rna_fasta,sequence, rna_id, dir):
    fig, ax = plt.subplots(figsize=(20, 20))
    cg = forgi.load_rna(rna_fasta, allow_many=False)
    ax, coords = plot_rna(cg,
                          ax=ax,
                          text_kwargs={"fontweight": "black", "fontsize": 12},
                          lighten=0.5,
                          backbone_kwargs={"color": 'blue', "linewidth": 1},
                          basepair_kwargs={"linewidth": 3, "color": 'blue'},
                          mismatch_kwargs={"color": "red", "zorder": 0, "linewidth": 1, "linestyle": "--"},
                          color=False,
                          annotations={}
                          )
    for i, (base, coord) in enumerate(zip(sequence, coords)):
        x, y = coord  # Coordinates for each base
        color = get_base_color(i)  # Get color based on position
        ax.text(x, y, base, color=color, fontsize=12, ha='center', va='center', fontweight='bold')

    # Save the plot as an SVG file
    plt.title(f"RNA Secondary Structure of {rna_id}")
    plt.savefig(f"{dir}/{rna_id}.svg", format="svg")


if __name__ == '__main__':
    dot_bracket = '.((..(((((.((....)).)...)))).((((((((((((....))).))......)))))))......(((.(((((.(((((....))))).))))).))).....)).'
    sequence = 'GGCAACGGCGGCGGCAACGGCGGAGCCGGCGGUGCCGGGGGAACGCCCACCGGCAGUGGCACCGAGGGGACCGGCGGCGACGGUGGAGAUGCCGGCGCCGGCGGCAACGGCG'  # Example sequence with "LLL" as the linker

    rna_structure_str = f">seq\n{sequence}\n{dot_bracket}"

    # Since forgi.load_rna works with files, we need to use a temporary file to load the structure
    with open("rna_temp.fx", "w") as f:
        f.write(rna_structure_str)
    print(rna_structure_str)
    rna_fasta, rna_id, dir = "rna_temp.fx", "rna_temp", "/home/hebeibei/Work/crispr/code"
    save_rna_svg(rna_fasta, rna_id, dir)
