#                                                                    #
#           PLOT ARROWS FOR GENE CLUSTER GIVEN A GenBank FILE        #
#            https://github.com/petercim/Arrower/tree/master         #
#                           Peter Cimermancic                        #
#                               April 2010                           #
#                heavily modified by He Beibei 2024                  #

import random, ast, cairosvg, os
import concurrent.futures
from Bio import SeqIO
import sys
import pandas as pd
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

PFAM_LIST = [
    "PF01867", "PF01507", "PF09704", "PF08798", "PF09707", "PF09344", "PF03460", "PF00009", "PF03144", "PF04389", "PF01583",
    "PF12084", "PF03819", "PF04055", "PF01077", "PF09827", "PF00113", "PF03952", "PF07690", "PF00117", "PF13395", "PF06418",
    "PF07722", "PF00005", "PF00370", "PF00667", "PF00175", "PF00258", "PF01242", "PF18541", "PF16593", "PF00753", "PF04536",
    "PF02782", "PF13561", "PF05107", "PF18019", "PF01565", "PF02913", "PF09614", "PF02597", "PF09615", "PF12832", "PF00083",
    "PF07731", "PF07732", "PF00994", "PF00394", "PF03606", "PF03453", "PF02391", "PF09611", "PF12799", "PF09709", "PF00664",
    "PF00072", "PF16592", "PF00278", "PF13457", "PF02784", "PF21384", "PF09711", "PF16595", "PF09618", "PF13304", "PF01381",
    "PF18470", "PF00270", "PF03787", "PF00271", "PF00501", "PF02518", "PF13193", "PF00589", "PF00165", "PF04011", "PF21802",
    "PF00528", "PF09485", "PF00512", "PF13411", "PF00004", "PF01469", "PF00990", "PF09278", "PF00419", "PF00550", "PF09481",
    "PF00746", "PF17210", "PF00400", "PF08447", "PF01391", "PF00668", "PF12833", "PF04122", "PF04270", "PF00248"
]

des = pd.read_csv('/home/hebeibei/Data/pfam/pfam_acc2des.txt',sep='\t')
des_dict = des.set_index(des.columns[0]).to_dict()[des.columns[1]]


def arrow(X, Y, L, H, strand, h, l, color):
    '''
    SVG code for arrow with improved visualization:
        - (X,Y) ... upper left (+) or right (-) corner of the arrow
        - L ... arrow length
        - H ... arrow height
        - strand
        - h ... arrow head edge width
        - l ... arrow head length
        - color ... color of the arrow
    '''
    if strand == '+':
        A = [X, Y]
        B = [X + L - l, Y]
        C = [X + L - l, Y - h]
        D = [X + L, Y + H / 2]
        E = [X + L - l, Y + H + h]
        F = [X + L - l, Y + H]
        G = [X, Y + H]

        if L < l:
            B = [X, Y]
            C = [X, Y - h]
            D = [X + L, Y + H / 2]
            E = [X, Y + H + h]
            F = [X, Y + H]

    elif strand == '-':
        A = [X + L, Y]
        B = [X + l, Y]
        C = [X + l, Y - h]
        D = [X, Y + H / 2]
        E = [X + l, Y + H + h]
        F = [X + l, Y + H]
        G = [X + L, Y + H]

        if L < l:
            B = [X + L, Y]
            C = [X + L, Y - h]
            D = [X, Y + H / 2]
            E = [X + L, Y + H + h]
            F = [X + L, Y + H]

    else:
        return 0

    line = f'''<polygon points="{A[0]},{A[1]} {B[0]},{B[1]} {C[0]},{C[1]} {D[0]},{D[1]} {E[0]},{E[1]} {F[0]},{F[1]} {G[0]},{G[1]}"
               style="fill:{color};fill-opacity:1;
               stroke:#000000;stroke-width:0">
               </polygon>'''
    return line


def color_from_locus_tag(locus_tag):
    '''
    Generate a consistent color based on the locus_tag.
    '''
    random.seed(locus_tag)
    return "#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])


def color_from_pfam(pfam):
    '''
    Generate a consistent color based on the pfam.
    '''
    if pfam == "":
        return "#cccccc"  # Gray color for empty pfam
    random.seed(pfam)
    return "#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])


def random_color():
    '''
    Generate a random color in HTML hex format.
    '''
    return "#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])


def repeat_line(X, Y, height):
    '''
    Draw a vertical line representing a repeat region.
    '''
    return f'<line x1="{X}" y1="{Y}" x2="{X}" y2="{Y + height}" style="stroke:#FF0000;stroke-width:0.5"/>'


def draw_connector(x1, y1, x2, y2):
    '''
    Draw a connector line between two points with a right-angle bracket shape.
    '''
    mid_y = (y1 + y2) / 2 + 30
    y1 += 3
    connector = f'<path d="M {x1},{y1} L {x1},{mid_y} L {x2},{mid_y} L {x2},{y1}" style="stroke:black; stroke-width:0.5; fill:none;"/>'
    text_x = (x1 + x2) / 2 - 60
    text_y = mid_y + 15
    text = f'<text x="{text_x}" y="{text_y}" style="font-family: Arial; font-size:10px;">Co-conserved array-protein</text>'
    return connector + text


def convert_svg_to_pdf(input_svg, output_pdf):
    cairosvg.svg2pdf(url=input_svg, write_to=output_pdf)
    print(f"PDF saved to {output_pdf}")


def save_svgs_to_pdf(svg_files, output_pdf):
    """
    Combine multiple SVG files into a single PDF, preserving text.
    """
    c = canvas.Canvas(output_pdf, pagesize=letter)
    width, height = letter

    for svg_file in svg_files:
        drawing = svg2rlg(svg_file)
        renderPDF.draw(drawing, c, 0, 0)
        c.showPage()  # Create a new page for each SVG file

    c.save()
    print(f"Saved to {output_pdf}")


def SVG(GenBankFile,output_file, ArrowHeight=20, HeadEdge=8, HeadLength=10, marginX=30, marginY=100, scaling=30.0, font=14,num_column=3):
    '''
    Create the main SVG document with enhanced visualization:
        - read in GenBank document
        - find genes, start and stop positions, and strands
        - write the SVG files with different colors for each gene
    '''
    GenBankFile_name = GenBankFile.split('/')[-1].replace('tracrRNA_False__','').replace('tracrRNA_True__','')
    ALL_TEXT = f"""<?xml version="1.0" standalone="no"?>
    <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" 
    "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">

    <svg width="210mm" height="297mm" xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'>
    <text x="{marginX}" y="{marginY - 10}" style="font-family: Arial; font-size:10px; font-weight: bold;">{GenBankFile_name}</text>
    """
    dr_added = False
    gene_label = 'A'
    gene_pfam_mapping = {}
    gene_highlight = {}
    array_position = None
    gene_positions = {}

    with open(GenBankFile, 'r') as file:
        for seq_record in SeqIO.parse(file, "genbank"):
            ClusterSize = len(seq_record.seq)
            scaling = ClusterSize * 30 / 20456
            ALL_TEXT += f'<line x1="{marginX}" y1="{marginY + ArrowHeight / 2}" x2="{marginX + ClusterSize / scaling}" y2="{marginY + ArrowHeight / 2}" style="stroke:rgb(99,99,99);stroke-width:1"/>'

            # Add scale for the length of the gene cluster
            scale_x = marginX + len(GenBankFile_name) * 0.5 * 10 + 280  # len(text) * avg_char_width * font_size
            scale_y = marginY - 20
            x1 = scale_x
            x2 = scale_x + 1000 / scaling
            ALL_TEXT += f'<line x1="{x1}" y1="{scale_y + 5}" x2="{x2}" y2="{scale_y + 5}" style="stroke:#808080;stroke-width:1"/>'
            ALL_TEXT += f'<text x="{(x1 + x2) / 2 - 8 }" y="{scale_y}" style="font-family: Arial; font-size:10px;">1 kb</text>'
            dot_radius = 2  # You can adjust the size of the dot
            ALL_TEXT += f'<circle cx="{x1}" cy="{scale_y + 5}" r="{dot_radius}" style="fill:#808080;"/>'
            ALL_TEXT += f'<circle cx="{x2}" cy="{scale_y + 5}" r="{dot_radius}" style="fill:#808080;"/>'
            previousStop = 0
            previousLevel = 1
            for feature in seq_record.features:
                if feature.type == 'CDS':
                    GeneName = feature.qualifiers.get('gene', [''])[0]
                    # locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
                    pfam = feature.qualifiers.get('pfam', [''])[0]
                    product = feature.qualifiers.get('product', [''])[0]
                    strand = '+' if feature.location == -1 else '-'
                    start = int(str(feature.location.start).replace('>', '').replace('<', '')) / scaling
                    stop = int(str(feature.location.end).replace('>', '').replace('<', '')) / scaling

                    # color = random_color()
                    # color = color_from_locus_tag(locus_tag)
                    color = color_from_pfam(pfam)
                    ALL_TEXT += arrow(start + marginX, marginY, stop - start, ArrowHeight, strand, HeadEdge, HeadLength, color)

                    if previousLevel == 4: previousLevel = 0
                    text_y_offset = marginY + 50 + (font + 6) * previousLevel
                    text_x = marginX + start + (stop - start) / 2 - len(GeneName) * font / 4
                    ALL_TEXT += f'<text x="{text_x}" y="{text_y_offset}" style="font-family: Arial; font-size:{font}; font-style: italic;">{GeneName}</text>'

                    # Add pfam number below the arrow
                    # if pfam:
                    gene_pfam_mapping[f'Gene {gene_label}'] = pfam if pfam else 'No Pfam domain'
                    gene_highlight[f'Gene {gene_label}'] = product
                    gene_positions[gene_label] = (marginX + start + (stop - start) / 2, marginY + ArrowHeight / 2)
                    gene_text_y_offset = marginY + ArrowHeight + 15
                    text_x = marginX + start + (stop - start) / 2 - len(gene_label) * font / 4
                    ALL_TEXT += f'<text x="{text_x}" y="{gene_text_y_offset}" style="font-family: Arial; font-size:{font-4};font-weight: bold;">{gene_label}</text>'
                    gene_label = chr(ord(gene_label) + 1)  # Increment the label to the next letter

                    previousStop = stop
                    previousLevel += 1
                elif feature.type == 'repeat_region':
                    start = int(str(feature.location.start).replace('>', '').replace('<', '')) / scaling
                    stop = int(str(feature.location.end).replace('>', '').replace('<', '')) / scaling
                    repeat_len = stop - start
                    ALL_TEXT += repeat_line(start + marginX, marginY, ArrowHeight)
                    # if not dr_added:
                    #     ALL_TEXT += f'<text x="{start + marginX - 5}" y="{marginY + ArrowHeight + 15}" style="font-family: Arial; font-size:10px;">dr</text>'
                    #     dr_added = True

                    product = feature.qualifiers.get('product', [''])[0]
                    if product == 'entire_repeat':
                        start = int(str(feature.location.start).replace('>', '').replace('<', '')) / scaling
                        stop = int(str(feature.location.end).replace('>', '').replace('<', '')) / scaling
                        text_x = marginX + start + (stop - start) / 2 - 5 * font / 4
                        ALL_TEXT += f'<text x="{text_x}" y="{marginY + ArrowHeight + 15}" style="font-family: Arial; font-size:10px; font-weight: bold;">Array</text>'
                        array_position = (marginX + start + (stop - start) / 2, marginY + ArrowHeight + 15)
    if array_position:
        for label, pfam in gene_pfam_mapping.items():
            if gene_highlight[label] == 'Putative Cas':
                gene_x, gene_y = gene_positions[label[-1]]
                array_x, array_y = array_position
                ALL_TEXT += draw_connector(array_x, array_y, gene_x, gene_y)

    mapping_x = marginX
    mapping_y = marginY + 100
    line_height = font + 5

    # # defined by number of items per column
    # num_items_per_column = 4
    # column_height = 0
    # mapping_y_gene_label = marginY + 100
    # for idx, (k, v) in enumerate(gene_pfam_mapping.items()):
    #     condition = True if gene_highlight[k] == 'Putative Cas' else False
    #     v = ast.literal_eval(v)
    #     annotations = [f'{pfam_id}: {des_dict[pfam_id]}' for pfam_id in v]
    #     column_offset = (idx // num_items_per_column) * 250  # Offset to switch columns
    #     current_mapping_x = mapping_x + column_offset
    #     current_mapping_y = mapping_y_gene_label + (idx % num_items_per_column) * line_height
    #
    #     ALL_TEXT += f'<text x="{current_mapping_x + 10}" y="{current_mapping_y}" style="font-family: Arial; font-size:{font - 6}; font-weight: bold;">{k + " (Protein con-conserved with array)" if condition else k}</text>'
    #     for annotation in annotations:
    #         current_mapping_y += (font - 4)
    #         is_normal_weight = True if gene_highlight[k] == 'Putative Cas' else False
    #         ALL_TEXT += f'<text x="{current_mapping_x + 10}" y="{current_mapping_y}" style="font-family: Arial; font-size:{font - 6}; {"font-weight: bold;" if is_normal_weight else ""}">{annotation}</text>'
    #     mapping_y_gene_label += font
    #     if idx % num_items_per_column == 0:
    #         column_height = max(column_height, mapping_y)
    #     if (idx + 1) % num_items_per_column == 0:
    #         mapping_y_gene_label = column_height  # + 10  # Reset y-coordinate for the next column

    # defined by number of columns
    total_items = len(gene_pfam_mapping)
    column_split = total_items // num_column + (1 if total_items % num_column != 0 else 0)  # Ensure enough rows for all items
    column_height = 0
    mapping_y_gene_label = marginY + 100
    for idx, (k, v) in enumerate(gene_pfam_mapping.items()):
        condition = True if gene_highlight[k] == 'Putative Cas' else False
        v = ast.literal_eval(v) if v != 'No Pfam domain' else v
        annotations = [f'{pfam_id}: {des_dict[pfam_id]}' for pfam_id in v] if v != 'No Pfam domain' else ['No Pfam domain']

        # Calculate the offset for the columns
        column_offset = (idx // column_split) * 250  # Offset to switch columns
        current_mapping_x = mapping_x + column_offset
        current_mapping_y = mapping_y_gene_label + (idx % column_split) * line_height

        ALL_TEXT += f'<text x="{current_mapping_x + 10}" y="{current_mapping_y}" style="font-family: Arial; font-size:{font - 6}; font-weight: bold;">{k + " (Protein co-conserved with array)" if condition else k}</text>'

        # Add annotations below the gene label
        for annotation in annotations:
            current_mapping_y += (font - 4)
            is_normal_weight = True if gene_highlight[k] == 'Putative Cas' else False
            ALL_TEXT += f'<text x="{current_mapping_x + 10}" y="{current_mapping_y}" style="font-family: Arial; font-size:{font - 6}; {"font-weight: bold;" if is_normal_weight else ""}">{annotation}</text>'

        mapping_y_gene_label += max([10 * len(annotations),font])

        # Update column height and position if needed
        if idx % column_split == 0:
            column_height = max(column_height, mapping_y)
        if (idx + 1) % column_split == 0:
            mapping_y_gene_label = column_height  # Reset y-coordinate for the next column

    ALL_TEXT += '</svg>'

    with open(output_file, 'w') as f:
        f.write(ALL_TEXT)
    print(f"SVG saved to {output_file}")

    # convert_svg_to_pdf(output_file, output_file.replace('.svg', '.pdf'))
    # print(ALL_TEXT)


if __name__ == "__main__":
    gbk_dir = '/Your/Path/To/Candidate'
    gbk_dirs = [os.path.join(gbk_dir, f) for f in os.listdir(gbk_dir) if f.endswith('.gb')]
    for gbk_file in gbk_dirs:
        gbk2svg = f'{gbk_dir}/{gbk_file.split('/')[-1]}.svg'
        print(gbk_file)
        print(gbk2svg)
        SVG(gbk_file, gbk2svg)


    # svgs = ['./gbk_svg_5.svg','./gbk_svg_4.svg','./gbk_svg_3.svg']
    # output_pdf = 'gbk_merged.pdf'
    # save_svgs_to_pdf(svgs,output_pdf)


