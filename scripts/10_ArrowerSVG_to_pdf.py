#                                                                    #
#           PLOT ARROWS FOR GENE CLUSTER GIVEN A GenBank FILE        #
#            https://github.com/petercim/Arrower/tree/master         #
#                           Peter Cimermancic                        #
#                               April 2010                           #
#                heavily modified by He Beibei 2024                  #

import random, ast, cairosvg, os, collections, shutil, re
from Bio import SeqIO
from PyPDF2 import PdfMerger
import sys
import pandas as pd

# first 100 pfam
PFAM_LIST = [
    "PF01867", "PF01507", "PF09704", "PF08798", "PF09707", "PF09344", "PF03460", "PF00009", "PF03144", "PF04389",
    "PF01583", "PF12084", "PF03819", "PF04055", "PF01077", "PF09827", "PF00113", "PF03952", "PF07690", "PF00117",
    "PF13395", "PF06418", "PF07722", "PF00005", "PF00370", "PF00667", "PF00175", "PF00258", "PF01242", "PF18541",
    "PF16593", "PF00753", "PF04536", "PF02782", "PF13561", "PF05107", "PF18019", "PF01565", "PF02913", "PF09614",
    "PF02597", "PF09615", "PF12832", "PF00083", "PF07731", "PF07732", "PF00994", "PF00394", "PF03606", "PF03453",
    "PF02391", "PF09611", "PF12799", "PF09709", "PF00664", "PF00072", "PF16592", "PF00278", "PF13457", "PF02784",
    "PF21384", "PF09711", "PF16595", "PF09618", "PF13304", "PF01381", "PF18470", "PF00270", "PF03787", "PF00271",
    "PF00501", "PF02518", "PF13193", "PF00589", "PF00165", "PF04011", "PF21802", "PF00528", "PF09485", "PF00512",
    "PF13411", "PF00004", "PF01469", "PF00990", "PF09278", "PF00419", "PF00550", "PF09481", "PF00746", "PF17210",
    "PF00400", "PF08447", "PF01391", "PF00668", "PF12833", "PF04122", "PF04270", "PF00248"
]

# pfam 101-200
pf_list_101_200 = [
    "PF00107", "PF07992", "PF00440", "PF00665", "PF08240", "PF13401", "PF00583", "PF13424",
    "PF01930", "PF09299", "PF13683", "PF00535", "PF06445", "PF02687", "PF13565", "PF13673",
    "PF17953", "PF00486", "PF13437", "PF09659", "PF13306", "PF13518", "PF02498", "PF00155",
    "PF00353", "PF00849", "PF20020", "PF02861", "PF01380", "PF01832", "PF13538", "PF13412",
    "PF00106", "PF00108", "PF00679", "PF06421", "PF00563", "PF01418", "PF00293", "PF14492",
    "PF03193", "PF07724", "PF16576", "PF10431", "PF00672", "PF04542", "PF00929", "PF18395",
    "PF03382", "PF00561", "PF07670", "PF02803", "PF00334", "PF02463", "PF17871", "PF00903",
    "PF19127", "PF00805", "PF02769", "PF00586", "PF13242", "PF01979", "PF03588", "PF12704",
    "PF01905", "PF07883", "PF06993", "PF13560", "PF12802", "PF13361", "PF11975", "PF02056",
    "PF13475", "PF00580", "PF00480", "PF08281", "PF00291", "PF07314", "PF10609", "PF03458",
    "PF01554", "PF13443", "PF00196", "PF00702", "PF13614", "PF11047", "PF02321", "PF00300",
    "PF18070", "PF01176", "PF00534", "PF19079", "PF02617", "PF08448", "PF00069", "PF10040",
    "PF05016", "PF00515", "PF13855", "PF14659"
]

pf_list_transposase = [
    "PF09299", "PF01609", "PF05598", "PF01610", "PF13340", "PF12784", "PF01548",
    "PF01527", "PF13586", "PF04754", "PF14690", "PF13751", "PF01710", "PF03400",
    "PF00872", "PF02371", "PF01797", "PF13701", "PF13612", "PF13006", "PF04986",
    "PF07282", "PF01526", "PF14319", "PF13542", "PF12762", "PF13022", "PF12760",
    "PF01498", "PF02914", "PF10551", "PF01385", "PF13005", "PF03221", "PF07592",
    "PF02281", "PF03050", "PF04693", "PF11427"
]

endonuclease = [
    "PF00665", "PF00929", "PF01042", "PF00565", "PF12705", "PF12784", "PF03372", "PF01934", "PF11867", "PF00730",
    "PF20469", "PF08411", "PF02342", "PF10576", "PF13613", "PF07510", "PF02870", "PF13359", "PF13546", "PF01541",
    "PF19580", "PF14622", "PF04383", "PF08797", "PF18451", "PF05685", "PF18813", "PF10298", "PF02646", "PF04471",
    "PF13930", "PF14279", "PF13366", "PF13091", "PF03432", "PF13358", "PF18335", "PF13558", "PF13166", "PF14529",
    "PF00614", "PF08378", "PF04493", "PF14436", "PF18555", "PF04231", "PF14040", "PF01939", "PF14412", "PF13391",
    "PF14414", "PF08721", "PF02577", "PF14588", "PF01935", "PF14338", "PF09565", "PF01612", "PF13643", "PF16902",
    "PF09810", "PF18819", "PF02732", "PF20600", "PF05901", "PF00445", "PF16786", "PF15518", "PF02075", "PF15648",
    "PF13392", "PF00825", "PF08696", "PF08774", "PF05866", "PF08011", "PF14890", "PF19215", "PF14281", "PF20704",
    "PF06958", "PF12639", "PF09545", "PF14082", "PF09568", "PF16473", "PF06023", "PF18729", "PF08929", "PF05066",
    "PF10108", "PF18497", "PF18731", "PF18739", "PF05203", "PF06479", "PF08722", "PF06877", "PF01713", "PF14411",
    "PF07000", "PF08928", "PF20731", "PF09553", "PF12183", "PF09062", "PF14410", "PF01949", "PF18738", "PF01076",
    "PF02720", "PF09225", "PF02601", "PF10150", "PF13396", "PF21386", "PF06044", "PF18334", "PF18026", "PF18062",
    "PF09562", "PF04857", "PF10107", "PF20796", "PF13017", "PF18810", "PF18740", "PF14448", "PF20473", "PF05872",
    "PF09019", "PF21210", "PF18305", "PF14130", "PF03104", "PF10996", "PF20466", "PF04599", "PF09281", "PF02778",
    "PF18165", "PF12106", "PF18136", "PF16413", "PF21315", "PF16953", "PF09378", "PF07102", "PF20472", "PF16784",
    "PF18812", "PF21598", "PF09517", "PF04308", "PF00636", "PF15650", "PF12008", "PF21693", "PF05076", "PF13156",
    "PF08461", "PF08579", "PF05204", "PF02609", "PF09126", "PF17146", "PF14511", "PF00961", "PF04556", "PF03265",
    "PF11977", "PF05991", "PF08567", "PF17726", "PF18618", "PF15649", "PF09494", "PF09491", "PF17411", "PF02923",
    "PF18737", "PF15565", "PF18631", "PF21434", "PF18156", "PF18157", "PF18616", "PF03161", "PF18081", "PF02265",
    "PF18614", "PF18809", "PF07328", "PF07788", "PF04555", "PF08745", "PF09376", "PF08463", "PF05876", "PF18501",
    "PF03354", "PF02130", "PF02945", "PF08573", "PF20441", "PF04411", "PF15515", "PF21169", "PF20947", "PF15517",
    "PF13032", "PF00604", "PF20467", "PF17427", "PF05367", "PF14562", "PF01868", "PF15636", "PF04381", "PF20987",
    "PF18561", "PF09573", "PF21368", "PF09570", "PF09567", "PF15637", "PF06250", "PF21202", "PF20748", "PF03851",
    "PF20976", "PF15647", "PF09217", "PF17728", "PF06616", "PF07924", "PF18814", "PF17770", "PF09552", "PF03159",
    "PF18780", "PF00545", "PF12117", "PF18479", "PF09462", "PF20952", "PF01348", "PF18154", "PF15653", "PF18732",
    "PF20951", "PF07460", "PF20464", "PF06319", "PF18867", "PF09665", "PF09292", "PF15657", "PF11463", "PF06064",
    "PF09554", "PF21003", "PF01620", "PF13691", "PF11634", "PF09566", "PF21111", "PF10141", "PF19778", "PF06630",
    "PF01138", "PF03725", "PF09504", "PF18516", "PF09499", "PF11718", "PF20454", "PF18779", "PF09412", "PF03184",
    "PF01974", "PF08857", "PF18350", "PF11663", "PF20950", "PF13651", "PF21123", "PF02017", "PF09631", "PF20376",
    "PF04152", "PF18733", "PF02980", "PF08724", "PF15652", "PF19500", "PF20948", "PF01986", "PF04466"
]

PFAM_LIST = PFAM_LIST.extend(pf_list_101_200)

des = pd.read_csv('/home/hebeibei/Data/pfam/pfam_acc2des.txt', sep='\t')
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


def SVG(GenBankFiles, output_file, meta_info=None,ArrowHeight=20, HeadEdge=8, HeadLength=10, marginX=30, marginY=100, scaling=30.0,
        font=14, num_column=3):
    '''
    Create the main SVG document with enhanced visualization:
        - read in GenBank document
        - find genes, start and stop positions, and strands
        - write the SVG files with different colors for each gene
    '''
    ALL_TEXT = f"""<?xml version="1.0" standalone="no"?>
    <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" 
    "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">

    <svg width="210mm" height="297mm" xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'>
    <text x="{marginX}" y="{marginY - 30}" style="font-family: Arial; font-size:14px; font-weight: bold;">{meta_info}</text>
    """
    for GenBankFile in GenBankFiles:
        GenBankFile_name = GenBankFile.split('/')[-1].replace('tracrRNA_False__', '').replace('tracrRNA_True__', '')
        repeat_num = GenBankFile_name.split('_repeats_')[0].split('__')[-1]
        repeat_num = int(repeat_num)
        ALL_TEXT += f'<text x="{marginX}" y="{marginY - 10}" style="font-family: Arial; font-size:10px; font-weight: bold;">{GenBankFile_name}</text>'
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
                ALL_TEXT += f'<text x="{(x1 + x2) / 2 - 8}" y="{scale_y}" style="font-family: Arial; font-size:10px;">1 kb</text>'
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
                        strand = '+' if feature.location.strand == -1 else '-'
                        start = int(str(feature.location.start).replace('>', '').replace('<', '')) / scaling
                        stop = int(str(feature.location.end).replace('>', '').replace('<', '')) / scaling

                        # color = random_color()
                        # color = color_from_locus_tag(locus_tag)
                        color = color_from_pfam(pfam)
                        ALL_TEXT += arrow(start + marginX, marginY, stop - start, ArrowHeight, strand, HeadEdge,
                                          HeadLength, color)

                        if previousLevel == 4: previousLevel = 0
                        text_y_offset = marginY + 50 + (font + 6) * previousLevel
                        text_x = marginX + start + (stop - start) / 2 - len(GeneName) * font / 4
                        ALL_TEXT += f'<text x="{text_x}" y="{text_y_offset}" style="font-family: Arial; font-size:{font}; font-style: italic;">{GeneName}</text>'

                        # Add pfam number below the arrow
                        # if pfam:
                        gene_pfam_mapping[f'Gene {gene_label}'] = pfam if pfam else 'No Pfam domain'
                        gene_highlight[f'Gene {gene_label}'] = product
                        gene_positions[gene_label] = (
                            marginX + start + (stop - start) / 2, marginY + ArrowHeight / 2)
                        gene_text_y_offset = marginY + ArrowHeight + 15
                        text_x = marginX + start + (stop - start) / 2 - len(gene_label) * font / 4
                        ALL_TEXT += f'<text x="{text_x}" y="{gene_text_y_offset}" style="font-family: Arial; font-size:{font - 4};font-weight: bold;">{gene_label}</text>'
                        gene_label = chr(ord(gene_label) + 1)  # Increment the label to the next letter

                        previousStop = stop
                        previousLevel += 1
                    elif feature.type == 'repeat_region':
                        # start = int(str(feature.location.start).replace('>', '').replace('<', '')) / scaling
                        # ALL_TEXT += repeat_line(start + marginX, marginY, ArrowHeight)
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

                            interval = (stop - start) / (repeat_num + 1)  # Calculate interval for repeat lines

                            # Draw the specified number of repeat lines
                            for i in range(repeat_num):
                                x_position = start + (i + 1) * interval  # Calculate the x-position of each repeat line
                                ALL_TEXT += repeat_line(x_position + marginX, marginY, ArrowHeight)
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
        column_split = total_items // num_column + (
            1 if total_items % num_column != 0 else 0)  # Ensure enough rows for all items
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

            mapping_y_gene_label += max([10 * len(annotations), font])

            # Update column height and position if needed
            if idx % column_split == 0:
                column_height = max(column_height, mapping_y)
            if (idx + 1) % column_split == 0:
                mapping_y_gene_label = column_height  # Reset y-coordinate for the next column

        marginY += 350
    ALL_TEXT += '</svg>'

    with open(output_file, 'w') as f:
        f.write(ALL_TEXT)
    print(f"SVG saved to {output_file}")

    convert_svg_to_pdf(output_file, output_file.replace('.svg', '.pdf'))
    # print(ALL_TEXT)


def sort_key(filename):
    match = re.match(r"(\d{3}_PF\d+)_Cluster_(\d+)\.pdf", filename)
    if match:
        prefix = match.group(1)
        cluster_number = int(match.group(2))  # Convert cluster number to integer for correct sorting
        return (prefix, cluster_number)
    return (filename, 0)


if __name__ == "__main__":
    # # ------------------ for testing ------------------ ##
    # gbks = [
    #     '/home/hebeibei/Work/crispr/code/cas_co_NCBI_v6_10kb/2_on_domian/part_1_before_11_12/001_PF01507_Pho_ad_sulfate_red/cluster_4_crispr_minimal_0_repeats/tracrRNA_True__DAAYBX010000006.1__162810_163265_#14__8_repeats_3482__gap.gb',
    #     '/home/hebeibei/Work/crispr/code/cas_co_NCBI_v6_10kb/2_on_domian/part_1_before_11_12/001_PF01507_Pho_ad_sulfate_red/cluster_3_crispr_minimal_0_repeats/tracrRNA_False__AARXCY010000166.1__1689_2082_#0__7_repeats_1310__gap.gb',
    #     '/home/hebeibei/Work/crispr/code/cas_co_NCBI_v6_10kb/2_on_domian/part_1_before_11_12/013_PF03952_Enolase_N-terminal/cluster_1_crispr_minimal_0_repeats/tracrRNA_False__AAEEFS010000006.1__203421_204669_#9__21_repeats_1104__gap.gb']
    #
    # try:
    #     SVG(gbks, 'gbk_svg_merged.svg')
    # except Exception as e:
    #     print(f'''Incorrect usage. Please, read manual again:\n{e}''')
    #     sys.exit(1)

    # gbk_dir = '/media/Data/hebeibei/minced_output/array_cds_pairs_ncbi_10kb/CRISPR_candidates/2_on_domain_customized/nuclease/'
    # destination_dir = '/media/Data/hebeibei/minced_output/array_cds_pairs_ncbi_10kb/CRISPR_candidates/2_on_domain_customized/nuclease_PDF_SI/'
    gbk_dir = "/home/hebeibei/Work/crispr/code/cas_co_NCBI_v6_10kb/3_no_pfam/mmseq_cluster/"
    destination_dir = "/home/hebeibei/Work/crispr/code/cas_co_NCBI_v6_10kb/3_no_pfam/mmseq_cluster/PDF_SI/"
    merged_pdf = f'{destination_dir}/merged_output.pdf'
    subdirs = [d for d in os.listdir(gbk_dir) if os.path.isdir(os.path.join(gbk_dir, d))]
    n = 0
    for subdir in subdirs:
        family_dir = os.path.join(gbk_dir, subdir)
        cluster_dirs = [d for d in os.listdir(family_dir) if os.path.isdir(os.path.join(family_dir, d)) and d.startswith('cluster')]
        if cluster_dirs:
            for cluster_dir in cluster_dirs:
                cluster2gbk = collections.defaultdict(list)
                for file in os.listdir(os.path.join(family_dir, cluster_dir)):
                    if file.endswith('.gb'):
                        cluster2gbk[cluster_dir].append(os.path.join(family_dir, cluster_dir, file))
                        if len(cluster2gbk[cluster_dir]) == 3:  # Break out of the loop if 3 items are stored
                            break
                if cluster2gbk.values():
                    dgk_dirs = list(cluster2gbk.values())[0]
                    cluster_num = cluster_dir.split("_")[1]
                    cluster_num = cluster_num.zfill(2)
                    family_num = '_'.join(subdir.split("_")[0:2])
                    output_pdf = f'{family_dir}/{family_num}_Cluster_{cluster_num}.pdf'

                    print(dgk_dirs)
                    print(output_pdf)
                    meta_info = f'{subdir} Cluster {cluster_dir.split("_")[1]}'
                    SVG(dgk_dirs,output_pdf,meta_info=meta_info)

                    shutil.copy(output_pdf, destination_dir)

    files = [f for f in os.listdir(destination_dir)]
    files = sorted(files, key=sort_key)
    file_dirs = [os.path.join(destination_dir, f) for f in files]

    # if os.path.isfile(merged_pdf):
    #     os.remove(merged_pdf)

    merger = PdfMerger()
    for pdf in file_dirs:
        merger.append(pdf)
    merger.write(merged_pdf)
    merger.close()