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

des = pd.read_csv('/media/Data/liudong/pfam_acc2des.txt',sep='\t')
des_dict = des.set_index(des.columns[0]).to_dict()[des.columns[1]]


def arrow(X, Y, L, H, strand, h, l, color, GeneName):
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

    line = f'''<polygon id="{GeneName}" points="{A[0]},{A[1]} {B[0]},{B[1]} {C[0]},{C[1]} {D[0]},{D[1]} {E[0]},{E[1]} {F[0]},{F[1]} {G[0]},{G[1]}"
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


def draw_connector(x1, y1, x2, y2, font):
    '''
    Draw a connector line between two points with a right-angle bracket shape.
    '''
    mid_y = (y1 + y2) / 2 + 30
    y1 += 3
    connector = f'<path d="M {x1},{y1} L {x1},{mid_y} L {x2},{mid_y} L {x2},{y1}" style="stroke:black; stroke-width:0.5; fill:none;"/>'
    text_x = (x1 + x2) / 2 - 60
    text_y = mid_y + 15
    text = f'<text x="{text_x}" y="{text_y}" style="font-family: Arial; font-size:{font - 4};">Co-conserved array-protein</text>'
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


def SVG(GenBankFile,output_file, ArrowHeight=20, HeadEdge=8, HeadLength=10, marginX=10, marginY=30, scaling=30.0, font=16):
    '''
    Create the main SVG document with enhanced visualization:
        - read in GenBank document
        - find genes, start and stop positions, and strands
        - write the SVG files with different colors for each gene
    '''
    GenBankFile_name = GenBankFile.split('/')[-1].replace('tracrRNA_False__','').replace('tracrRNA_True__','')
    ALL_TEXT = ""
    dr_added = False
    gene_label = 'A'
    gene_pfam_mapping = {}
    gene_highlight = {}
    array_position = None
    gene_positions = {}
    MAX_Y = 0; # 动态调整svg图片的最大高度

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
                    strand = '-' if feature.location.strand == -1 else '+'
                    start = int(str(feature.location.start).replace('>', '').replace('<', '')) / scaling
                    stop = int(str(feature.location.end).replace('>', '').replace('<', '')) / scaling

                    # color = random_color()
                    # color = color_from_locus_tag(locus_tag)
                    color = color_from_pfam(pfam)
                    ALL_TEXT += arrow(start + marginX, marginY, stop - start, ArrowHeight, strand, HeadEdge, HeadLength, color, gene_label)

                    if previousLevel == 4: previousLevel = 0
                    text_y_offset = marginY + 50 + (font + 6) * previousLevel
                    text_x = marginX + start + (stop - start) / 2 - len(GeneName) * font / 4
                    text_id = "index1-" + gene_label
                    ALL_TEXT += f'<text id="{text_id}" x="{text_x}" y="{text_y_offset}" style="font-family: Arial; font-size:{font}; font-style: italic;">{GeneName}</text>'

                    # Add pfam number below the arrow
                    # if pfam:
                    gene_pfam_mapping[f'Gene {gene_label}'] = pfam if pfam else 'No Pfam domain'
                    gene_highlight[f'Gene {gene_label}'] = product
                    gene_positions[gene_label] = (marginX + start + (stop - start) / 2, marginY + ArrowHeight / 2)
                    gene_text_y_offset = marginY + ArrowHeight + 15
                    text_x = marginX + start + (stop - start) / 2 - len(gene_label) * font / 4
                    text_id = "index2-" + gene_label
                    ALL_TEXT += f'<text id="{text_id}" x="{text_x}" y="{gene_text_y_offset}" style="font-family: Arial; font-size:{font-4};font-weight: bold;">{gene_label}</text>'
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
                        text_x = text_x if text_x > 1 else 1
                        ALL_TEXT += f'<text x="{text_x}" y="{scale_y + 10}" style="font-family: Arial; font-size:10px; font-weight: bold;">Array</text>'
                        array_position = (marginX + start + (stop - start) / 2, marginY + ArrowHeight + 15)
    if array_position:
        for label, pfam in gene_pfam_mapping.items():
            if gene_highlight[label] == 'Putative Cas':
                gene_x, gene_y = gene_positions[label[-1]]
                array_x, array_y = array_position
                ALL_TEXT += draw_connector(array_x, array_y, gene_x, gene_y, font)

    mapping_x = marginX + 5
    mapping_y = marginY + 105
    MAX_Y = mapping_y
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
    max_chars = 40 # 设定最大显示字符数（不包括省略号）
    num_column = 3 #最大列数
    offset = 237;
        
    #每列的行数
    column_split = total_items // num_column + (1 if total_items % num_column != 0 else 0)  # Ensure enough rows for all items
    column_height = 0
    mapping_y_gene_label = mapping_y
    for idx, (k, v) in enumerate(gene_pfam_mapping.items()):
        last_k = k[-1] #截取字母
        condition = True if gene_highlight[k] == 'Putative Cas' else False
        v = ast.literal_eval(v) if v != 'No Pfam domain' else v
        annotations = [f'{pfam_id}: {des_dict[pfam_id]}' for pfam_id in v] if v != 'No Pfam domain' else ['No Pfam domain']

        # Calculate the offset for the columns
        current_col = idx // column_split  #列号
        column_offset = current_col * offset  # Offset to switch columns
        current_mapping_x = mapping_x + column_offset
        current_mapping_y = mapping_y_gene_label + (idx % column_split) * line_height

        show_k = k + " (Protein co-conserved with array)" if condition else k
        text_id = "info-" + last_k
        if (len(show_k) > (max_chars)):
            display_text = show_k[:max_chars - 3] + '...'
            ALL_TEXT += f'<text id="{text_id}" x="{current_mapping_x}" y="{current_mapping_y}" style="font-family: Arial; font-size:{font - 4}; font-weight: bold;"><title>{show_k}</title>{display_text}</text>'
        else:
            ALL_TEXT += f'<text id="{text_id}" x="{current_mapping_x}" y="{current_mapping_y}" style="font-family: Arial; font-size:{font - 4}; font-weight: bold;">{show_k}</text>'

        # Add annotations below the gene label
        for annotation in annotations:
            current_mapping_y += (font - 4)
            is_normal_weight = True if gene_highlight[k] == 'Putative Cas' else False

            # 如果文本超过最大长度，则截断并添加省略号，否则使用原文
            if (len(annotation) > max_chars):
                display_text = annotation[:max_chars - 3] + '...'  # 保留前17个字符，加上三个点
                ALL_TEXT += f'<text id="{text_id}" x="{current_mapping_x}" y="{current_mapping_y}" lengthAdjust="spacing" style="font-family: Arial; font-size:{font - 4}; {"font-weight: bold;" if is_normal_weight else ""}"><title>{annotation}</title>{display_text}</text>'
            else:
                ALL_TEXT += f'<text id="{text_id}" x="{current_mapping_x}" y="{current_mapping_y}" lengthAdjust="spacing" style="font-family: Arial; font-size:{font - 4}; {"font-weight: bold;" if is_normal_weight else ""}">{annotation}</text>'

        mapping_y_gene_label += max([(font - 4) * len(annotations),font])
        MAX_Y = max(current_mapping_y, MAX_Y)

        # Update column height and position if needed
        if idx % column_split == 0:
            column_height = max(column_height, mapping_y)
        if (idx + 1) % column_split == 0:
            mapping_y_gene_label = column_height  # Reset y-coordinate for the next column

    ALL_TEXT += '</svg>'
    HEAD = f"""<?xml version="1.0" standalone="no"?>
    <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" 
    "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">

    <svg width="195mm" height="{MAX_Y + 10}px" xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'>
    """
    ALL_TEXT = HEAD + ALL_TEXT

    with open(output_file, 'w') as f:
        f.write(ALL_TEXT)
    print(f"SVG saved to {output_file}")

    # convert_svg_to_pdf(output_file, output_file.replace('.svg', '.pdf'))
    # print(ALL_TEXT)

def extract_protein_id(filename):
    """
    从.gb文件名中提取proteins_id。
    格式示例：tracrRNA_False__DTQM01000097.1__10610_10821_#5__4_repeats_3478__gap.gb
    提取规则：取第一个"__"和第三个"__"之间的部分（如DTQM01000097.1__10610_10821_#5）。
    """
    # 分割文件名：使用"__"作为分隔符，得到各部分列表
    parts = filename.split('__')
    # 检查分割后是否有足够的部分（至少3个部分）
    if len(parts) < 3:
        raise ValueError(f"文件名 '{filename}' 格式无效，无法提取proteins_id")
    # 提取索引1和索引2的部分（第一个__后和第三个__前），并用"__"连接
    protein_id = f"{parts[1]}__{parts[2]}"
    return protein_id
    
def get_gb_files_absolute_paths(root_folder):
    """
    获取根目录下所有以"result_"开头的子文件夹中.gb文件的绝对路径
    :param root_folder: 根目录路径 (例如: "/test")
    :return: .gb文件的绝对路径列表
    """
    gb_files = []
    
    # 遍历根目录下的所有项目
    for entry in os.listdir(root_folder):
        entry_path = os.path.join(root_folder, entry)
        
        # 检查是否为"result_"前缀的文件夹
        if os.path.isdir(entry_path) and entry.startswith("result_"):
            # 递归遍历子文件夹中的文件 
            for foldername, subfolders, filenames in os.walk(entry_path):
                for filename in filenames:
                    # 检查.gb扩展名
                    if filename.endswith('.gb'):
                        # 获取绝对路径
                        abs_path = os.path.join(foldername, filename)
                        gb_files.append(os.path.abspath(abs_path))
    
    return gb_files


if __name__ == "__main__":
    gbk_dir = '/media/Data/hebeibei/minced_output/array_cds_pairs_ncbi_10kb/BioPrinCRISPR_WebData/06_data'
    svg_dir = '/media/Data/hebeibei/minced_output/array_cds_pairs_ncbi_10kb/BioPrinCRISPR_WebData/svg'
    gb_file_paths = get_gb_files_absolute_paths(gbk_dir)
    for gbk_file in gb_file_paths:
        protein_id = extract_protein_id(gbk_file.split("/")[-1])
        gbk2svg = f'{svg_dir}/{protein_id}.svg'
        print(gbk_file)
        print(gbk2svg)
        SVG(gbk_file, gbk2svg)

