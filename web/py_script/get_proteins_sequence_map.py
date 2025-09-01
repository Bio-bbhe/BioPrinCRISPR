import ast
import json
import os

def parse_text_file(file_path, svg_path):
    """
    读取文本文件，解析PfamAnnotation和IDs列，返回去重后的PFAM代码与IDs字典
    
    参数:
        file_path (str): 文件路径
    返回:
        dict: 键为PFAM代码，值为去重后的IDs列表
    """
    id_dict = {}  # 主字典存储PFAM-ID映射
    json_arr = []
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        # 跳过标题行
        for line in lines[1:]:
            columns = line.strip().split('\t')
            if len(columns) < 8:
                continue
                
            # 分割ID列表
            id_list = columns[7].split(';')
            seq_list = columns[8].split(';')
            
            # 更新字典（含去重）
            for index, id in enumerate(id_list):
                # 初始化该PFAM的ID集合
                if id not in id_dict:
                    id_dict[id] = seq_list[index]
                    file_name = f"{id}.svg"
                    full_path = os.path.join(svg_path, file_name)
                    svg_status = os.path.exists(full_path)
                    entry = {
                        'id': id,
                        'sequence': seq_list[index],
                        'svg_status': svg_status
                    }
                    json_arr.append(entry)
    
    return json_arr

if __name__ == "__main__":
    svg_path = '/media/Data/hebeibei/minced_output/array_cds_pairs_ncbi_10kb/BioPrinCRISPR_WebData/svg'
    result = parse_text_file("v1_rep_info_add_mapping_clean.tab", svg_path)
    
    # 写入JSON文件
    with open("proteins_sequence_map.txt", 'w', encoding='utf-8') as out_f:
        json.dump(result, out_f, indent=2)
