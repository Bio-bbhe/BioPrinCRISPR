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
    pfam_dict = {}  # 主字典存储PFAM-ID映射
    id_tracker = {}  # 辅助字典用于高效去重（引用[[8]]）
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        # 跳过标题行
        for line in lines[1:]:
            columns = line.strip().split('\t')
            if len(columns) < 8:
                continue
                
            # 解析PFAM列表
            try:
                pfam_list = ast.literal_eval(columns[1])
            except (SyntaxError, ValueError):
                continue
                
            # 分割ID列表
            id_list = columns[7].split(';')
            
            # 更新字典（含去重）
            if len(pfam_list) > 1:
                for index, pfam in enumerate(pfam_list):
                    s_index = index - 1;
                    for i in range(2):
                        s_index += i;
                        if(s_index >= 0 and s_index + 1 < len(pfam_list)):
                            key = pfam_list[s_index] + ',' + pfam_list[s_index + 1]
                            # 初始化该PFAM的ID集合
                            if key not in pfam_dict:
                                pfam_dict[key] = []
                                id_tracker[key] = set()  # 辅助集合用于O(1)查重
                            
                            # 添加新ID（自动去重）
                            for id_val in id_list:
                                file_name = f"{id_val}.svg"
                                full_path = os.path.join(svg_path, file_name)
                                if os.path.exists(full_path) and (id_val not in id_tracker[key]):
                                    pfam_dict[key].append(id_val)
                                    id_tracker[key].add(id_val)
    
    return pfam_dict

if __name__ == "__main__":
    svg_path = '/media/Data/hebeibei/minced_output/array_cds_pairs_ncbi_10kb/BioPrinCRISPR_WebData/svg'
    result = parse_text_file("v1_rep_info_add_mapping_clean.tab", svg_path)
    
    # 写入JSON文件
    with open("proteins_ids_map.txt", 'w', encoding='utf-8') as out_f:
        json.dump(result, out_f, indent=2)
