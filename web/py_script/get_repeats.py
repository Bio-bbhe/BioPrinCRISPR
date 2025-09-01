import os
import json

# 设置目标文件夹路径（请根据实际情况修改）
folder_path = "/media/Data/hebeibei/minced_output/array_cds_pairs_ncbi_10kb/BioPrinCRISPR_WebData/forSpacers/01.mincedout"  # 替换为您的文件夹路径
#folder_path = "./"
# 初始化结果列表
result = []

# 遍历文件夹中的所有文件 [[6, 7, 12]]
for filename in os.listdir(folder_path):
    file_path = os.path.join(folder_path, filename)
    
    # 检查是否为TXT文件
    if os.path.isfile(file_path) and filename.endswith(".txt"):
        # 读取文件内容并分割成行 
        with open(file_path, "r") as file:
            lines = file.read().splitlines()
        
        # 检查行数：如果小于6，跳过此文件
        if len(lines) < 6:
            continue
        
        # 提取第6行到倒数第5行（索引从0开始：第6行索引5，倒数第5行索引len(lines)-5）
        # 使用切片 lines[5:-4] 包括索引5到索引len(lines)-5（不包括len(lines)-4）
        data_lines = lines[5:-6]
        
        # 初始化items列表存储当前文件的数据
        items = []
        
        # 处理每一行数据
        for line in data_lines:
            # 跳过空行
            if not line.strip():
                continue
            
            # 分割行：最多分割成四部分（position, repeat, spacer_rest）
            parts = line.split(None, 3)  # split(None, 2) 使用空白分割，最多四部分
            
            # 解析列：position转换为整数，repeat和spacer作为字符串
            position = None
            repeat = ""
            spacer = ""
            
            if len(parts) >= 1:
                try:
                    position = int(parts[0])  # position应为整数
                except ValueError:
                    position = None  # 如果转换失败，跳过此行
            
            if position is not None:
                if len(parts) >= 2:
                    repeat = parts[1] + "(" + str(len(parts[1]))+ "bp)"
                if len(parts) >= 3:
                    spacer = parts[2].strip() + "(" + str(len(parts[2]))+ "bp)" # 移除spacer首尾空白，但保留内部格式
                else:
                    spacer = ""  # 如果spacer列缺失，设为空字符串
                
                # 添加字典到items
                items.append({
                    "position": position,
                    "repeat": repeat,
                    "spacer": spacer
                })
                
        #解析最后一行（总结行）
        last_line = lines[-5].strip()  # 获取最后一行并去除首尾空白
        repeat_count = 0
        avg_repeat_len = 0
        avg_spacer_len = 0
        
        # 使用split()分割单词并提取数值
        words = last_line.split()
        if len(words) >= 8:  # 确保有足够的分词
            try:
                # 提取三个数值（位置固定）
                repeat_count = int(words[1])          # "Repeats:"后的第一个数字
                avg_repeat_len = int(words[4])         # 第一个"Average Length:"后的数字
                avg_spacer_len = int(words[7])         # 第二个"Average Length:"后的数字
            except (ValueError, IndexError):
                # 如果解析失败，保持默认值0
                pass
        
        # 处理文件名：移除 .minced.txt 后缀以获取name
        base_name = os.path.basename(filename)  # 获取基本文件名（不带路径）
        name_without_txt = os.path.splitext(base_name)[0]  # 移除 .txt 后缀
        name = os.path.splitext(name_without_txt)[0]  # 移除 .minced 后缀
        
        # 添加当前文件的数据到结果列表
        result.append({
            "name": name,
            "items": items,
            "repeat_count": repeat_count,
            "avg_repeat_len": avg_repeat_len,
            "avg_spacer_len": avg_spacer_len
        })

# 写入结果文件
output_file="repeats.json"
with open(output_file, 'w') as out_f:
    json.dump(result, out_f, indent=4)
