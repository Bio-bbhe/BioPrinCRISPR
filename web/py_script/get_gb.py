# -*- coding: utf-8 -*-
import os
import shutil

# 1. 定义源目录（当前文件夹）和目标目录
source_dir = '/media/Data/hebeibei/minced_output/array_cds_pairs_ncbi_10kb/BioPrinCRISPR_WebData/06_data'  # 当前工作目录
target_dir = '/media/Data/hebeibei/minced_output/array_cds_pairs_ncbi_10kb/BioPrinCRISPR_WebData/gbk'  # 指定目标文件夹，可修改为实际路径

# 创建目标文件夹（如果不存在）
os.makedirs(target_dir, exist_ok=True)  # exist_ok=True避免重复创建错误 [[11]]

# 2. 遍历所有子文件夹及其文件
for root, dirs, files in os.walk(source_dir):  # os.walk递归遍历子目录 [[2,5,7,18]]
    for file in files:
        # 3. 仅处理.gb文件
        if file.endswith('.gb'):
            # 获取文件完整路径
            src_path = os.path.join(root, file)  # 构建源文件路径 [[4,11]]
            
            # 4. 提取新文件名（基于规则）
            try:
                # 移除扩展名（如.gb）
                base_name = os.path.splitext(file)[0]  # 分离文件名和扩展名 [[13]]
                
                # 用'__'分割文件名
                parts = base_name.split('__')  # 分割字符串以定位关键部分 [[8,13]]
                
                # 检查分割后部分数量（至少需3部分：前缀、目标段1、目标段2）
                if len(parts) >= 3:
                    # 取索引1和索引2（即两个'__'之间的部分），并用'__'连接
                    new_prefix = parts[1] + '__' + parts[2]  # 例如：'RFFG01000001.1__203677_203801_#12'
                    new_name = new_prefix + '.gb'  # 添加扩展名
                else:
                    # 如果不符合规则，保留原文件名（避免错误）
                    new_name = file
                    print(f"警告：文件 '{file}' 不符合命名规则，使用原文件名")
            except Exception as e:
                # 错误处理（如分割失败）
                new_name = file
                print(f"错误处理文件 '{file}': {e}, 使用原文件名")
            
            # 5. 构建目标路径并复制文件
            dest_path = os.path.join(target_dir, new_name)  # 目标文件完整路径 [[1,9,16]]
            shutil.copy2(src_path, dest_path)  # 复制文件（保留元数据） [[1,9,16]]

print("操作完成！所有.gb文件已复制并重命名到目标文件夹。")
