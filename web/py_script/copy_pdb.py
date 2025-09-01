import os
import shutil

# 定义源文件夹和目标文件夹路径
father_folder = "/media/Data/hebeibei/minced_output/array_cds_pairs_ncbi_10kb/BioPrinCRISPR_WebData/PDB_500-900"  # 源文件夹
target_folder = "/media/Data/hebeibei/minced_output/array_cds_pairs_ncbi_10kb/BioPrinCRISPR_WebData/PDB"         # 目标文件夹

# 确保目标文件夹存在（如果不存在则创建）
if not os.path.exists(target_folder):
    os.makedirs(target_folder)

# 获取源文件夹中的所有子文件夹
subfolders = os.listdir(father_folder)  # 

# 过滤出以"part"开头的子文件夹
part_subfolders = [sub for sub in subfolders if sub.startswith("part")]

# 遍历每个符合条件的子文件夹
for subfolder in part_subfolders:
    subfolder_path = os.path.join(father_folder, subfolder)
    
    # 递归遍历子文件夹中的所有文件（包括嵌套子文件夹）
    for root, dirs, files in os.walk(subfolder_path):  # 
        for file in files:
            source_path = os.path.join(root, file)       # 源文件完整路径
            destination_path = os.path.join(target_folder, file)  # 目标路径
            
            # 复制文件到目标文件夹（使用copy2保留元数据）
            if not os.path.exists(destination_path):
                shutil.copy2(source_path, destination_path)  
                print(f"已复制: {source_path} -> {destination_path}")

print("所有文件复制完成！")
