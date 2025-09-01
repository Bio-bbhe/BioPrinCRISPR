import os

# 定义输出文件名
output_file = 'output.txt'

# 使用with语句打开输出文件，确保文件正确关闭（参考Evidence 11和16关于文件写入的安全性）
with open(output_file, 'w') as f:
    # 遍历当前目录下的所有项（文件和文件夹），使用os.listdir('.')获取列表（参考Evidence 4和9）
    for item in os.listdir('/media/Data/hebeibei/minced_output/array_cds_pairs_ncbi_10kb/BioPrinCRISPR_WebData/PDB'):
        file_path = os.path.join('/media/Data/hebeibei/minced_output/array_cds_pairs_ncbi_10kb/BioPrinCRISPR_WebData/PDB', item)
        # 检查项是否为文件（而非文件夹），使用os.path.isfile()（参考Evidence 4和6）
        if os.path.isfile(file_path):
            # 检查文件是否为空（文件大小为0字节），使用os.path.getsize()
            if os.path.getsize(file_path) == 0:
                # 将文件名写入输出文件，每个文件名后添加换行符\n（参考Evidence 2、11和16）
                f.write(item + '\n')

# 打印完成消息
print(f"当前文件夹下的空文件名已输出到 {output_file}")

