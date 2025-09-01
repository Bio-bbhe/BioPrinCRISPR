#连接 MongoDB
import json
from pymongo import MongoClient

client = MongoClient("localhost",27017)  # 默认地址
db = client["cytoscape_db"]  # 数据库名
collection_nodes = db["nodes"]  # 集合名

# 读取蛋白质ID映射文件
with open('./mongo/proteins_ids_map.txt', 'r') as f:
    proteins_map = json.load(f)

# 更新nodes集合
for pfam_id, protein_ids in proteins_map.items():
    result = collection_nodes.update_many(
        {"pfam_accession": pfam_id},
        {"$set": {"protein_ids": protein_ids}}
    )
    print(f'Updated {result.modified_count} documents for {pfam_id}')