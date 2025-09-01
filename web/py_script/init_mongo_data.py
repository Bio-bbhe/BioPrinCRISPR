from pymongo import MongoClient
import json

#连接 MongoDB
client = MongoClient("172.19.196.58",27017)  # 默认地址
db_info = client["protein_info"]  # 数据库名
collection_proteins = db_info["proteins"]  # 蛋白质序列集合
collection_repeats = db_info["repeats"]  # 重复序列集合

def insert_network(db_name, cyjs_path, proteins_map):
     # 解析 Cytoscape 数据文件
    with open(cyjs_path, "r") as f:
        cytoscape_data = json.load(f)

    db = client[db_name]  # 数据库名
    collection_nodes = db["nodes"]  # 节点集合
    collection_edges = db["edges"]  # 边集合

    #清空数据库
    client.drop_database(db.name)
    print(f"数据库 {db.name} 已清空")

    #提取并重组节点数据
    nodes_to_insert = []
    for node in cytoscape_data["elements"]["nodes"]:
        node_data = node["data"]
        node_position = node["position"]

        # 合并 data 和 position 字段
        combined_data = {
            **node_data,  # 展开 data 所有字段
            "x": node_position["x"],  # 坐标单独存储
            "y": node_position["y"],
        }
        nodes_to_insert.append(combined_data)

    #批量插入 MongoDB
    if nodes_to_insert:
        result = collection_nodes.insert_many(nodes_to_insert)
        print(f"成功插入 {len(result.inserted_ids)} 个节点")
    else:
        print("未找到有效节点数据")

    #提取并重组边数据
    edges_to_insert = []
    for edge in cytoscape_data["elements"]["edges"]:
        edge_data = edge["data"]
        edge_data["source_pfam_accession"] = edge_data["pfam_accession"].split(",")[0]
        edge_data["target_pfam_accession"] = edge_data["pfam_accession"].split(",")[1]
        edge_key = f"{edge_data['source_pfam_accession']},{edge_data['target_pfam_accession']}"

        # 合并 data 和 position 字段
        combined_data = {
            **edge_data,  # 展开 data 所有字段
            "protein_ids": proteins_map.get(edge_key, [])
        }
        edges_to_insert.append(combined_data)

    #批量插入 MongoDB
    if edges_to_insert:
        result = collection_edges.insert_many(edges_to_insert)
        print(f"成功插入 {len(result.inserted_ids)} 个边")
    else:
        print("未找到有效边数据")

try:
    # 读取蛋白质ID映射文件
    with open("./proteins_ids_map.txt", 'r') as f:
        proteins_map = json.load(f)

    #加载网络数据到mongo
    insert_network("network_big", "./domain_p-0.0_m-5.graphml.cyjs", proteins_map)
    insert_network("network_small", "./domain_p-0.0_m-5_1e-10.graphml.cyjs", proteins_map)

    # 解析 sequence 数据文件
    with open("./proteins_sequence_map.txt", "r") as f:
        seq_map = json.load(f)

    # 解析重复序列文件
    with open("./repeats.json", "r") as f:
        repeat_map = json.load(f)

    #清空数据库
    client.drop_database(db_info.name)
    print(f"数据库 {db_info.name} 已清空")

    result = collection_proteins.insert_many(seq_map)
    print(f"成功插入 {len(result.inserted_ids)} 个蛋白质序列")

    result = collection_repeats.insert_many(repeat_map)
    print(f"成功插入 {len(result.inserted_ids)} 个重复序列")


finally:
    # 关闭连接
    client.close()
