from flask import Flask, jsonify, request
from pymongo import MongoClient
from flask_cors import CORS
import os
from config import get_config

app = Flask(__name__)
config = get_config()  # 获取当前环境配置
CORS(app, resources={r"/*": {"origins": config.ALLOWED_ORIGINS}})

# 连接MongoDB
client = MongoClient(config.MONGO_HOST, config.MONGO_PORT)  # 根据实际情况修改主机和端口
db_network = client[config.DB_NETWORK]  # 数据库名
db_info = client[config.DB_INFO]  # 数据库名
SVG_PATH = config.SVG_PATH
PDB_PATH = config.PDB_PATH
GBK_PATH = config.GBK_PATH

graph_data = None

@app.route('/load_data', methods=['GET'])
def load_data():
    global graph_data
    if(graph_data is not None):
        return jsonify({
            'status': 'success',
            'data': graph_data
        })
    # 查询节点数据
    nodes = list(db_network.nodes.find({}, {
        '_id': 0,  # 排除MongoDB的_id字段
        'id' : 1,
        'name': 1,
        'number_of_proteins': 1,
        'pfam_accession': 1,
        'x': 1,
        'y': 1,
    }))

    # 查询边数据
    edges = list(db_network.edges.find({}, {
        '_id': 0,
        'id': 1,
        'name': 1,
        'source': 1,
        'target': 1,
        'number_of_proteins': 1,
        'presence_status': 1
    }))

    # 转换数据格式
    graph_data = {
        'nodes': nodes,
        'edges': edges
    }

    return jsonify({
        'status': 'success',
        'data': graph_data
    })

@app.route('/sequence', methods=['GET'])
def get_seq():
    id = request.args.get('id')
    seq = db_info.proteins.find_one({'id': id},{'_id': 0, 'sequence': 1})
    if not seq:
        return jsonify({'status': 'error', 'message': '未找到匹配记录'}), 404
    return jsonify({
        'status': 'success',
        'data': seq['sequence']
    })

@app.route('/svg', methods=['GET'])
def get_svg():
    #获取并验证ID参数
    id = request.args.get('id')
    if not id:
        return jsonify({'status': 'error', 'message': 'ID参数缺失'}), 400
    
    #构建文件路径
    file_name = f"{id}.svg"
    full_path = os.path.join(SVG_PATH, file_name)
    
    #验证文件存在性
    if not os.path.exists(full_path):
        return jsonify({'status': 'error', 'message': f'未找到{file_name}文件'}), 404
    
    #读取文件内容
    with open(full_path, 'r', encoding='utf-8') as f:
        svg_content = f.read()
    
    return jsonify({
        'status': 'success',
        'data': svg_content
    })
    
@app.route('/svg/page', methods=['GET'])
def get_svg_page():
    # 获取并验证参数
    nodeId = request.args.get('id')
    pageNum = request.args.get('pageNum', type=int, default=1)
    pageSize = request.args.get('pageSize', type=int, default=10)
    
    if not nodeId:
        return jsonify({'status': 'error', 'message': 'ID参数缺失'}), 400

    # 查询节点数据
    proteinIds = list(db_network.edges.find(
        {
            '$or': [
                {'source': nodeId},
                {'target': nodeId}
            ]
        },
        {'_id': 0, 'protein_ids': 1}
    ))

    # 打平并去重（保持顺序）
    seen = set()
    merged_protein_ids = []
    for doc in proteinIds:
        for pid in doc.get('protein_ids', []):
            if pid not in seen:
                seen.add(pid)
                merged_protein_ids.append(pid)

    # 分页处理
    start = (pageNum - 1) * pageSize
    end = start + pageSize
    page_data = merged_protein_ids[start:end]
    
    svgs = []
    proteinIds = []
    total = len(merged_protein_ids)
    for proteinId in page_data:
        #构建文件路径
        file_name = f"{proteinId}.svg"
        full_path = os.path.join(SVG_PATH, file_name)
        
        #验证文件存在性
        if os.path.exists(full_path):
            #读取文件内容
            with open(full_path, 'r', encoding='utf-8') as f:
                svg_content = f.read()
                proteinIds.append(proteinId)
                svgs.append(svg_content)

    return jsonify({
        'status': 'success',
        'total': total,
        'pageNum': pageNum,
        'pageSize': pageSize,
        'proteinIds': proteinIds,
        'svgs': svgs
    })

    
    

@app.route('/pdb', methods=['GET'])
def get_pdb():
    #获取并验证ID参数
    id = request.args.get('id')
    if not id:
        return jsonify({'status': 'error', 'message': 'ID参数缺失'}), 400
    
    #构建文件路径
    file_name = f"{id}.pdb"
    full_path = os.path.join(PDB_PATH, file_name)
    
    #验证文件存在性
    if not os.path.exists(full_path):
        return jsonify({'status': 'error', 'message': f'未找到{file_name}文件'})
    
    #读取文件内容
    with open(full_path, 'r', encoding='utf-8') as f:
        pdb_content = f.read()
    
    return jsonify({
        'status': 'success',
        'data': pdb_content
    })
    
@app.route('/gbk', methods=['GET'])
def get_gbk():
    #获取并验证ID参数
    id = request.args.get('id')
    if not id:
        return jsonify({'status': 'error', 'message': 'ID参数缺失'}), 400
    
    #构建文件路径
    file_name = f"{id}.gb"
    full_path = os.path.join(GBK_PATH, file_name)
    
    #验证文件存在性
    if not os.path.exists(full_path):
        return jsonify({'status': 'error', 'message': f'未找到{file_name}文件'})
    
    #读取文件内容
    with open(full_path, 'r', encoding='utf-8') as f:
        gbk_content = f.read()
    
    return jsonify({
        'status': 'success',
        'data': gbk_content
    })

@app.route('/repeat', methods=['GET'])
def get_repeats():
    #获取并验证ID参数
    id = request.args.get('id')
    if not id:
        return jsonify({'status': 'error', 'message': 'ID参数缺失'}), 400

    repeats = db_info.repeats.find_one({'name': id},{'_id': 0})
    return jsonify({
        'status': 'success',
        'data': repeats
    })

    

@app.route('/load_data/<node_id>', methods=['GET'])
def load_data_by_node_id(node_id):
    # 查询节点
    node = db_network.nodes.find_one(
        {'id': node_id},
        {'_id': 0, 'id': 1, 'name': 1, 'number_of_proteins': 1, 'pfam_accession': 1}
    )
    
    if not node:
        return jsonify({'status': 'error', 'message': '未找到匹配记录'}), 404
        
    # 查询关联边
    edges = list(db_network.edges.find({
        '$or': [
            {'source': node_id},
            {'target': node_id}
        ]
    }, {
        '_id': 0,
        'id': 1,
        'name': 1,
        'source': 1,
        'target': 1,
        'number_of_proteins': 1,
        'presence_status': 1
    }))

    # 获取所有关联节点ID
    neighbor_ids = {*({edge['source'] for edge in edges} | {edge['target'] for edge in edges})} - {node['id']}
    
    # 查询邻居节点
    neighbors = []
    if neighbor_ids:
        neighbors = list(db_network.nodes.find(
            {'id': {'$in': list(neighbor_ids)}},
            {'_id': 0, 'id': 1, 'name': 1, 'number_of_proteins': 1, 'pfam_accession': 1}
        ))

    return jsonify({
        'status': 'success',
        'data': {
            'node': node,
            'edges': edges,
            'neighbors': neighbors
        }
    })

@app.route('/load_data/batch', methods=['POST'])
def load_data_batch_ids():
    data = request.get_json()
    if not data or 'ids' not in data:
        return jsonify({'status': 'error', 'message': 'Missing ID array'}), 400

    # 查询节点数据
    edges = list(db_network.edges.find(
        {
            '$or': [  # 逻辑 OR 操作符，满足任一条件即可
                {'source': {'$in': data['ids']}},  # source 字段值在 data['ids'] 中
                {'target': {'$in': data['ids']}}   # target 字段值在 data['ids'] 中
            ]
        },
        {'_id': 0, 'id': 1, 'source': 1, 'target': 1, 'source_pfam_accession': 1, 'target_pfam_accession': 1, 'protein_ids': 1}
    ))

    # 处理蛋白质域数据
    id_domains = {}
    for edge in edges:
        if edge['source'] not in id_domains:
                id_domains[edge['source']] = []
        if edge['target'] not in id_domains:
                id_domains[edge['target']] = []
        proteinIds = edge['protein_ids']
        if len(proteinIds) > 0:
            domains = {
                'protein_id': proteinIds[0],
                'source': edge['source_pfam_accession'],
                'target': edge['target_pfam_accession'],
            }
            id_domains[edge['source']].append(domains)
            if edge['target'] != edge['source']:
                id_domains[edge['target']].append(domains)

    return jsonify({
        'status': 'success',
        'data': id_domains,
    })

@app.errorhandler(404)  # 统一处理未找到资源
def handle_404(e):
    return jsonify({'status': 'error', 'message': f'资源不存在: {str(e)}'}), 404

@app.errorhandler(Exception)  # 捕获所有未处理异常
def handle_generic_error(e):
    return jsonify({'status': 'error', 'message': f'服务器内部错误: {str(e)}'}), 500

if __name__ == '__main__':
    app.run(debug=config.DEBUG, host=config.HOST, port=config.PORT)  # 启动Flask服务
