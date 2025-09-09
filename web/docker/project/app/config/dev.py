class Config:
    DEBUG = True
    # MongoDB 配置
    MONGO_HOST = "172.19.196.58"
    MONGO_PORT = 27017
    DB_NETWORK = "network_big"  # 测试数据库
    DB_INFO = "protein_info"
    
    # 文件路径
    SVG_PATH = '/media/Data/hebeibei/minced_output/array_cds_pairs_ncbi_10kb/BioPrinCRISPR_WebData/svg'
    PDB_PATH = '/media/Data/hebeibei/minced_output/array_cds_pairs_ncbi_10kb/BioPrinCRISPR_WebData/PDB'
    GBK_PATH = '/media/Data/hebeibei/minced_output/array_cds_pairs_ncbi_10kb/BioPrinCRISPR_WebData/gbk'
    
    # CORS 配置
    ALLOWED_ORIGINS = "http://172.19.97.72:8080"  # 本地测试
    
    # Flask 服务配置
    HOST = '0.0.0.0'
    PORT = 5001
