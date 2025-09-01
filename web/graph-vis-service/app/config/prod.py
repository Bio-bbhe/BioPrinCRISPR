class Config:
    DEBUG = False
    # MongoDB 配置
    MONGO_HOST = "mongo" #ip:mongo容器名
    MONGO_PORT = 27017
    DB_NETWORK = "network_big"
    DB_INFO = "protein_info"
    
    # 文件路径
    SVG_PATH = '/svg'
    PDB_PATH = '/PDB'
    GBK_PATH = '/gbk'
    
    # CORS 配置
    ALLOWED_ORIGINS = "http://vue-app:8080"  # 生产前端地址，ip:vue-app容器名
    
    # Flask 服务配置
    HOST = '127.0.0.1'
    PORT = 5000
