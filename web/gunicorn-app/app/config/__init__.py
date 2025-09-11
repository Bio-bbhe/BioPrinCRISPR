import os

def get_config():
    env = os.environ.get('FLASK_ENV', 'dev').lower() #默认dev环境
    
    if env == 'prod':
        from .prod import Config
    else:  # 默认使用开发配置
        from .dev import Config
    
    return Config
