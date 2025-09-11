bind = "0.0.0.0:5000"    # 监听所有网络接口的5000端口
workers = 8               # 进程数（推荐值为CPU核数*2+1）
worker_class = "gevent"   # 使用异步处理[[3]][[7]][[11]]
timeout = 30              # 超时时间（秒）
accesslog = "-"           # 访问日志输出到控制台
errorlog = "-"            # 错误日志输出到控制台

