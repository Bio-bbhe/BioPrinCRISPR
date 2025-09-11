# 使用官方Python轻量镜像
FROM python:3.8-slim

# 设置工作目录
WORKDIR /app

# 复制依赖文件并安装
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 复制应用代码
COPY . .

# 设置环境变量
ENV FLASK_ENV=prod
ENV PYTHONPATH=/app/app

# 暴露端口（需与Gunicorn配置一致）
EXPOSE 5000

# 启动命令（核心！）
CMD ["gunicorn", "-c", "gunicorn.conf.py", "app:app"] 

