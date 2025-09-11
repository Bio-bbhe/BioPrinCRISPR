# 基于官方Nginx镜像（可指定版本，如 nginx:1.25）
FROM nginx:latest

# 删除Nginx默认配置（避免冲突）
RUN rm /etc/nginx/conf.d/default.conf 

# 复制自定义Nginx配置到容器
COPY nginx.conf /etc/nginx/conf.d/

# 复制Vue静态资源到容器
COPY dist/ /usr/share/nginx/html/ 

# 暴露端口（Nginx默认监听80）
EXPOSE 8080

