module.exports = {
  configureWebpack: {
    devtool: 'source-map' // 关键配置，否则无法映射原始代码 
  },
  transpileDependencies: true,
}
