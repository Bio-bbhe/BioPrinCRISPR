// src/utils/request.js
import axios from 'axios'

// 统一设置基础 URL
const service = axios.create({
  baseURL: process.env.VUE_APP_BASE_URL,
  timeout: 5000
})

export default service
