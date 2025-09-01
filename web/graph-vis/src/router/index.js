import Vue from 'vue'
import VueRouter from 'vue-router'  // Vue 2 必须显式调用 install 
import GlobalNetworkView from '../views/GlobalNetwork.vue'
import SubNetworkView from '../views/SubNetwork.vue'
import HomePageView from '../views/HomePage.vue'

Vue.use(VueRouter)  // 关键注册步骤！

const routes = [
  {
    path: '/',              // URL路径
    name: 'HomePage',           // 路由名称（可选）
    component: HomePageView     // 绑定的Vue组件
  },
  {
    path: '/network',              // URL路径
    name: 'GlobalNetwork',           // 路由名称（可选）
    component: GlobalNetworkView     // 绑定的Vue组件
  },
  {
    path: '/network/info',
    name: 'SubNetwork',
    component: SubNetworkView
  }
];
export default new VueRouter({ routes })
