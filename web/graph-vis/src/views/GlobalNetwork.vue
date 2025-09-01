<template>
  <div class="sigma-container">
    
    <!-- 网络图容器 -->
    <div ref="container" class="graph-container" v-loading="networkLoading">
      <!-- 新增交互组件 -->
      <div class="interact-box">
        <div style="margin-bottom: 2px">
          number of proteins
        </div>
        <div style="display: flex; flex-direction: row;">
          <el-image
            :src="'/gep.svg'"
            style="margin-right: 8px;flex: 1;"
            fit="contain"
          />
          <el-input
            v-model.number="inputNumber"
            placeholder="0"
            size="small"
            style="margin-right: 10px;flex: 5;"
            type="number"
          />
          <el-button size="small" icon="el-icon-search" circle @click="handleNopInput" style="flex: 0.5;"></el-button>
        </div>
      </div>
      <div class="select-box">
        <el-select
          v-model="selectedPfam"
          filterable
          placeholder="请选择pfam_accession"
          @change="handlePfamSelect"
          clearable
          style="flex: 8;margin-right: 10px;"
        >
          <el-option
            v-for="item in pfamOptions"
            :key="item.value"
            :label="item.label"
            :value="item.value"
          />
        </el-select>
        <el-button size="small" plain @click="handleNodeDetail" style="flex: 1;">查看详情</el-button>
      </div>
      <el-button
        class="reset-btn"
        size="small"
        circle
        icon="el-icon-refresh-right"
        @click="handleReset"
      />
    </div>
  </div>
</template>

<script setup>
import { ref, onMounted, onBeforeUnmount } from 'vue';
import Graph from 'graphology';
import Sigma from '@/libs/sigma/sigma' //@ 是 Vue CLI / Vite 默认配置的别名，指向 src/
import { circular } from 'graphology-layout';
import { useRouter } from 'vue-router/composables'
import request from '@/utils/request';
import { Message } from 'element-ui';

const router = useRouter()
// 使用Composition API声明响应式变量
const container = ref(null);
const sigmaInstance = ref(null);
const graph = ref(new Graph());
let highlightNodes = [];
let selfEdges = {}; //记录自连接的节点id（id : number of proteins）
let edgesMap = {}; //记录所有起点边的终点数组
let doubleEdges = {}; //记录双向边（其中一条边）
let nodeNumEdgeMap = {}; //记录每个节点关联的边数（节点是边的起点或终点的次数，自连接边记1次）
const selectedPfam = ref('');
const pfamOptions = ref([]);
const inputNumber = ref(0);
let commonEdgeList = [] //普通边（不包含自连接边和双向边）
let selectedNodeId = null;
let curRatio = 1;
const networkLoading = ref(true);
const baseNodeSize = 3; // 初始尺寸
const baseEdgeSize = 1; // 初始尺寸
const baseLabelSize = 8;

// 初始化图表结构
const initGraph = () => {
  // 清空现有图表
  graph.value.clear();
  fetchGraphData();
  
  // 应用圆形布局
  circular.assign(graph.value);
}

// 从后端API获取数据
const fetchGraphData = async () => {
  networkLoading.value = true;
  try {
    const response = await request.get('/load_data');
    const { nodes, edges } = response.data.data;
    
    // 批量添加节点
    nodes.forEach(node => {
      const lines = node.name.split('\n');
      const color = getNodeColor(node.number_of_proteins);
      nodeNumEdgeMap[node.id] = 0;
      graph.value.addNode(node.id, {
        label: lines[0],
        x: node.x ,
        y: -node.y, //解决坐标系不匹配问题
        size: baseNodeSize,
        nop: lines[1],
        bpsl: lines[2],
        color: color,
        default_color: color,
        num_of_proteins: node.number_of_proteins,
        pfam_accession: node.pfam_accession
      });
      pfamOptions.value.push({
        value: node.id,
        label: node.pfam_accession
      });
    })

    for (let i = 0; i < edges.length; i++) {
      const edge = edges[i];
      nodeNumEdgeMap[edge.source] ++;
      if (edge.source === edge.target) {
        selfEdges[edge.source] = {
          nop: edge.number_of_proteins,
          color: getEdgeColor(edge.presence_status),
          hidden: false
        };
      } else {
        nodeNumEdgeMap[edge.target] ++;
        let targetArr = edgesMap[edge.source];
        if(!targetArr) {
          targetArr = [];
          edgesMap[edge.source] = targetArr
        }
        targetArr.push(edge.target);
        const sourceArr = edgesMap[edge.target];
        let flag = true;
        if(sourceArr){
          for (let i = 0; i < sourceArr.length; i++) {
            if (sourceArr[i] === edge.source) {
              let doubleEdgeArr = doubleEdges[edge.source];
              if(!doubleEdgeArr) {
                doubleEdgeArr = [];
                doubleEdges[edge.source] = doubleEdgeArr;
              }
              doubleEdgeArr.push({
                target: Number(edge.target),
                nop: edge.number_of_proteins,
                status: edge.presence_status,
                color: getEdgeColor(edge.presence_status),
                hidden: false
              });
              flag = false;
              break;
            }
          }
        }
        if(flag){
          graph.value.addEdge(edge.source, edge.target, {
            size: baseEdgeSize,
            label: edge.number_of_proteins,
            status: edge.presence_status,
            type: 'arrow',
            color: getEdgeColor(edge.presence_status),
          });
          commonEdgeList.push(edge);
        }
      }
    }
    // 刷新渲染器
    if (sigmaInstance.value) {
      sigmaInstance.value.refresh()
    }
    
  } catch (error) {
    console.error('数据加载失败:', error)
    Message.error({
      showClose: true,
      message: '数据加载失败',
    });
  }
  networkLoading.value = false;
}

// 初始化Sigma渲染器
const initSigma = () => {
  if (container.value) {
    // 销毁现有实例
    if (sigmaInstance.value) {
      sigmaInstance.value.kill()
    }
    
    // 创建新实例
    sigmaInstance.value = new Sigma(
      graph.value, 
      container.value,
      {
        renderer: {
          type: "canvas", // 使用Canvas渲染器
        },
        renderLabels : false,
        defaultNodeColor: '#FFEDA0',
        defaultEdgeColor: '#4563F8', //#444bf8
        labelFont : 'lighter',
        labelSize: 8,
        maxCameraRatio : 1,
        labelDensity: 10, //禁用标签密度过滤，强制显示所有节点的标签
        selfEdges: selfEdges,
        doubleEdges: doubleEdges,
        edgeLabelColor: {//边标签默认颜色为黑色
          color: "#000",
        },
        edgeLabelSize: 10,
      }
    );
    curRatio = sigmaInstance.value.getCamera().ratio;

    sigmaInstance.value.setSetting('getEdgeColor', getEdgeColor);

    console.log(sigmaInstance.value);

    //关键：空实现circle节点的悬停渲染器，避免节点内部的label被遮挡
    sigmaInstance.value.nodeHoverPrograms['circle'].render = () =>{};

    //关键：空实现悬停渲染器，在labelRenderer完成悬停高亮效果
    sigmaInstance.value.setSetting('hoverRenderer', () => {

    });

    // sigmaInstance.value.on('beforeRender', () => {
    // });

    //添加相机缩放监听
    sigmaInstance.value.getCamera().on('updated', () => {
      const ratio = sigmaInstance.value.getCamera().ratio;
      if(curRatio == ratio)  //相机缩放比例未发生变化，直接返回
        return;
      curRatio = ratio;
      // 动态调整节点尺寸（ratio越小尺寸越大）
      graph.value.updateEachNodeAttributes((node, attr) => {
        return {
          ...attr,
          size: baseNodeSize * (1 + (1 - ratio)), // 缩放因子公式
        };
      });

      // 设置可见性阈值（缩放比例小于0.3时显示边标签）
      sigmaInstance.value.setSetting('renderEdgeLabels', ratio < 0.09);
      sigmaInstance.value.setSetting('renderLabels', ratio < 0.09);
      sigmaInstance.value.setSetting('labelSize', baseLabelSize * (1 + (1 - ratio) ** 10));
      const edgeLabelSize = sigmaInstance.value.getSetting('edgeLabelSize')
      if(ratio < 0.06){
        if(edgeLabelSize != 15) sigmaInstance.value.setSetting('edgeLabelSize', 15);
      }else if(ratio < 0.15){
        if(edgeLabelSize != 12) sigmaInstance.value.setSetting('edgeLabelSize', 12);
      }else if(ratio < 0.3){
        if(edgeLabelSize != 10) sigmaInstance.value.setSetting('edgeLabelSize', 10);
      }
    });

    // 自定义node标签渲染器
    sigmaInstance.value.setSetting('labelRenderer', (context, data, settings) => {
      const size = settings.labelSize;
      const font = settings.labelFont;
      const weight = settings.labelWeight;
      context.textAlign = 'center';
      const needHover = data.hovered || data.selected;

      //高亮显示：增加白色边框和阴影
      if (needHover || data.highlighted) { 
        context.lineWidth = 5;
        const PADDING = 1.5;
        context.beginPath();
        context.strokeStyle = "#FFF";
        context.shadowOffsetX = 0;
        context.shadowOffsetY = 0;
        context.shadowBlur = 6; // 增大阴影扩散范围
        context.shadowColor = "rgba(0,0,0,0.80)"; // 调整阴影透明度
        context.arc(data.x, data.y, data.size + PADDING, 0, Math.PI * 2);
        context.closePath();
        context.stroke(); // 仅保留描边
      }

      //绘制节点实心圆，防止节点被自连接边和双向边遮挡
      context.beginPath();
      context.fillStyle = data.color;
      context.shadowColor = "rgba(0,0,0,0.0)"; // 调整阴影透明度
      context.arc(data.x, data.y, data.size, 0, Math.PI * 2);
      context.closePath();
      context.fill(); 
  
      // 显示主标签
      const verticalPos = size / 3;
      context.fillStyle = "#000"; //黑字
      context.font = `${weight} ${size}px ${font}`;
      let offset = needHover ? verticalPos * 2 : 0;
      const label = data.label;
      const length = label.length;
      if(length > data.size / 3.85){ //若标签过长，分两行显示
        let label1 = label.substring(0, length / 2);
        let label2 = label.substring(length / 2, length);
        context.fillText(label1, data.x, data.y - verticalPos - offset);
        context.fillText(label2, data.x, data.y + size - offset);
        offset = size - offset;
      }else{
        context.fillText(data.label, data.x, data.y + verticalPos - offset);
        offset = verticalPos - offset;
      }

      if (needHover) {    
        context.fillStyle = "#000"; //黑字       
        // 第二行信息（nop）
        context.font = `${weight} ${size * 0.8}px ${font}`;
        context.fillText(data.nop, data.x, data.y + offset + size * 1.2);
        
        // 第三行信息（bpsl）
        context.fillText(data.bpsl, data.x, data.y + offset + size * 2.2);
      } 
    });
    
    // 鼠标单击节点高亮
    sigmaInstance.value.on('clickNode', (event) => {
      const nodeId = event.node;
      console.log("clickNode", nodeId)
      if(selectedNodeId && nodeId !== selectedNodeId){
        clearHighlightNodes();
        graph.value.setNodeAttribute(selectedNodeId, 'selected', false);
      }
      selectedNodeId = nodeId;
      selectedPfam.value = selectedNodeId;
      graph.value.setNodeAttribute(nodeId, 'selected', true);
      highlightNode(nodeId)
    });

    // 鼠标双击进入节点详情页
    sigmaInstance.value.on('doubleClickNode', (event) => {
      event.preventSigmaDefault(); //阻止Sigma内置的双击缩放行为
      console.log('doubleClickNode', event.node);
      openNodeDetail(event.node);
    });

    sigmaInstance.value.on('enterNode', (event) => {
      //forceLabel（系统属性），关键：强制显示label，避免节点缩放太大，标签不可见，无法触发自定义标签渲染器
      graph.value.setNodeAttribute(event.node, 'forceLabel', true);
      //hovered（自定义属性），关键：鼠标悬停时显示多行label
      graph.value.setNodeAttribute(event.node, 'hovered', true);
    });

    sigmaInstance.value.on('leaveNode', (event) => {
      graph.value.setNodeAttribute(event.node, 'forceLabel', false);
      graph.value.setNodeAttribute(event.node, 'hovered', false);
    });

    // 点击空白画布，取消所有节点高亮
    sigmaInstance.value.on('clickStage', () => {
      clearHighlightNodes();
      selectedPfam.value = '';
      selectedNodeId = null;
    });
  }
}

const clearHighlightNodes = () => {
  highlightNodes.forEach(nodeId => {
    graph.value.setNodeAttribute(nodeId, "color", graph.value.getNodeAttribute(nodeId, "default_color"));
    graph.value.setNodeAttribute(nodeId, "highlighted", false);
    graph.value.setNodeAttribute(nodeId, "selected", false);
  });
  highlightNodes.length = 0;
}

const openNodeDetail = (nodeId) => {
  const label = graph.value.getNodeAttribute(nodeId, "label");
      const routeData = router.resolve({
        path: '/network/info',
        query: { 
          id: nodeId,
          label
        }
      })
      window.open(routeData.href, '_blank')
}

const handleNodeDetail = () => {
  if(selectedNodeId){
    openNodeDetail(selectedNodeId);
  }
}

// 高亮选中节点
const highlightNode = (nodeId) => {
  highlightNodes.push(nodeId);
  const neighbors = graph.value.neighbors(nodeId);
  // graph.value.setNodeAttribute(nodeId, "color", "#59f652");
  graph.value.setNodeAttribute(nodeId, "highlighted", true);
  neighbors.forEach(neighbor => {
      highlightNodes.push(neighbor);
      // graph.value.setNodeAttribute(neighbor, "color", "#59f652");
      graph.value.setNodeAttribute(neighbor, "highlighted", true);
    });
}

//根据number_of_proteins大小设置节点颜色
const getNodeColor = (number_of_proteins) => {
  if (number_of_proteins < 300) {
    return '#FFEA9C'; 
  } else if (number_of_proteins <= 500) {
    return '#FFDC87'; 
  }else if (number_of_proteins <= 1000) {
    return '#FFD178'; 
  }else if (number_of_proteins <= 2000) {
    return '#FEB754'; 
  }else if (number_of_proteins <= 5000) {
    return '#FEA447'; 
  }else if (number_of_proteins <= 10000) {
    return '#FD8F40'; 
  }else if (number_of_proteins <= 20000) {
    return '#FD6A33'; 
  } else {
    return '#F7482A'; 
  }
}

const getEdgeColor = (status) => {
  if (status == 'Found') {
    return '#F7482A'; 
  }else {
    return '#4A6FF8';
  }
}

// 处理选择事件
const handleNopInput = () => {
  const tempNodeNumEdgeMap = {...nodeNumEdgeMap};
  if (inputNumber.value >= 0) {
    commonEdgeList.forEach(edge => {
      const nop = edge.number_of_proteins;
      if(nop < inputNumber.value){
        graph.value.setEdgeAttribute(edge.source, edge.target, "hidden", true);
        tempNodeNumEdgeMap[edge.source]--;
        tempNodeNumEdgeMap[edge.target]--;
      }else{
        graph.value.setEdgeAttribute(edge.source, edge.target, "hidden", false);
      }
    });
    Object.entries(selfEdges).forEach(([source, selfEdge]) => {
      if(selfEdge.nop < inputNumber.value){
        selfEdge.hidden = true;
        tempNodeNumEdgeMap[source]--;
      }else{
        selfEdge.hidden = false;
      }
    })
    Object.entries(doubleEdges).forEach(([source, doubleEdgeArr]) => {
      doubleEdgeArr.forEach(doubleEdge => {
        if(doubleEdge.nop < inputNumber.value){
          doubleEdge.hidden = true;
          tempNodeNumEdgeMap[source]--;
          tempNodeNumEdgeMap[doubleEdge.target]--;
        }else{
          doubleEdge.hidden = false;
        }
      })
    })
    pfamOptions.value = [];
    delete tempNodeNumEdgeMap['undefined'];
    Object.entries(tempNodeNumEdgeMap).forEach(([key, value]) => {
      if(value <= 0){
        graph.value.setNodeAttribute(key, "hidden", true);
      }else{
        graph.value.setNodeAttribute(key, "hidden", false);
        pfamOptions.value.push({
          value: key,
          label: graph.value.getNodeAttribute(key, "pfam_accession")
        });
      }
    });
    if(selectedNodeId && graph.value.getNodeAttribute(selectedNodeId, "hidden")){
      selectedNodeId = null;
      selectedPfam.value = '';
      graph.value.setNodeAttribute(selectedNodeId, 'selected', false);
    }
  }
};

const handleReset = () => {
  highlightNodes = [];
  selfEdges = {};
  edgesMap = {};
  doubleEdges = {};
  nodeNumEdgeMap = {};
  pfamOptions.value = [];
  inputNumber.value = 0;
  commonEdgeList = [];
  selectedPfam.value = '';
  selectedNodeId = null;
  initGraph();
  sigmaInstance.value.getCamera().animate({
    x: 0.5,
    y: 0.5,
    ratio: 1
  }, {
    duration: 500
  });
  if (sigmaInstance.value) {
      sigmaInstance.value.refresh();
  }
};

const handlePfamSelect = (nodeId) => {
  if(selectedNodeId){
    graph.value.setNodeAttribute(selectedNodeId, 'selected', false);
  }
  if (!nodeId) {
    selectedNodeId = null;
    return;
  }
  selectedNodeId = nodeId;
  graph.value.setNodeAttribute(nodeId, 'selected', true);
  const node = sigmaInstance.value.nodeDataCache[nodeId]
  sigmaInstance.value.getCamera().animate({
    x: node.x,
    y: node.y,
    ratio: 0.08
  }, {
    duration: 500
  });
};

// 生命周期钩子
onMounted(() => {
  document.title = 'Network';
  initGraph();
  initSigma();
})

onBeforeUnmount(() => {
  if (sigmaInstance.value) {
    sigmaInstance.value.kill()
  }
})
</script>

<style scoped>
.sigma-container {
  position: relative;
  width: 100%;
  height: 100%;
  /* max-width: 1200px; */
  margin: 0 auto;
}

.graph-container {
  position: relative;
  width: 100%;
  height: 97vh;
  background-color: white;
  border-radius: 8px;
  /* box-shadow: 0 4px 6px rgba(0,0,0,0.1); */
}

.interact-box {
  position: absolute;
  top: 20px; 
  z-index: 1000;  /* 位于顶层，底层是网络图 */
  background: rgba(255,255,255,0.9);
  padding: 8px;
  border-radius: 4px;
  box-shadow: 0 2px 4px rgba(0,0,0,0.1);
  display: flex; 
  flex-direction: column;
  left: 1%; 
  width: 12%;
}

.select-box {
  position: absolute;
  top: 20px; 
  z-index: 1000;  /* 位于顶层，底层是网络图 */
  background: rgba(255,255,255,0.9);
  padding: 8px;
  border-radius: 4px;
  box-shadow: 0 2px 4px rgba(0,0,0,0.1);
  display: flex; 
  flex-direction: row;
  left: 15%; 
  width: 18%;
}

.reset-btn {
  position: absolute;
  left: 35%; 
  margin-left: -3px; 
  margin-top: 12px; 
  top: 20px; 
  background: rgba(255,255,255,0.9);z-index: 1000;
}


/* :deep(.el-select) {
  width: 240px;
} */
</style>
