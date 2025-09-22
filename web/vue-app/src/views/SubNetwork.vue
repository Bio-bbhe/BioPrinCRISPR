<template>
  <div>
    <div class="main-container">
      <el-card class="left-container">
        <network
          class="network"
          ref="networkRef"
          :nodes="nodes"
          :edges="edges"
          :options="options"
        ></network>
        <div class="network-info">
          <p style="margin: 0px;margin-bottom: 3px;">Click node to view relevant information on the right</p>
          <p style="margin: 0px;">Double click node to expand or collapse it</p>
        </div>
      </el-card>
      <el-card class="right-container" :body-style="{ padding: '7px' }">
        <div class="clearfix">
          <div slot="header" style="display: flex;flex-direction: row;">
            <el-select v-model="selectedOption" filterable placeholder="" class="select-container" size="small" @change="handleSelectChange">
              <el-option
                v-for="(item, index) in proteinIdMap.get(selectedNodeId) || [selectedNodeId || '']"
                :key="index"
                :label="item.protein_id + '    ( '+ item.source + '  =>  '+ item.target +' )'"
                :value="selectedNodeId + '@' + item.protein_id">
              </el-option>
            </el-select>
            <el-button v-if="selectedOption" type="text" icon="el-icon-document-copy" title="Copy Sequence" style="padding: 0px;" @click="handleCopySequence(false, $event)">Sequence</el-button>
            <el-button v-if="selectedOption" type="text" icon="el-icon-download" title="Download GBK" style="padding: 0px;" @click="handleDownloadGBK(false)">Gbk</el-button>
          </div>
        </div>
        <el-card  class="top-pane" shadow="never" v-loading="svgLoading" >
          <object class="svg-container" :data="svgObjectUrl" type="image/svg+xml"></object>
        </el-card>
        <div class="button-pane">
          <el-card class="left-pane" shadow="never" v-loading="repeatLoading" :body-style="{ padding: '0px'}">
            <div v-if="repeatTable.length > 0">
              <el-table
                id = "repeatTable"
                :data="repeatTable"
                max-height="365"
                border>
                <el-table-column
                  prop="repeat"
                  label="REPEAT"
                  min-width="150">
                </el-table-column>
                <el-table-column
                  prop="spacer"
                  label="SPACER"
                  min-width="170">
                </el-table-column>
              </el-table>
              <div class="repeat-info"><p style="margin: 0px;">Repeat Count: {{repeatInfo.repeat_count}}, Average Repeat Length: {{repeatInfo.avg_repeat_len}}</p><p style="margin: 0px;">Average Spacer Length: {{repeatInfo.avg_spacer_len}}</p></div>
            </div>
          </el-card>
          <el-card class="right-pane" shadow="never" v-loading="pdbLoading">
            <el-button class="download-pdb-btn" v-if="pdbExists" type="primary" plain icon="el-icon-download" title="Download PDB" @click="downloadPDB" ></el-button>
            <div id="molstar-viewer"></div>
          </el-card>
        </div>
      </el-card>
    </div>
    <div class="page-container" v-loading="pageLoading">
      <div class="page-info">
        <span class="page-title">Representative operons</span>
        <el-pagination
          class = "pagination"
          background
          layout="prev, pager, next"
          :current-page.sync="pageNum"
          :total="totalPage"
          :page-size="10"
          @current-change="handleSvgPageNumChange">
        </el-pagination>
      </div>
      <div class="page-pane" v-for="i in Math.ceil(pageProteinIds.length / 2)" :key="i">
        <!-- 使用v-for迭代索引范围[0,1]-->
        <template v-for="j in [0, 1]">
          <el-card v-if="pageProteinIds[(i - 1) * 2 + j]" :key="j" class="svg-page">
            <div class="svg-page-tag" :style="{ marginBottom: j === 0 ? '10px' : '5px' }">
              <el-tag type="info" size="small" style="margin-right: 10px;">
                {{ pageProteinIds[(i - 1) * 2 + j] }}
              </el-tag>
              <!-- 按钮组复用，避免重复声明 -->
              <el-button :id="pageProteinIds[(i - 1) * 2 + j]" type="text" icon="el-icon-document-copy" title="Copy Sequence" style="padding: 0px;" @click="handleCopySequence(true, $event)">
                Sequence
              </el-button>
              <el-button :id="pageProteinIds[(i - 1) * 2 + j]" type="text" icon="el-icon-download" title="Download GBK" style="padding: 0px;" @click="handleDownloadGBK(true, $event)">
                Gbk
              </el-button>
            </div>
            <object class="svg-container" :data="pageSvgs[(i - 1) * 2 + j]" type="image/svg+xml"></object>
          </el-card>
        </template>
      </div>
    </div>
    <!-- <div class="watermark-overlay"></div> -->
  </div>
</template>


<script setup>
import { ref, onMounted, watch } from 'vue';
import { Network, DataSet } from 'vue-vis-network';
import { useRoute } from 'vue-router/composables'
import copy from 'copy-to-clipboard';
import request from '@/utils/request'
import { Message } from 'element-ui';

const route = useRoute();
const nodes = ref(new DataSet([]));  //节点数据（响应式：使用ref处理基本类型数组）
const edges = ref(new DataSet([]));
const networkRef = ref(null);
const showedNodes = new Set();
const openedNodes = new Set();
const showedEdges = new Set();
const proteinIdMap = new Map();
const selectedNodeId = ref(null);
const selectedOption = ref('');
const svgObjectUrl = ref('');
const fontSize = 19;
const svgCache = new Map();
const svgLoading = ref(true);
const isCopied = ref(false);
let viewerInstance = null;
const pdbCache = new Map();
const pdbLoading = ref(true);
const repeatTable = ref([]);
const repeatInfo = ref({});
const repeatLoading = ref(true);
const repeatCache = new Map();
const pageProteinIds = ref([]);
const pageSvgs = ref([]);
let pageNum = ref(1);
const totalPage = ref(0);
const pageLoading = ref(true);
let svgExtraContentText = '';
const pdbExists = ref(false);
let showProteinId = '';
const nextNodeIdsMap = new Map(); // 存储每个节点的下一级节点，用于撤销节点扩展
const nextEdgeIdsMap = new Map(); // 存储每个节点的下一级边，用于撤销节点扩展

// 响应式配置对象（使用reactive处理复杂对象）
const options = {
  interaction: {
    hover: true,
    selectConnectedEdges: false,
    navigationButtons: false,
  },
  physics: {
    enabled: true
  },
  edges: {
    length: 120,
    arrows: {
      to: { 
        enabled: true,
        scaleFactor: 0.8 
      }
    },
    smooth: true,
    color: '#444bf8',
  },
  nodes: {
    shape: 'circle',
    borderWidth: 2,
    fixed: false,
    margin: 8,
    font: {
      color: '#343434', //标签颜色默认黑色
      size: fontSize,
      strokeWidth: 0,
    },
    widthConstraint: { maximum: 70, minimum: 70 }, // 根据新字体大小调整最小宽度
    shapeProperties: { borderDashes: [8, 4] },
  },
};

const molstarOptions = {
  alphafoldView: true,
  bgColor: { r: 255, g: 255, b: 255 },
  // hideControls: true,
  hideCanvasControls: [
    // 'selection',
    // 'animation',
    // 'controlToggle',
    // 'controlInfo',
  ],
  sequencePanel: true,
  // landscape: true,
};

const exampleNodes = { // 示例节点(起始节点：后续节点列表),依次扩展节点
  '4830' : ['4833', '231'],
  '3927' : ['2340'],
  '6039' : ['2307'],
  '888' : ['4125', '672'],
};

watch(selectedNodeId, (newVal) => {
  const options = proteinIdMap.get(newVal) || [];
  if(options.length > 0){
    const firstProteinId = options[0].protein_id;
    selectedOption.value = newVal + '@' + firstProteinId || '';
    showProteinId = firstProteinId;
  }else {
    selectedOption.value = "";
    showProteinId = '';
  }
  fetchSvg(showProteinId);
  renderMolstarPDB(showProteinId);
  fetchRepeat(showProteinId);
  handleSvgPageNumChange(1);
});

const handleSelectChange = async (val) => {
  const proteinId = val.split('@')[1];
  fetchSvg(proteinId);
  renderMolstarPDB(proteinId);
  fetchRepeat(proteinId);
  showProteinId = proteinId;
}

// 获取网络数据
const fetchGraphData = async (nodeId) => {
  try {
    const response = await request.get('/load_data/' + nodeId);
    openedNodes.add(nodeId);
    return response.data.data;
  } catch (error) {
    console.error('数据获取失败:', error);
    Message.error({
      showClose: true,
      message: '数据加载失败',
    });
    return null;
  }
}

const fetchProteinData = async (proteinIds) => {
  try {
    const response = await request.post('/load_data/batch', {
      ids: Array.from(proteinIds)
    });
    // 遍历键值对并添加到Map
    Object.entries(response.data.data).forEach(([key, value]) => {
      proteinIdMap.set(key, value); 
    });
  } catch (error) {
    console.error('数据获取失败:', error);
    Message.error({
      showClose: true,
      message: '数据加载失败',
    });
  }
}

const addNode = (node, opened, x, y) => {
  const lines = node.name.split('\n');
  const borderWidth = opened ? 0 : 2;
  nodes.value.add({
    id: node.id,
    label: lines[0],
    color: getNodeColor(node.number_of_proteins),
    x: x,
    y: y,
    pfam_accession: node.pfam_accession,
    borderWidth
  });
}

// 添加邻居节点
const addNeighborNodes = (nodeId, neighbors, x = 0, y = 0) => {
  const nextNodeIds = new Set();
  neighbors.forEach(neighbor => {
    if (showedNodes.has(neighbor.id)) return;
    showedNodes.add(neighbor.id);
    nextNodeIds.add(neighbor.id);
    addNode(neighbor, false, x, y);
  })
  nextNodeIdsMap.set(nodeId, nextNodeIds);
  return nextNodeIds;
}

// 添加边
const addEdges = (nodeId, edgesData) => {
  const nextEdgeIds = new Set();
  const countEdge = edgesData.length;
  const length = Math.min(300, Math.max(120, countEdge * 10)); //根据边数，动态调整边长度，避免节点太密集
  let offset = 0;
  edgesData.forEach(edge => {
    if (showedEdges.has(edge.id)) return;
    showedEdges.add(edge.id);
    nextEdgeIds.add(edge.id); 
    const newEdge = {
      id: edge.id,
      from: edge.source,
      to: edge.target,
      label: edge.number_of_proteins + '',
      color: getEdgeColor(edge.presence_status),
      length: length,
    }
    if(countEdge > 20){ //边数过多时，调整边类型，避免重叠
      newEdge['length'] = length - 50 + offset
      newEdge['smooth'] = {
          type: 'straightCross',
          enabled: true,
          roundness: 0.8
      }
    }
    offset += 3;
    edges.value.add(newEdge)
  })
  nextEdgeIdsMap.set(nodeId, nextEdgeIds);
}

const deleteNodes = (nodeId) => {
  const nextNodeIds = nextNodeIdsMap.get(nodeId);
  // 删除下一级节点
  if(nextNodeIds) {
    openedNodes.delete(nodeId);
    nextNodeIds.forEach(neighborNodeId => {
      deleteNodes(neighborNodeId);
      nodes.value.remove(neighborNodeId);
      showedNodes.delete(neighborNodeId);
    })
  }
  const nextEdgeIds = nextEdgeIdsMap.get(nodeId);
  // 删除下一级边
  if(nextEdgeIds) {
    openedNodes.delete(nodeId);
    nextEdgeIds.forEach(neighborEdgeId => {
      edges.value.remove(neighborEdgeId);
      showedEdges.delete(neighborEdgeId);
    })
  }
}

const expandNode = async (nodeId, x, y) => {
  nodes.value.update({ //去除虚线边框
      id: nodeId,
      borderWidth: 0,
    })
  const graphData = await fetchGraphData(nodeId);
  const neighborNodeIds = addNeighborNodes(nodeId, graphData.neighbors, x, y);
  addEdges(nodeId, graphData.edges);
  fetchProteinData(neighborNodeIds);
}

const handleDblClick = async (params) => {
  console.log('双击节点:', params);
  if (params.nodes.length > 0) {
    const nodeId = params.nodes[0];
    if(openedNodes.has(nodeId)) { //撤销扩展节点
      nodes.value.update({ //增加虚线边框
        id: nodeId,
        borderWidth: 2,
      })
      deleteNodes(nodeId);
    }else{ //扩展节点
      const canvas = params.pointer.canvas;
      expandNode(nodeId, canvas.x, canvas.y);
    }
  }
};

const fetchSvgExtraContent = async () => {
  const svgExtraContent = await fetch('/svg_extra_content.txt');
  svgExtraContentText = await svgExtraContent.text();
}

const handleExampleNodes = async(nodeId, example) => {
  if(example == 1){
    networkRef.value.network.body.view.scale = 0.9;
    const nodes = exampleNodes[nodeId];
    if(nodes) {
      for(let i = 0; i < nodes.length; i++){
        await expandNode(nodes[i]);
      }
    }
  }
}

// 初始化数据
const initData = async (nodeId, example) => {
  networkRef.value.network.on('click', (params) =>{ //鼠标单击节点或选中节点
    if (params.nodes.length > 0) {
      selectedNodeId.value = params.nodes[0];
    }
  });
  networkRef.value.network.on('doubleClick', (params) =>{ //鼠标双击节点，扩展节点
    if (params.nodes.length > 0) {
      handleDblClick(params);
    }
  });
  console.log(networkRef.value);
  networkRef.value.network.on('afterDrawing', (ctx) =>{
    const network = networkRef.value.network;
    const nodeIndices = network.getPositions(network.body.nodeIndices);
    for (let id in nodeIndices){
      const position = nodeIndices[id];
      const pfam_accession = networkRef.value.getNode(id).pfam_accession;
      ctx.font = `${fontSize - 7}px arial`;
      ctx.fillStyle = '#5e6063';
      ctx.fillText(pfam_accession, position.x - fontSize * 1.2, position.y + fontSize * 1.6);
    }
  });
  const graphData = await fetchGraphData(nodeId);
  if (!graphData) return;
  addNode(graphData.node, true);
  showedNodes.add(nodeId);
  const neighborNodeIds = addNeighborNodes(nodeId, graphData.neighbors);
  const nodeIds = new Set(neighborNodeIds);
  nodeIds.add(nodeId);
  addEdges(nodeId, graphData.edges);
  // networkRef.value.selectNodes([graphData.node.id]); //选中起始节点
  handleExampleNodes(nodeId, example);
  fetchSvgExtraContent();
  await fetchProteinData(nodeIds);
  if(viewerInstance){ //molstarJS加载完时，才渲染pdb文件
    selectedNodeId.value = nodeId;
  }
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

const renderMolstarPDB = async (proteinId) => {
  pdbExists.value = false;
  pdbLoading.value = true;
  const viewerContainer = document.getElementById('molstar-viewer');
  viewerContainer.innerHTML = '';
  if(!proteinId) {
    pdbLoading.value = false;
    return;
  }
  let blob = null;
  if(pdbCache.has(proteinId)) {
    blob = pdbCache.get(proteinId);
  }else{
    const response = await request.get('/pdb', {
      params: {
        id: proteinId,
      }
    });
    if(response.data.status == 'success' && response.data.data){ //pdb文件可能为0kb
      const pdbData = response.data.data;
      blob = new Blob([pdbData]);
      pdbCache.set(proteinId, blob);
    }
  }
  if(blob && viewerInstance) {
    molstarOptions.customData = {
      url: URL.createObjectURL(blob),
      format: 'pdb'
    };
    viewerInstance.render(viewerContainer, molstarOptions);
    pdbExists.value = true;
  }
  pdbLoading.value = false;
}

const downloadPDB = () => {
  if(pdbCache.has(showProteinId)) {
    const blob = pdbCache.get(showProteinId);
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${showProteinId}.pdb`;
    a.click();
    URL.revokeObjectURL(url);
  }
}

const fetchRepeat = async (protein) =>{
  repeatLoading.value = true;
  repeatTable.value = [];
  repeatInfo.value = {};
  if(!protein) {
    repeatLoading.value = false;
    return;
  }
  if(repeatCache.has(protein)){
    repeatTable.value = repeatCache.get(protein).items;
    repeatInfo.value = repeatCache.get(protein).info;
    repeatLoading.value = false;
    return;
  }
  const proteinId = protein.split('_#')[0];
  try {
    const response = await request.get('/repeat', {
      params: {
        id: proteinId,
      }
    });
    if(response.data.status == 'success'){
      const repeatData = response.data.data;
      if(repeatData){
        repeatTable.value = repeatData.items;
        const info = {
          repeat_count: repeatData.repeat_count,
          avg_repeat_len: repeatData.avg_repeat_len,
          avg_spacer_len: repeatData.avg_spacer_len,
        };
        repeatInfo.value = info;
        repeatCache.set(protein, {
          items: repeatData.items,
          info: info
        });
      }
      repeatLoading.value = false;
    }
  } catch (error) {
    console.error('数据获取失败:', error);
    Message.error({
      showClose: true,
      message: '数据加载失败',
    });
  }
}

const getEdgeColor = (status) => {
  if (status == 'Found') {
    return '#F7482A'; 
  }else {
    return '#444bf8';
  }
}

const copySequence = async (protein_id) => {
  try {
    const response = await request.get('/sequence', {
      params: {
        id: protein_id,
      }
    });
    copy(response.data.data);
    isCopied.value = true;
    setTimeout(() => {
      isCopied.value = false;
    }, 5000);
  } catch (error) {
    console.error('数据获取失败:', error);
    Message.error({
      showClose: true,
      message: '数据加载失败',
    });
  }
}

const handleCopySequence = (isPage, e) => {
  if(!isPage) {
    if(!selectedNodeId.value) return;
    copySequence(showProteinId);
  }else{
    copySequence(e.target.parentNode.id);
  }
  e.target.innerHTML = 'Copied';
  setTimeout(() => {
    e.target.innerHTML = 'Sequence';
  }, 5000);
}

const addSvgExtraContent = (svgText, hidden) => {
  //定位 </svg> 标签位置
  const svgCloseTagIndex = svgText.indexOf('</svg>');
  const hiddenText = hidden ? '<div id ="hidden" style="display: none;"></div>' : '';
  //在闭合标签前插入额外内容
  const modifiedSvg = svgText.slice(0, svgCloseTagIndex) + hiddenText + svgExtraContentText + svgText.slice(svgCloseTagIndex);
  return new Blob([modifiedSvg], { type: 'image/svg+xml' });
}

const fetchSvg = async (proteinId) => {
  svgLoading.value = true;
  if(!proteinId) {
    svgObjectUrl.value = URL.createObjectURL(new Blob());
    svgLoading.value = false;
    return;
  }
  let blob;
  if(svgCache.has(proteinId)) {
    blob = svgCache.get(proteinId);
  }else {
    try {
      const response = await request.get('/svg', {
        params: {
          id: proteinId,
        }
      });
      const svgText = response.data.data;
      blob = addSvgExtraContent(svgText, false);
      svgCache.set(proteinId, blob);
    } catch (error) {
      svgObjectUrl.value = URL.createObjectURL(new Blob());
      svgLoading.value = false;
      console.error('数据获取失败:', error);
      Message.error({
      showClose: true,
      message: '数据加载失败',
    });
    }
  }
  svgObjectUrl.value = URL.createObjectURL(blob);
  svgLoading.value = false;
};

const handleSvgPageNumChange = async (val) =>{
  pageNum.value = val;
  pageLoading.value = true;
  try {
    const response = await request.get('/svg/page', {
      params: {
        id: selectedNodeId.value,
        pageNum: pageNum.value
      }
    });
    pageProteinIds.value = response.data.proteinIds;
    totalPage.value = response.data.total;
    pageSvgs.value = [];
    for(let i = 0; i < response.data.svgs.length; i++){
      const svgText = response.data.svgs[i];
      const blob = addSvgExtraContent(svgText, true);
      pageSvgs.value.push(URL.createObjectURL(blob));
    }
    pageLoading.value = false;
  } catch (error) {
    console.error('数据获取失败:', error);
    Message.error({
      showClose: true,
      message: '数据加载失败',
    });
  }
}

const initMolstarJS = () => {
  const script = document.createElement('script');
  script.src = '/js/pdbe-molstar-plugin.js';
  script.onload = () => {
    // JS 加载完成后初始化
    if (window.PDBeMolstarPlugin) {
      viewerInstance = new window.PDBeMolstarPlugin();

      if(proteinIdMap.has(route.query.id)){ //初始节点的proteins准备好时，才渲染pdb文件
        selectedNodeId.value = route.query.id;
      }
    }
  };
  document.body.appendChild(script);

  const cssUrl = '/css/pdbe-molstar-light.css'; // 可动态生成路径
  const link = document.createElement('link');
  link.rel = 'stylesheet';
  link.type = 'text/css';
  link.href = cssUrl;
  link.onerror = () => console.error('CSS 加载失败');
  document.head.appendChild(link);
}

const downloadGBK = async (protein_id) => {
  try {
    const response = await request.get('/gbk', {
      params: {
        id: protein_id,
      }
    });
    if(response.data.status == 'success'){
      const blob = new Blob([response.data.data]);
      const a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = protein_id + '.gb';
      a.click();
      URL.revokeObjectURL(a.href);
    }
  } catch (error) {
    console.error('数据获取失败:', error);
    Message.error({
      showClose: true,
      message: '数据加载失败',
    });
  }
}

const handleDownloadGBK = (isPage, e) => {
  if(isPage) {
    downloadGBK(e.target.parentNode.id);
  }else{
    if(!selectedOption.value) return;
    downloadGBK(showProteinId);
  }
}


// 生命周期钩子
onMounted(() => {
  document.title = 'Network (' + route.query.label + ')';
  initMolstarJS();
  initData(route.query.id, route.query.example);
})
</script>

<style>
.main-container {
  display: flex;
  flex-direction: row;
  height: 96vh;
  gap: 10px;
  padding: 5px;
}

.network-info {
  position: absolute; 
  bottom: 3.5vh; 
  left: 25px; 
  z-index:1;
  font-size: 15px;
  color: #b4b3b3;
}

/* .watermark-overlay {
  position: fixed;
  top: 0;
  left: 0;
  width: 100vw;
  height: 100vh;
  background-image: url('data:image/svg+xml;utf8,<svg xmlns="http://www.w3.org/2000/svg" width="200" height="100"><text x="50%" y="50%" fill="rgba(0,0,0,0.1)" font-size="10" text-anchor="middle" dominant-baseline="middle" transform="rotate(-30, 100, 50)">华大BGI</text></svg>');
  background-repeat: repeat;
  background-size: 350px 250px; 
  pointer-events: none;
  z-index: 9999;
  opacity: 0.5;
} */

/* .el-card {
  box-shadow: none !important;
} */

.left-container {
  flex: 1;
  height: 100%;
}

.network {
  position: relative;
  height: 94vh; /* 继承容器高度 */
  margin: -15px;
}

.right-container {
  flex: 1;
  display: flex;
  flex-direction: column;
  /* gap: 10px; */
  height: 100%;
}

.select-container {
  width: 70%;
  margin-right: 10px;
}

.clearfix {
  flex: 1;
  height: 3%;
  padding-bottom: 20px;
  margin-top: -3px;
  width: 100%;
  min-width: 763px;
  padding-top: 5px;
}

.top-pane {
  flex: 1;
  height: 33vh;
  max-height: 50vh;
  overflow-y: auto;
  box-sizing: border-box;
  width: 100%;
}

.svg-container {
  margin-right: -40px;
  margin-left: -10px;
  margin-bottom: -20px;
  margin-top: -10px;
  /* transform: scale(1.1);
  transform-origin: 0 0; */
}


.button-pane {
  flex: 1;
  height: 55vh;
  display: flex;
  margin-top: 5px;
  width: 100%;
}

#repeatTable {
  border: 1px solid #e4e7ed;
}

.left-pane{
  width: 45%;
  box-sizing: border-box;
  position: relative;
  margin-right: 5px;
}

.right-pane {
  width: 55%;
  box-sizing: border-box;
  position: relative;
}


#molstar-viewer {
  width: 100%;
  height: 100px;
}

.repeat-info {
  margin: 5px;
  color: #696a6a;
  /* font-weight: bold; */
  font-size: 15px;
}

.download-pdb-btn {
  width : 32px; 
  height: 32px; 
  padding: 0px; 
  color: #93621b; 
  border: 0px; 
  background-color: rgb(233 230 224 / 36%); 
  position: absolute; 
  top: 218px; 
  right: 10px; 
  z-index:1;
}

::v-deep .download-pdb-btn:hover {
  color: #734602; 
  border: 1px solid #c8c4be; 
  background-color: rgb(233 230 224 / 36%); 
}

::v-deep .download-pdb-btn:focus {
  color: #734602; 
  border: 1px solid #c8c4be; 
  background-color: rgb(233 230 224 / 36%); 
}

/* 解决molstar全屏时，表格边框显示异常的问题 */
.el-table::before {
    height: 0px;
}
.el-table--border::after, .el-table--group::after {
    width: 0px;
}
.page-container {
  flex: 1;
  height: 100%;
  display: flex;
  flex-direction: column;
  padding: 5px;
  margin-top: 20px;
}

.page-info {
  display: flex;
  flex-direction: row;
  margin-bottom: 10px;
  margin-left: 5px;
}

.page-title {
  font-size: large;
  font-weight: bold;
  padding-top: 3px; 
  margin-right: 5px;
}

.pagination {
  padding: 0px;
}

.page-pane {
  flex: 1;
  display: flex;
  flex-direction: row;
  gap: 10px;
}

.svg-page {
  height: 22.8vh;
  max-height: 22.8vh;
  overflow-y: auto;
  box-sizing: border-box;
  width: 50%;
  margin-bottom: 10px;
}

.el-tag.el-tag--info {
    background-color: #f4f4f5;
    border-color: #e9e9eb;
    color: #343538;
}

.svg-page-tag {
  display: flex;
  flex-direction: row;
}

</style>
