import scanpy as sc
import pysal as ps
import numpy as np
import matplotlib.pyplot as plt

# 加载示例单细胞 RNA-seq 数据
adata = sc.datasets.pbmc3k()  # PBMC 3k 单细胞数据
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
adata = adata[:, adata.var['highly_variable']]

# 计算 PCA 并进行拟时序分析
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.dpt(adata)  # 拟时序分析

# 差异基因分析
sc.tl.rank_genes_groups(adata, 'dpt_pseudotime', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)

from pysal.lib import weights
from esda.moran import Moran

# 获得基因表达值
gene_expression = adata.raw.X[:, 0]  # 第一个基因表达的值
coords = adata.obsm['X_umap']  # 使用 UMAP 或空间坐标

# 构建空间加权矩阵
w = weights.KNN.from_array(coords, k=5)  # 最近邻加权矩阵
w.transform = 'r'

# 计算莫兰指数
moran = Moran(gene_expression, w)
print(f"莫兰指数: {moran.I}, p值: {moran.p_sim}")

sc.pl.dpt_groups_pseudotime(adata)
sc.pl.dpt_timeseries(adata, gene_names=['GeneA', 'GeneB', 'GeneC'])
