import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pysal.lib as ps
from esda.moran import Moran

# 读取10x单细胞RNA-seq数据
data_path = "filtered_gene_bc_matrices/hg19"
adata = sc.read_10x_mtx(data_path, var_names='gene_symbols', cache=True)

# 创建AnnData对象
adata.var_names_make_unique()  # 确保基因名唯一
adata = sc.AnnData(adata.X, obs=adata.obs, var=adata.var)

# 数据预处理：标准化、寻找高变基因、PCA降维、聚类、UMAP可视化
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
sc.pp.scale(adata, max_value=10)

# PCA降维
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
sc.pl.pca_loadings(adata)

# 最近邻搜索并进行聚类
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
sc.tl.leiden(adata, resolution=0.5)

# UMAP降维和可视化
sc.tl.umap(adata)
sc.pl.umap(adata, color=['leiden'])

# 差异基因表达分析
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# 使用PAGA学习细胞轨迹
sc.tl.paga(adata, groups='leiden')
sc.pl.paga_compare(adata)

# 设置伪时间并进行排序
adata.uns['iroot'] = 0  # 手动设置伪时间的起点
sc.tl.dpt(adata)

# 可视化伪时间
sc.pl.umap(adata, color=['dpt_pseudotime'])

# 莫兰指数伪时间差异基因分析
# 生成邻接矩阵
adj_matrix = adata.obsp['connectivities']


# 莫兰指数计算
def calculate_moran_for_gene(adata, gene, adj_matrix):
    # 获取基因表达值
    gene_expression = adata[:, gene].X.toarray().flatten()

    # 使用 pysal 计算莫兰指数
    w = ps.weights.util.full2W(adj_matrix)
    moran = Moran(gene_expression, w)
    return moran.I, moran.p_sim


# 计算所有高变基因的莫兰指数
moran_results = []
for gene in adata.var.highly_variable:
    moran_i, p_value = calculate_moran_for_gene(adata, gene, adj_matrix)
    moran_results.append({'gene': gene, 'moran_I': moran_i, 'p_value': p_value})

# 转换为DataFrame
moran_df = pd.DataFrame(moran_results)

# 筛选显著差异基因（p_value < 0.05）
significant_genes = moran_df[moran_df['p_value'] < 0.05]
print(significant_genes.head())

# 保存差异基因结果到CSV文件
output_file = "pseudotime_significant_genes_moran.csv"
significant_genes.to_csv(output_file)

# 展示UMAP图
plt.show()
