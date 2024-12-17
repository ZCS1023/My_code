import os

from utils import *
from reducedimensions import *
from ordercells import *
from select_genes import *
from ordering_genes import *

import warnings
warnings.simplefilter("ignore")
ROOT_DIR = os.path.abspath("/usr/share/dgh/")

''' 5k_cells '''
OUTPUT_DIR = os.path.join(ROOT_DIR, "monocle", "result", "5k_cells")

adata = sc.read_h5ad("/usr/share/dgh/3k_monocle2_test.h5ad")

adata.var['gene_short_name'] = adata.var_names
# data = adata.X
# pd = adata.obs.copy()
# gene_names = adata.var_names
# fData = pd.DataFrame({'gene_short_name': gene_names}, index=gene_names)

adata = estimateSizefactors(adata, round_exprs=True, use_layer=None, inplace=True)
adata = estimateDispersions(adata, min_cells_detected=1)
disp_table = dispersion_table(adata)

# 筛选出符合条件的基因
ordering_genes = disp_table[disp_table['mean_expression'] >= 0.1]
adata = set_ordering_filter(adata, ordering_genes)

# adata = detect_genes(adata, min_expr=None)

# 选取至少在10个细胞中表达的基因
# num_cells_expressed = np.sum(adata.X > 0, axis=0).A1  
# expressed_genes = adata.var_names[num_cells_expressed >= 10]

# diff_test_res = differential_gene_test(adata[:,expressed_genes],full_model_formula_str="~Media")
# ordering_genes = diff_test_res.loc[diff_test_res['qval'] < 0.01].index.tolist()

# def set_ordering_filter(adata, ordering_genes):

#     # 初始化False
#     adata.var['use_for_ordering'] = False
    
#     # 将ordering_genes里的标记为True
#     adata.var.loc[ordering_genes, 'use_for_ordering'] = True
    
#     return adata

# adata = set_ordering_filter(adata, ordering_genes)

adata = reduceDimension(adata)
# adata = orderCells(adata)