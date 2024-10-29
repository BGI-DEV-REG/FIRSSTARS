import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import squidpy as sq
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import hypergeom

parser = argparse.ArgumentParser()
parser.add_argument('-I', '--input')
parser.add_argument('-M', '--meta')
parser.add_argument('-F', '--save_filename')
parser.add_argument('-K', '--k', type=int, default=20)

args = parser.parse_args()
print(args)



def count_adjacent_categories(A, categories):
    # 获取类别的唯一值和它们的索引
    unique_categories = np.unique(categories)
    
    # 转换类别为一个one-hot编码的稀疏矩阵
    # 每一行对应一个节点，每一列对应一个类别
    n_nodes = len(categories)
    n_categories = len(unique_categories)
    rows = np.arange(n_nodes)
    cols = [np.where(unique_categories == category)[0][0] for category in categories]
    data = np.ones(n_nodes)
    category_matrix = csr_matrix((data, (rows, cols)), shape=(n_nodes, n_categories))
    
    # 计算邻接节点的类别计数
    adjacent_category_matrix = A.dot(category_matrix)
    
    category_counts_df = pd.DataFrame(0, index=unique_categories, columns=unique_categories)
    # 对于每个类别，统计相邻节点的类别计数
    for category in unique_categories:
        category_index = np.where(unique_categories == category)[0][0]
        # 获取该类别所有节点的相邻类别计数
        category_counts_df.loc[category] = adjacent_category_matrix[rows[categories == category]].sum(axis=0)
    
    return category_counts_df


def hypergeom_test(category_counts_df):
    M = category_counts_df.to_numpy().sum()
    results_df = pd.DataFrame(columns=['From','To','Counts','Log2FC','Pvalue'])
    
    for i, row in category_counts_df.iterrows():
        N = row.sum()
        for j in category_counts_df.columns:
            n = category_counts_df[j].sum()
            x = row[j]
            foldchange = np.log2((x / N) / (n / M))
            p_value = 1 - hypergeom(M, n, N).cdf(x - 1)
            results_df = results_df._append({'From': i, 'To': j, 'Counts': x, 'Log2FC': foldchange, 'Pvalue': p_value, 'Percentage': x/N}, ignore_index=True)
    
    return results_df


adata = sc.read_h5ad(args.input)
meta = pd.read_csv(args.meta, index_col=0)
adata = adata[meta.index,:]
adata.obs = meta


sq.gr.spatial_neighbors(adata, n_neighs=args.k)
category_counts = count_adjacent_categories(adata.obsp['spatial_connectivities'], adata.obs['pred_celltype'])

hypergeom_results = hypergeom_test(category_counts)
hypergeom_results.to_csv(args.save_filename)

