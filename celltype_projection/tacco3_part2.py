# /jdfssz1/ST_SUPERCELLS/P21Z10200N0171/USER/wangshuai3/miniconda/miniconda3/envs/Tacco/bin/python3.9

import os
import sys

import anndata as ad
import argparse
import pandas as pd
import numpy as np
import tacco as tc
import scanpy as sc


def parse_modifications(modifications_str):
    """
    Parse a string of modifications in the format 'key1=value1,key2=value2' into a dictionary.
    """
    modifications = {}
    for item in modifications_str.split(','):
        key, value = item.split('=')
        modifications[key] = float(value)
    return modifications


parser = argparse.ArgumentParser(
                    prog='Iter Tacco 3',
                    description='Transferring annotations from single cell to spatial transcriptomics data',
                    )
parser.add_argument('-R', '--reference')
parser.add_argument('-S', '--spatial')
parser.add_argument('-P', '--save_path', type=str, default='result')
parser.add_argument('-A', '--annotation', type=str)
parser.add_argument('--meta_sc', type=str)
parser.add_argument('--meta_sp', type=str)
parser.add_argument('--celltype_sp', type=str)
parser.add_argument('--selected', nargs='+')
parser.add_argument('--markers', type=str)
parser.add_argument('--min_cell', type=int, default=20)
parser.add_argument('--modifications', type=parse_modifications)
args = parser.parse_args()
print(args)

if not os.path.exists(args.save_path):
    os.mkdir(args.save_path)
os.chdir(args.save_path)


def cluster_small_multiples(adata, clust_key, size=10, frameon=False, **kwargs):
    tmp = adata.copy()
    for i, clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype('category')
        tmp.uns[clust+'_colors'] = ['#d3d3d3', adata.uns[clust_key+'_colors'][i]]

    sc.pl.scatter(tmp, groups=tmp.obs[clust].cat.categories[1:].values, color=adata.obs[clust_key].cat.categories.tolist(), size=size, frameon=frameon, **kwargs)


def markers_plot(adata, clust_key, df, **kwargs):
    adata = adata.copy()
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)
    markers = {i['celltype']: [j for j in i['marker'].split(',') if j in adata.var_names] for (_, i) in df.iterrows()}
    sc.pl.dotplot(adata, markers, groupby=clust_key, **kwargs)



def update_celltype_proportions2(df, modifications, ct='celltype'):
    # 计算每个粗注释的细胞比例
    celltype_proportions = df[ct].value_counts(normalize=True)

    if modifications is None:
        return celltype_proportions

    # 计算未被修改的细胞类型的当前总比例
    current_remaining_proportion = celltype_proportions.drop(modifications.keys()).sum()
    # 计算修改后剩余的比例
    remaining_proportion = 1 - sum(modifications.values())

    # 更新未被修改的细胞类型的比例
    celltype_proportions.update(
        celltype_proportions.drop(modifications.keys()) * remaining_proportion / current_remaining_proportion)

    # 更新指定的粗注释细胞比例
    celltype_proportions.update(pd.Series(modifications))

    return celltype_proportions


reference = ad.read(args.reference)
meta_sc = pd.read_csv(args.meta_sc, index_col=0)
cells = reference.obs_names.intersection(meta_sc.index)
reference = reference[cells, :]
reference.obs[args.annotation] = meta_sc.loc[reference.obs.index, args.annotation]


spatial_raw = ad.read(args.spatial)
meta_sp = pd.read_csv(args.meta_sp, index_col=0)
spatial_raw.obs[args.celltype_sp] = meta_sp.loc[spatial_raw.obs.index, args.celltype_sp]
spatial_raw.X = spatial_raw.layers['counts']
spatial_raw.obsm['X_spatial'] = spatial_raw.obsm['spatial']
spatial = spatial_raw[spatial_raw.obs[args.celltype_sp].isin(args.selected), :].copy()

table = reference.obs[args.annotation].value_counts()
removed = table[table < args.min_cell].index.tolist()
reference = reference[~reference.obs[args.annotation].isin(removed), :]
print(f"Cell types with fewer than {args.min_cell} cells removed:\n{table[table < args.min_cell]}")

raw_prob = reference.obs[args.annotation].value_counts()/reference.obs[args.annotation].size
print("reference annotation_prior:")
print(raw_prob)
print()

# mean_prob = (raw_prob > 0).astype('float')
# mean_prob = mean_prob/mean_prob.sum()
celltype_proportions = update_celltype_proportions2(reference.obs, args.modifications, ct=args.annotation)
print("modified annotation_prior:")
print(celltype_proportions)
print()

# prob = raw_prob.copy()  if args.raw_prob else mean_prob.copy()
prob = celltype_proportions.copy()

markers_df = pd.read_csv(args.markers, index_col=0)
markers_plot(reference, args.annotation, markers_df, show=False, save=f"_reference.png")

print(f"input annotation_prior:")
print(prob)
print()

prob.to_csv(f'anno_prior.csv')

tc.tl.annotate(spatial, reference, args.annotation,
               result_key='pred_celltype', annotation_prior=prob, verbose=False)

prediction = spatial.obsm['pred_celltype'].copy()
prediction['pred_celltype'] = prediction.idxmax(1)
prediction['max_score'] = prediction.iloc[:, :-1].max(1)
pred_prob = prediction['pred_celltype'].value_counts(normalize=True)

mse = ((pred_prob - prob) ** 2).mean()
print(f"mse {mse:.5f}")
print(pred_prob)
print()
pred_prob.to_csv(f'pred_prob.csv')

spatial.obs['pred_celltype'] = prediction['pred_celltype']
spatial_raw.obs['pred_celltype'] = spatial_raw.obs[args.celltype_sp]
spatial_raw.obs.loc[spatial.obs_names, 'pred_celltype'] = prediction['pred_celltype']
prediction.to_csv(f'pred.csv.gz')

sc.pl.scatter(spatial, color='pred_celltype', basis="spatial", show=False, save=f"pred_celltype.png")
cluster_small_multiples(spatial, 'pred_celltype', basis="spatial", show=False, save=f"pred_celltype_split.png")
markers_plot(spatial, 'pred_celltype', markers_df, show=False, save=f"spatial.png")
sc.pl.scatter(spatial_raw, color='pred_celltype', basis="spatial", show=False, save=f"pred_celltype_all.png")
cluster_small_multiples(spatial_raw, 'pred_celltype', basis="spatial", show=False, save=f"pred_celltype2_all_split.png")



print("Done!")
