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
parser.add_argument('-B', '--annotation2', type=str)
parser.add_argument('--meta', type=str)
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


def update_celltype_proportions(df, modifications, ct_lowres='celltype_lowres', ct_hires="celltype_hires"):
    # 计算每个粗注释的细胞比例
    celltype_proportions = df[ct_lowres].value_counts(normalize=True)

    # 对每个粗注释分组，然后计算每个粗注释中其细分群的细胞比例
    hires_proportions = df.groupby(ct_lowres)[ct_hires].value_counts(normalize=True)
    grouped_proportions = hires_proportions.mul(celltype_proportions, level=0)

    if modifications is None:
        return celltype_proportions, grouped_proportions

    # 计算未被修改的细胞类型的当前总比例
    current_remaining_proportion = celltype_proportions.drop(modifications.keys()).sum()
    # 计算修改后剩余的比例
    remaining_proportion = 1 - sum(modifications.values())

    # 更新未被修改的细胞类型的比例
    celltype_proportions.update(
        celltype_proportions.drop(modifications.keys()) * remaining_proportion / current_remaining_proportion)

    # 更新指定的粗注释细胞比例
    celltype_proportions.update(pd.Series(modifications))

    # 调整每个细分群的细胞比例
    grouped_proportions = hires_proportions.mul(celltype_proportions, level=0)

    return celltype_proportions, grouped_proportions


reference = ad.read(args.reference)
if args.meta:
    meta = pd.read_csv(args.meta, index_col=0)
    reference.obs[args.annotation] = meta.loc[reference.obs.index, args.annotation]
    reference.obs[args.annotation2] = meta.loc[reference.obs.index, args.annotation2]

spatial = ad.read(args.spatial)
spatial.X = spatial.layers['counts']
spatial.obsm['X_spatial'] = spatial.obsm['spatial']

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
celltype_proportions, grouped_proportions = update_celltype_proportions(reference.obs,
                                                                        args.modifications,
                                                                        ct_lowres=args.annotation,
                                                                        ct_hires=args.annotation2)
print("modified annotation_prior:")
print(celltype_proportions)
print()
print(grouped_proportions)
print()

# prob = raw_prob.copy()  if args.raw_prob else mean_prob.copy()
prob = grouped_proportions.reset_index(level=0, drop=True)
ct_map = grouped_proportions.reset_index(level=0)
ct_map = ct_map.iloc[:, :-1]

markers_df = pd.read_csv(args.markers, index_col=0)
markers_plot(reference, args.annotation, markers_df, show=False, save=f"_reference.png")
markers_plot(reference, args.annotation2, markers_df, show=False, save=f"_reference2.png")

print(f"input annotation_prior:")
print(prob)
print()

celltype_proportions.to_csv(f'anno_prior_1.csv')
prob.to_csv(f'anno_prior.csv')

tc.tl.annotate(spatial, reference, args.annotation2,
               result_key='pred_celltype2', annotation_prior=prob, verbose=False)

prediction = spatial.obsm['pred_celltype2'].copy()
prediction['pred_celltype2'] = prediction.idxmax(1)
prediction['max_score'] = prediction.iloc[:, :-1].max(1)
prediction['pred_celltype'] = ct_map.loc[prediction['pred_celltype2'], args.annotation].values
pred_prob = prediction['pred_celltype2'].value_counts(normalize=True)
pred_prob_1 = prediction['pred_celltype'].value_counts(normalize=True)

mse = ((pred_prob - prob) ** 2).mean()
print(f"mse {mse:.5f}")
print(pred_prob)
print()
print(pred_prob_1)
print()
pred_prob_1.to_csv(f'pred_prob_1.csv')

spatial.obs['pred_celltype'] = prediction['pred_celltype']
spatial.obs['pred_celltype2'] = prediction['pred_celltype2']
prediction.to_csv(f'pred.csv.gz')

sc.pl.scatter(spatial, color='pred_celltype', basis="spatial", show=False, save=f"pred_celltype.png")
cluster_small_multiples(spatial, 'pred_celltype', basis="spatial", show=False, save=f"pred_celltype_split.png")
markers_plot(spatial, 'pred_celltype', markers_df, show=False, save=f"spatial.png")
sc.pl.scatter(spatial, color='pred_celltype2', basis="spatial", show=False, save=f"pred_celltype2.png")
cluster_small_multiples(spatial, 'pred_celltype2', basis="spatial", show=False, save=f"pred_celltype2_split.png")
markers_plot(spatial, 'pred_celltype2', markers_df, show=False, save=f"_spatial_2.png")


print("Done!")
