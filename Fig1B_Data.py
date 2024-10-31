import argparse
import os
from datetime import datetime
import pandas as pd
import anndata
from functools import reduce
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import scanpy.external as sce
from harmony import harmonize


parser = argparse.ArgumentParser(description="scRNA-seq data integration")
parser.add_argument('-I', '--input', help='input path of qc result directory')
parser.add_argument('-D', '--dataframe', help='Dataframe of QC result ')
parser.add_argument('-N', '--name', help='The prefix of integrated object')
parser.add_argument('-O', '--out', help='Output directory')


args = parser.parse_args()
print(args)
print(datetime.now())
os.makedirs(args.out, exist_ok=True)
os.chdir(args.out)


df = pd.read_csv(args.dataframe)
select_rows = df[(df['Technology'] == 'RNA') & (df['Species'] == 'r') & (df['QCpassed'] == 'Pass')]
AnalysisName = select_rows["AnalysisName"]


obj_list = []
for i in AnalysisName:
    file_path = os.path.join(args.input, i, f"{i}.h5ad")
    obj = anndata.read_h5ad(file_path)
    obj.obs['AnalysisName'] = i
    Aname_split  = i.split("_")
    obj.obs["Age"] = Aname_split[2]
    obj.obs["TimePostInjury"] = Aname_split[3]
    obj.obs["Age_TimePostInjury"] = Aname_split[2] + "_" + Aname_split[3]
    obj.obs["SampleIndex"] = Aname_split[4]
    obj_list.append(obj)


merge_obj = anndata.concat(obj_list)


s_genes = pd.read_csv('/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/03.RNA/02.Integrate/30.Rat.h5ad/s_genes.csv')['x'].tolist()
s_genes = [gene.title() for gene in s_genes]
g2m_genes = pd.read_csv('/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/03.RNA/02.Integrate/30.Rat.h5ad/g2m_genes.csv')['x'].tolist()
g2m_genes = [gene.title() for gene in g2m_genes]


sc.pp.normalize_total(merge_obj, target_sum=1e4)
sc.pp.log1p(merge_obj)
sc.pp.filter_genes_dispersion(merge_obj, n_top_genes=3000)
sc.tl.score_genes_cell_cycle(merge_obj, s_genes=s_genes, g2m_genes=g2m_genes)
sc.pp.regress_out(merge_obj, keys=['S_score', 'G2M_score'])
sc.pp.pca(merge_obj)
sce.pp.harmony_integrate(merge_obj, ["AnalysisName", "phase"])
sc.pp.neighbors(merge_obj, use_rep = "X_pca_harmony")
sc.tl.louvain(merge_obj,resolution = 1.5, key_added = 'louvain1')
sc.tl.louvain(merge_obj,resolution = 2.0, key_added = 'louvain2')
sc.tl.louvain(merge_obj,resolution = 3.0, key_added = 'louvain3')
sc.tl.umap(merge_obj)

merge_obj.write_h5ad(f"{args.name}_louvain.h5ad")
merge_obj.obs.to_csv(f"{args.name}_obs.csv", index= True)

umap_data = merge_obj.obsm['X_umap']
row_name = pd.DataFrame({'RowNames': merge_obj.obs.index})
umap_table = pd.concat([row_name, pd.DataFrame(umap_data, columns=['UMAP1', 'UMAP2'])], axis=1)
umap_table.set_index('RowNames', inplace=True)
umap_table = umap_table[['UMAP1', 'UMAP2']]
umap_table.to_csv(f"{args.name}_umap.csv", index= True)

obs_data = pd.DataFrame(merge_obj.obs)
merged_df = pd.concat([obs_data, umap_table], axis = 1)
merged_df.to_csv(f"{args.name}_obs_umap.csv", index= True)
