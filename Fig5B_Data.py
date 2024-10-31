import scanpy as sc
import pandas as pd
import numpy as np
from itertools import product
import matplotlib.pyplot as plt

adata = sc.read_h5ad("E:/CodePrograms/Immu/Immu_0524/20240510_scvi.h5ad")

Mono_celltype = ["Mono_classical","Macr_remodeling","Macr_neovascularization","Mono_noclassical","Mono_intermediate"]
Immu_df = pd.read_csv("E:/CodePrograms/Immu/Immu_0703/Immu_meta_0703.csv")
Mono_df = Immu_df[Immu_df["celltype"].isin(Mono_celltype) & Immu_df["Age"].isin(["P5","Adult"])]
adata_Immu = adata[adata.obs.index.isin(Mono_df["X"])]
adata_Immu.obs["celltype"] = adata_Immu.obs.index.map(Mono_df.set_index('CellID')['celltype'])

sc.pp.neighbors(adata_Immu, use_rep='X_scVI')

MIN_DISTS = np.arange(0.1,10.1,0.1)
SPREADS = np.arange(0.5,5.5,0.5)

for (i, min_dist), (j, spread) in product(enumerate(MIN_DISTS), enumerate(SPREADS)):
    param_str = " ".join(["min_dist =", str(min_dist), "and spread =", str(spread)])
    sc.tl.umap(adata_Immu, min_dist=min_dist, spread=spread)
    sc.pl.umap(adata_Immu, color = "celltype", legend_loc='on data', 
              title=param_str, legend_fontsize=10, legend_fontoutline=2, frameon=False)
    plt.savefig(f"./set_seed/Immu_seed_{str(min_dist)}_{str(spread)}_celltype_with_names.png")
    plt.close()
    
    umap_data = adata_Immu.obsm['X_umap']
    row_name = pd.DataFrame({'RowNames': adata_Immu.obs.index})
    umap_table = pd.concat([row_name, pd.DataFrame(umap_data, columns=['UMAP1', 'UMAP2'])], axis=1)
    umap_table.set_index('RowNames', inplace=True)
    umap_table = umap_table[['UMAP1', 'UMAP2']]
    obs_data = pd.DataFrame(adata_Immu.obs)
    merged_df = pd.concat([obs_data, umap_table], axis = 1)
    merged_df.to_csv("".join(["./set_seed/All_Immu_Seed_",str(min_dist),"_",str(spread),"_obs_umap.csv"]), index= True)
