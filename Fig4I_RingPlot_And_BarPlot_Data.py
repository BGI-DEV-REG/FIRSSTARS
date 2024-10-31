import scanpy as sc
import numpy as np
import pandas as pd

FB_l = ["Cxcl2", "Ccl7", "Il1b", "Cxcl3", "Ccl3"]
Mono_r = ["Cxcr2", "Ccr2", "Il1rap", "Il1r1", "Ccr5"]

FB_myo = ["FB_fascia_1_Adult", "FB_inflammatory_Adult", "FB_oxidative_stress_Adult", "FB_myo_Adult", "FB_fascia_EP", "FB_inflammatory_EP", "FB_oxidative_stress_EP"]
Mono_celltypes = ["Mono_Ly6c"]


adata = sc.read("E:/CodePrograms/commot/Adult_1D_C01830C1D1_LR_commot_0702.h5ad")
adata_df = pd.read_csv("E:/CodePrograms/commot/Adult_1D_C01830C1D1_detailed_celltype.csv.gz")


adata.obs['CellID'] = adata.obs.index
# 用CellID匹配并替换celltype
adata.obs['celltype'] = adata.obs['CellID'].map(adata_df.set_index('X')['celltype'])
selected_cellid = adata_df["X"]
adata = adata[adata.obs['CellID'].isin(selected_cellid)]


adata_FB = adata[adata.obs["celltype"].isin(FB_myo)]
adata_Mono = adata[adata.obs["celltype"].isin(Mono_celltypes)]
print(upper_id)

for fb_gene, mono_gene in zip(FB_l, Mono_r):
  ct_LR = "-".join([fb_gene, mono_gene])
if f"s-{ct_LR}" in adata.obsm["commot-CellChat-sum-sender"].columns and f"r-{ct_LR}" in adata.obsm["commot-CellChat-sum-receiver"].columns:
  FB_s_df = adata_FB.obsm["commot-CellChat-sum-sender"][f"s-{ct_LR}"]
FB_s_id = FB_s_df > 0
print(ct_LR)
print("FB_prop:")
sum(FB_s_id) / len(FB_s_df)
FB_s_id_FT = FB_s_id.to_numpy().flatten()
FB_s_exp = FB_s_df[FB_s_id].to_numpy().flatten()
FB_sum_exp = sum(FB_s_exp)
FB_ave_exp = FB_sum_exp / len(FB_s_df)
print("FB_ave_exp:")
print(FB_ave_exp)

Mono_r_df = adata_Mono.obsm["commot-CellChat-sum-receiver"][f"r-{ct_LR}"]
Mono_r_id = Mono_r_df > 0
print("Mono_prop:")
sum(Mono_r_id) / len(Mono_r_df)
Mono_r_id_FT = Mono_r_id.to_numpy().flatten()
Mono_r_exp = Mono_r_df[Mono_r_id].to_numpy().flatten()
Mono_sum_exp = sum(Mono_r_exp)
Mono_ave_exp = Mono_sum_exp / len(Mono_r_df)
print("Mono_ave_exp:")
print(Mono_ave_exp)
