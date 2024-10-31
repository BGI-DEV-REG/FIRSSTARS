import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

Mono_l = ["Col4a1", "Col6a1", "Col1a1", "Col1a2"]
FB_r = ['Sdc1', "Sdc1", "Sdc1", "Sdc1"]

FB_myo = ["FB_fascia_1_Adult", "FB_inflammatory_Adult", "FB_oxidative_stress_Adult", "FB_myo_Adult", "FB_fascia_EP", "FB_inflammatory_EP", "FB_oxidative_stress_EP"]
Mono_celltypes = ["Macr_Cd93", "Macr_Mmp14"]

adata = sc.read("Adult_5D_C01830E6F6_LR_commot_0702.h5ad")
adata_df = pd.read_csv("Adult_5D_C01830E6F6_detailed_celltype.csv.gz")

if "spatial_leiden_e30_s8_colors" in adata.uns:
    del adata.uns["spatial_leiden_e30_s8_colors"]

adata.obs['CellID'] = adata.obs.index
adata.obs['celltype'] = adata.obs['CellID'].map(adata_df.set_index('X')['celltype'])
selected_cellid = adata_df["X"]
adata = adata[adata.obs['CellID'].isin(selected_cellid)]

umap_coords = adata.obsm['spatial']
x_range = np.max(umap_coords[:, 0]) - np.min(umap_coords[:, 0])
y_range = np.max(umap_coords[:, 1]) - np.min(umap_coords[:, 1])
aspect_ratio = x_range / y_range

adata_FB = adata[adata.obs["celltype"].isin(FB_myo)]
adata_Mono = adata[adata.obs["celltype"].isin(Mono_celltypes)]

combined_FB_s_df = pd.DataFrame(0, index=adata_FB.obs.index, columns=[None])
combined_Mono_r_df = pd.DataFrame(0, index=adata_Mono.obs.index, columns=[None])

for fb_gene, mono_gene in zip(FB_r, Mono_l):
    ct_LR = "-".join([mono_gene, fb_gene])
    
    if f"s-{ct_LR}" in adata.obsm["commot-CellChat-sum-sender"].columns and f"r-{ct_LR}" in adata.obsm["commot-CellChat-sum-receiver"].columns:
        FB_s_df = adata_FB.obsm["commot-CellChat-sum-sender"][f"s-{ct_LR}"]
        FB_s_df = FB_s_df.to_frame(name=f"s-{ct_LR}")
        combined_FB_s_df = combined_FB_s_df.add(FB_s_df, fill_value=0)
        Mono_r_df = adata_Mono.obsm["commot-CellChat-sum-receiver"][f"r-{ct_LR}"]
        Mono_r_df = Mono_r_df.to_frame(name=f"r-{ct_LR}")
        combined_Mono_r_df = combined_Mono_r_df.add(Mono_r_df, fill_value=0)

combined_FB_s_df['row_sum'] = combined_FB_s_df.sum(axis=1)
combined_Mono_r_df['row_sum'] = combined_Mono_r_df.sum(axis=1)


FB_s_id = combined_FB_s_df["row_sum"] > 0
FB_s_id_FT = FB_s_id.to_numpy().flatten()
FB_s_exp = combined_FB_s_df[FB_s_id].to_numpy().flatten()
FB_s_umap = adata_FB.obsm["spatial"][FB_s_id_FT]

Mono_r_id = combined_Mono_r_df["row_sum"] > 0
Mono_r_id_FT = Mono_r_id.to_numpy().flatten()
Mono_r_exp = combined_Mono_r_df[Mono_r_id].to_numpy().flatten()
Mono_r_umap = adata_Mono.obsm["spatial"][Mono_r_id_FT]

    exp_combined_coordinates = np.vstack((FB_s_umap, Mono_r_umap))
    exp_combined_set = set(map(tuple, exp_combined_coordinates))
        
    non_combined_list = [coord for coord in umap_coords if tuple(coord) not in exp_combined_set]
    non_combined_coordinates = np.array(non_combined_list)
    print(f"Number of coordinates not in combined set: {non_combined_coordinates.shape[0]}")

for pt_size in [10,20,30,40,50,80,90,100]:

                plt.style.use('dark_background')                                                                                                                                                                                                                  
                fig, ax = plt.subplots(figsize=(30, 30 / aspect_ratio))

                scatter3 = ax.scatter(
                    non_combined_coordinates[:, 0], non_combined_coordinates[:, 1],
                    c='#484848', label='Coordinates not in combined set', s=3
                )

                scatter1 = ax.scatter(
                    FB_s_umap[:, 0], FB_s_umap[:, 1],
                    c="#fdd28a",
                    cmap=custom_cmap_fb, label=fb_gene, s=pt_size, edgecolors='none'
                )

                scatter2 = ax.scatter(
                    Mono_r_umap[:, 0], Mono_r_umap[:, 1],
                    c="#ff0038",
                    cmap=custom_cmap_mono, label=mono_gene, s=pt_size, edgecolors='none' 
                )

                plt.savefig(f"./Adult_5D_Macr/Col_and_Sdc1_commot_colorbar_{pt_size}.png", dpi=300)
                plt.close()
