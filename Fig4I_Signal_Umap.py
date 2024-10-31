import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

FB_l = ["Cxcl2", "Ccl7", "Il1b", "Cxcl3", "Ccl3"]
Mono_r = ["Cxcr2", "Ccr2", "Il1rap", "Il1r1", "Ccr5"]

FB_myo = ["FB_fascia_1_Adult", "FB_inflammatory_Adult", "FB_oxidative_stress_Adult", "FB_myo_Adult", "FB_fascia_EP", "FB_inflammatory_EP", "FB_oxidative_stress_EP"]
Mono_celltypes = ["Mono_Ly6c"]

Chip_data = pd.read_csv("/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/03.RNA/02.Integrate/50.Rat_All/Immu/05.st_featureplot/chipID_Early.csv")

for index, (lower_id, sct, upper_id) in Chip_data[['lower_ID', 'SCT', 'upper_ID']].iterrows():
    adata = sc.read(f"/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/03.RNA/02.Integrate/50.Rat_All/Immu/01.commot/11.commot_maincelltype_LR_0702/{upper_id}_LR_commot_0702.h5ad")
    adata_df = pd.read_csv(f"/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/03.RNA/02.Integrate/50.Rat_All/Immu/05.st_featureplot/detailed_ST/{upper_id}_detailed_celltype.csv.gz")
  
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
  
    for fb_gene, mono_gene in zip(FB_l, Mono_r):
        ct_LR = "-".join([fb_gene, mono_gene])
        if f"s-{ct_LR}" in adata.obsm["commot-CellChat-sum-sender"].columns and f"r-{ct_LR}" in adata.obsm["commot-CellChat-sum-receiver"].columns:
            FB_s_df = adata_FB.obsm["commot-CellChat-sum-sender"][f"s-{ct_LR}"]
            FB_s_id = FB_s_df > 0
            FB_s_id_FT = FB_s_id.to_numpy().flatten()
            FB_s_exp = FB_s_df[FB_s_id].to_numpy().flatten()
            FB_s_umap = adata_FB.obsm["spatial"][FB_s_id_FT]
            
            Mono_r_df = adata_Mono.obsm["commot-CellChat-sum-receiver"][f"r-{ct_LR}"]
            Mono_r_id = Mono_r_df > 0
            Mono_r_id_FT = Mono_r_id.to_numpy().flatten()
            Mono_r_exp = Mono_r_df[Mono_r_id].to_numpy().flatten()
            Mono_r_umap = adata_Mono.obsm["spatial"][Mono_r_id_FT]
            
            exp_combined_coordinates = np.vstack((FB_s_umap, Mono_r_umap))
            exp_combined_set = set(map(tuple, exp_combined_coordinates))
            non_combined_list = [coord for coord in umap_coords if tuple(coord) not in exp_combined_set]
            non_combined_coordinates = np.array(non_combined_list)
            print(f"Number of coordinates not in combined set: {non_combined_coordinates.shape[0]}")
        
            for pt_size in [80, 140, 160, 180]:
                plt.style.use('dark_background')
                fig, ax = plt.subplots(figsize=(15, 15 / aspect_ratio))
        
                scatter3 = ax.scatter(
                    non_combined_coordinates[:, 0], non_combined_coordinates[:, 1],
                    c='#484848', label='Coordinates not in combined set', s=5
                )
        
                scatter1 = ax.scatter(
                    FB_s_umap[:, 0], FB_s_umap[:, 1],
                    c="#8fd1e1",
                    cmap=custom_cmap_fb, label=fb_gene, s=pt_size, edgecolors='none'
                )
        
                scatter2 = ax.scatter(
                    Mono_r_umap[:, 0], Mono_r_umap[:, 1],
                    c="#fedc5e",
                    cmap=custom_cmap_mono, label=mono_gene, s=pt_size, edgecolors='none'
                )
                
                plt.savefig(f"./{upper_id}_{fb_gene}_and_{mono_gene}_commot_{pt_size}.png", dpi=600)
                plt.close()
