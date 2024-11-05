import scanpy as sc
from matplotlib.backends.backend_pdf import PdfPages
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

os.chdir("../fig1_sup/")
markers = pd.read_csv("Fig1_All_sc_celltype_marker.csv")
markers["marker"] = markers["marker"].str.split(',')
marker_genes = markers["marker"].sum()

adata = sc.read("change_celltype_0520.h5ad")

adata.obs.celltype = adata.obs.celltype.cat.reorder_categories(['FB_Injury','FB_Papillary','FB_Reticular','FB_Fascia','FB_DP','FB_DS','Immu_Myeloid','Immu_Neu','Immu_pDC_B','Immu_T_NK','KC_Basal','KC_Suprabasal','KC_SG','KC_HF','KC_HFSC','KC_IRS','KC_ORS','KC_TAC','Melanocyte','Muscle','Pericyte','Schwann','Endothelial','Erythrocyte'],ordered=True)

sc.pl.dotplot(adata, marker_genes, groupby='celltype',cmap='OrRd',
        save='Fig1_celltype_markers_heatmap_white.pdf',vmin=0,vmax=2,figsize=(20,10),swap_axes=False)