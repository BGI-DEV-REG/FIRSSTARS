import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import time
import os
import gc
import cv2
import argparse
import seaborn as sns

os.environ["OPENCV_SHRT_MAX"] = str(pow(2, 40))

sys.path.append('/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/Axolotl_Brain_Spatial/42.Registration')
from Tools.Spatial import *


parser = argparse.ArgumentParser()
parser.add_argument('-I', '--input') 
parser.add_argument('-M', '--meta') 
parser.add_argument('-P', '--save_path')
args = parser.parse_args()
print(args)

inh5ad = args.input
outpath = args.save_path
if not os.path.exists(outpath):
    os.mkdir(outpath)

incsv = args.meta

meta = pd.read_csv(incsv, index_col=0)
meta.loc[meta['pred_celltype'] == 'Macr_Mmp14', 'pred_celltype'] = 'Macr_Cd93'
meta.loc[meta['pred_celltype'] == 'Macr_Ms4a7', 'pred_celltype'] = 'Macr_Il10'

data = sc.read(inh5ad)

data.obs = data.obs.join(meta)

color_df = pd.read_csv("/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/SkinFullyInjury/03.RNA/02.Integrate/50.Rat_All/Immu/02.mapping/color_select_7.csv")
ct = set(data.obs['pred_celltype'])
color_df_sample = color_df[color_df["celltype"].isin(ct)]

color_dict = {row['celltype']: row['color'] for _, row in color_df_sample.iterrows()}

data.obs['pred_celltype'] = data.obs['pred_celltype'].astype('category')

na_color = "#484848"

data.obs['pred_celltype'] = data.obs['pred_celltype'].cat.add_categories(['NA'])
data.obs['pred_celltype'] = data.obs['pred_celltype'].fillna('NA')
color_dict['NA'] = na_color

colors = [color_dict[ct] for ct in data.obs['pred_celltype'].cat.categories]

featureplot_slices_discrete(
    obj=data,
    feature='pred_celltype',
    fname=os.path.join(outpath, 'pred_celltype_sep0.Spatial.pdf'),
    show=False,
    scale=True,
    legend_size=6,
    slices=None,
    sep=0,
    colors=colors,
    order=data.obs['pred_celltype'].cat.categories.tolist(),
    angle_dict=data.uns['angle_dict'],
    nrow=1,
    ncol=1,
    compress_factor=False,
    raw=False
)
