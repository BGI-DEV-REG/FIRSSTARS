import argparse
import scanpy as sc
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import sys
import time
import os
import gc

os.environ["OPENCV_SHRT_MAX"] = str(pow(2, 40))
import cv2
import ast

sys.path.append('/42.Registration')
from Tools.Spatial import *

parser = argparse.ArgumentParser()
parser.add_argument('-I', '--input')
parser.add_argument('-M', '--meta')
parser.add_argument('-P', '--save_path')
parser.add_argument('-D', '--idnames')
parser.add_argument('-A', '--angle')
args = parser.parse_args()
print(args)


inh5ad = args.input
outpath = args.save_path
if not os.path.exists(outpath):
    os.mkdir(outpath)


incsv = args.meta
id = args.idnames

meta = pd.read_csv(incsv, index_col=0)
data = sc.read(inh5ad)
data.obs['pred_celltype']=meta['pred_celltype2']
print(data.obs['pred_celltype'].value_counts())

color_df = pd.read_csv('Fig1_color_seq_20240512.csv')
ct = set(data.obs.pred_celltype)
color_df_sample = color_df[color_df["celltype"].isin(ct)]
data.uns['angle_dict']=ast.literal_eval(args.angle)

featureplot_single_discrete(obj = data,
                            feature = 'pred_celltype',
                            fname=os.path.join(outpath, id+'.pred_celltype_split.Spatial.pdf'),
                            show = True,
                            scale = True,
                            legend_size = 6,
                            order = color_df_sample['celltype'].tolist(),
                            colors = color_df_sample['color'].tolist(),
                            slice = str(id),
                            nrow = 7,
                            ncol = 4,
                            compress=True,
                            raw=False
)


