# plot tacco FB
import argparse
import scanpy as sc
import matplotlib.pyplot as plt
from skimage import (color, feature, filters, measure, segmentation, io)
import numpy as np
import pandas as pd
import sys
import time
import os
import gc

os.environ["OPENCV_SHRT_MAX"] = str(pow(2, 40))
import cv2

sys.path.append('/42.Registration')
from Tools.Spatial import *

parser = argparse.ArgumentParser()
parser.add_argument('-I', '--input')
parser.add_argument('-M', '--meta')
parser.add_argument('-P', '--save_path')
parser.add_argument('-C', '--cell_type')
parser.add_argument('-D', '--idnames')
parser.add_argument('-A', '--angle')
parser.add_argument('--select', nargs='*')

args = parser.parse_args()
print(args)


inh5ad = args.input
outpath = args.save_path
if not os.path.exists(outpath):
    os.mkdir(outpath)


incsv = args.meta
meta = pd.read_csv(incsv, index_col=0)
data = sc.read(inh5ad)
angle = int(args.angle)
idnames = args.idnames

if args.select:
    meta = meta.loc[meta["pred_celltype"].isin(args.select)]
data.obs = data.obs.join(meta)

color_df = pd.read_csv("Fig2_color_seq_20240521.csv")
color_df.columns = ['ct', 'colors']
color_df_sample = color_df[color_df['ct'].isin(set(data.obs["pred_celltype"]))]

npy=data.uns[idnames]['seg_cell']
npy = segmentation.expand_labels(npy, distance=15)
data.uns[idnames]['seg_cell']=npy
data.uns['angle_dict']={idnames: angle}
featureplot_slices_discrete(obj=data,
                            feature='pred_celltype',
                            fname=os.path.join(outpath, 'pred_celltype.Spatial.pdf'),
                            blank_color=(48,48,48),
                            line_color=(0, 0, 0),
                            show=False,
                            scale=True,
                            colors = color_df_sample['colors'].to_list(),
                            order = color_df_sample['ct'].to_list(),
                            sep = 2,
                            legend_size=6,
                            slices=None,
                            angle_dict=data.uns['angle_dict'],
                            nrow=1,
                            ncol=1,
                            compress_factor=False,
                            raw=False)

