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

sys.path.append('/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/Axolotl_Brain_Spatial/42.Registration')
from Tools.Spatial import *

parser = argparse.ArgumentParser()
parser.add_argument('-I', '--input')
parser.add_argument('-F', '--feature')
parser.add_argument('-A', '--angle')
parser.add_argument('-M', '--meta')
parser.add_argument('-P', '--save_path')
parser.add_argument('--prefix', type=str, default="plot_")

args = parser.parse_args()
print(args)

angle = int(args.angle)
inh5ad = args.input
outpath = args.save_path
if not os.path.exists(outpath):
    os.mkdir(outpath)

data = sc.read_h5ad(inh5ad)
meta = pd.read_csv(args.meta, index_col=0)
data.obs[args.feature] = meta[args.feature]
data.uns['angle_dict'] = {args.prefix:angle}
featureplot_slices_continuous(obj=data,
                              feature=args.feature,
			      #cmap='plasma',
                              fname=os.path.join(outpath, args.prefix + args.feature + '.Spatial.pdf'),
                              show=False,
                              scale=True,
                              slices=None,
                              angle_dict=data.uns['angle_dict'],
                              nrow=1,
                              ncol=1,
                              compress_factor=False,
                              raw=False)
