import sys
sys.path.append('/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/Axolotl_Brain_Spatial/42.Registration')
from Tools.Segmentation import Segmentation
import matplotlib.pyplot as plt
from Tools.Spatial import *
import numpy as np
import pandas as pd
import scanpy as sc
import os
import warnings
warnings.filterwarnings('ignore')

inputfile = sys.argv[1]
outputfile = sys.argv[2]
idname = sys.argv[3]
bs = sys.argv[4]
ot = sys.argv[5]
md = sys.argv[6]
ed = sys.argv[7]

bs = int(bs)
ot = float(ot)
md = int(md)
ed = int(ed)

if not os.path.exists(outputfile):
    os.mkdir(outputfile)

sobj = Segmentation()
sobj.load(
    img_path = f'{inputfile}/{idname}/{idname}_matched_ssDNA.png', 
    mRNA_path = f'{inputfile}/{idname}/{idname}.gem.gz',
)

plt.figure(figsize=(16,16))
plt.imshow(sobj.raw_img, 'gray')

path = f'{outputfile}/{idname}_matched_ssDNA.png'
cv2.imwrite(path, sobj.raw_img , [cv2.IMWRITE_PNG_COMPRESSION, 0])

# pre processing   sobj.raw_img --> sobj.img
sobj.pre_process(threshold='auto')

path = f'{outputfile}/{idname}_process_matched_ssDNA.png'
cv2.imwrite(path, sobj.img)

# watershed
sobj.watershed(
    block_size=bs,
    offset=ot,
    min_distance=md,
    expand_distance=ed
)

path = f'{outputfile}/{idname}_label_ssDNA.png'
cv2.imwrite(path, sobj.label)

# save
save_path = os.path.join(outputfile)
sobj.save_scGEM(save_path = save_path, name = idname)

