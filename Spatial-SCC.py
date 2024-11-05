import scanpy as sc
import matplotlib.pyplot as plt
import squidpy as sq
import numpy as np
import pandas as pd
import sys
import time
import os
os.environ["OPENCV_SHRT_MAX"] = str(pow(2,40))
import cv2
#import dynamo as dyn
sys.path.append('/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/Axolotl_Brain_Spatial/42.Registration')
from Tools.Spatial import *

inh5ad = sys.argv[1]
inputpath = sys.argv[2]
outpath = sys.argv[3]
idnames = sys.argv[4]
e_neigh = sys.argv[5]
s_neigh = sys.argv[6]
resolution = sys.argv[7]
Seurat_meta = sys.argv[8]

e_neigh = int(e_neigh)
s_neigh = int(s_neigh)
resolution = float(resolution)
data = sc.read(inh5ad)
data.var_names_make_unique()

sc.pp.pca(data,n_comps=30)

#### Consider Spatial distance Clustering

def spatially_constrained_cluster(adata, e_neigh=30, s_neigh=6, binsize=None,resolution = 1):
#     sc.pp.filter_genes(adata,min_cells=1)
#     sc.pp.normalize_total(adata)
#     sc.pp.log1p(adata)
#     sc.pp.pca(adata,n_comps=30)
    sc.pp.neighbors(adata, n_neighbors=e_neigh)
    sq.gr.spatial_neighbors(adata, n_neighs=s_neigh)
    conn = adata.obsp['connectivities'].copy()
    conn.data[conn.data > 0] = 1
    adj = conn + adata.obsp['spatial_connectivities']
    adj.data[adj.data  > 0] = 1
    sc.tl.leiden(adata, adjacency=adj, key_added='spatial_leiden_e' + str(e_neigh) + '_s' + str(s_neigh),resolution = resolution)

spatially_constrained_cluster(data,e_neigh=e_neigh, s_neigh=s_neigh, binsize=None,resolution = resolution)
#dyn.pl.space(data, genes='spatial_leiden_e' + str(e_neigh) + '_s' + str(s_neigh), show_legend=True, pointsize=0.2)
#plt.savefig(outpath+'/'+idnames+'.dyn.pdf')

#### Spatial plot

a = np.load(inputpath + '/' + str(idnames) + '_mask.npy')
#a= a[~(a==0).all(1)]

data.uns[str(idnames)]={}
data.uns[str(idnames)]['seg_cell']=a
data.obs['Batch'] = str(idnames)

#name.columns = 'CellID'
data.obs['cell_id'] = data.obs['CellID'].apply(lambda x:int(x.split('.')[1])).values
data.obs['cell_id']= data.obs['cell_id'].astype('category')

data.uns['angle_dict']={}
data.uns['angle_dict'][str(idnames)]=0

color_df = pd.read_csv('/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Axolotl_Brain_Spatial_Atlas/01.Segmentation/00.data/Color.seq.txt', sep='\t')
color_df = color_df.sort_values(by='order')
color_df['order'] = color_df['order'].astype('str')
ct = list(set(data.obs['spatial_leiden_e' + str(e_neigh) + '_s' + str(s_neigh)]))
color_df_sample = color_df[color_df['order'].isin(ct)]

featureplot_slices_discrete(obj = data,
    feature = 'spatial_leiden_e' + str(e_neigh) + '_s' + str(s_neigh),
    fname = os.path.join(outpath,str(idnames)+'.Segmentation.Spatial.pdf'),
    show = False,
    scale = True,
    legend_size = 6,
    order = color_df_sample['order'].tolist(),
    colors = color_df_sample['Color'].tolist(),
    slices = None,
    angle_dict = data.uns['angle_dict'],
    nrow = 1,
    ncol = 1,
    compress_factor=False,
    raw=False)

Seu = pd.read_csv(Seurat_meta, index_col = 0)
data.obs['seurat_clusters'] = Seu['seurat_clusters']

data.obs['seurat_clusters'] = data.obs['seurat_clusters'].astype('str')
ct = list(set(data.obs['seurat_clusters']))
color_df_sample = color_df[color_df['order'].isin(ct)]

featureplot_slices_discrete(obj = data,
    feature = 'seurat_clusters',
    fname = os.path.join(outpath,str(idnames)+'.seurat_clusters.Spatial.pdf'),
    show = False,
    scale = True,
    legend_size = 6,
    order = color_df_sample['order'].tolist(),
    colors = color_df_sample['Color'].tolist(),
    slices = None,
    angle_dict = data.uns['angle_dict'],
    nrow = 1,
    ncol = 1,
    #sep = 0,
    compress_factor=False,
    raw=False)


data.write(os.path.join(outpath,str(idnames)+'.Spatial.distance.cluster.h5ad'))
### find Markers
Marker = outpath+'/Markers'
if not os.path.isdir(Marker):
    os.mkdir(Marker)

#os.mkdir(outpath+'/Markers')
marker_path =os.path.join(outpath+'/Markers')

sc.tl.rank_genes_groups(data, 'spatial_leiden_e' + str(e_neigh) + '_s' + str(s_neigh), method='wilcoxon')
pvals = pd.DataFrame(data.uns['rank_genes_groups']['pvals'])
pvals.to_csv(os.path.join(marker_path,str(idnames)+'.Spatial.distance.cluster.pvals.csv'))

pvals_adj = pd.DataFrame(data.uns['rank_genes_groups']['pvals_adj'])
pvals_adj.to_csv(os.path.join(marker_path,str(idnames)+'.Spatial.distance.cluster.pvals_adj.csv'))

logfoldchanges = pd.DataFrame(data.uns['rank_genes_groups']['logfoldchanges'])
logfoldchanges.to_csv(os.path.join(marker_path,str(idnames)+'.Spatial.distance.cluster.logfoldchanges.csv'))

names = pd.DataFrame(data.uns['rank_genes_groups']['names'])
names.to_csv(os.path.join(marker_path,str(idnames)+'.Spatial.distance.cluster.names.csv'))
names.columns = 'cluster' + names.columns

axolotl_gene = pd.read_csv('/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/Axolotl_Brain_Spatial_V2/1.Figure1_Data_0711/Summary_Gene_Annotaion_0707.xls',sep='\t',header=0)

for index, row in names.iteritems():
    print(index)
    an = names[index]
    an = an.to_frame()
    ano = pd.merge(an, axolotl_gene, how='left', left_on = index, right_on='Axolotl_ID')
    ano.to_csv(os.path.join(marker_path, str(idnames) + '.Spatial.distance.cluster.'+ str(index) + '.csv'))
