#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys
import scipy
import matplotlib
from matplotlib import rcParams
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all" # to show output from all the lines in a cells
pd.set_option('display.max_column',None) # display all the columns in pandas
pd.options.display.max_rows = 100
rcParams['pdf.fonttype'] = 42
sc.settings.figdir = './figures/processing/'
sc.settings.set_figure_params(dpi = 150, color_map = 'RdPu', dpi_save = 600, vector_friendly = True, format = 'pdf')


# In[2]:


def Barplot(which_var, adata, var='clusters', height=3, color = False):
    plotdata = pd.crosstab(adata.obs[var], adata.obs[which_var], normalize='index') * 100
    if 'category' in plotdata.index.dtype.name:
        plotdata.index.reorder_categories(adata.obs[var].cat.categories[::-1])

    if not color:
        ax1 = plotdata.plot.barh(stacked = True, edgecolor = 'none', zorder = 3, figsize = (6,height), fontsize = 14, grid = False)
    else:
        ax1 = plotdata.plot.barh(stacked = True, edgecolor = 'none', zorder = 3, figsize = (6,height), fontsize = 14, grid = False, color = color)
    ax1.set_title(which_var+' %')
    ax1.set_ylabel(var)
    horiz_offset = 1
    vert_offset = 1.
    ax1 = ax1.legend(bbox_to_anchor = (horiz_offset, vert_offset))
#     ax1.figure.savefig(str(sc.settings.figdir)+'/barplot_'+var+'_proportions_'+which_var+'.pdf', bbox_inches='tight',
#                        dpi=300, orientation='landscape', format= 'pdf', optimize=True)


# In[3]:


cd F:\PROJECTS\PROJECT_HUMAN_FETAL_COCHLEAE\WORKPLACE\R\sc_pipeline\


# In[4]:


data_dir = '/PROJECTS/PROJECT_HUMAN_FETAL_COCHLEAE/WORKPLACE/R/sc_pipeline/'
meta = pd.read_csv(data_dir+'metadata.csv',index_col=0)
meta


# In[5]:


adata = sc.read(data_dir + 'all_rawcounts.h5ad')


# In[6]:


adata


# In[7]:


#Identify HVGs
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    subset=False,
    flavor="seurat_v3",
    batch_key="gw"
)


# In[8]:


# subset object for latter
bdata = adata[:, adata.var['highly_variable']]
bdata.layers["counts"] = bdata.X.copy() # preserve counts


# In[9]:


#Normalized
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)


# In[10]:


#Load latent space scVI
X_scVI = pd.read_csv(data_dir+'all_XscVI_latent_space.csv', index_col=0)
adata.obsm["X_scVI"] = X_scVI.to_numpy()


# In[11]:


adata


# In[12]:


# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep="X_scVI",n_neighbors=100)
sc.tl.umap(adata,min_dist=0.4)#, min_dist=0.4


# In[16]:


sc.pl.umap(
    adata,
    color=['sample',  "gender",'gw'], 
    frameon=True, wspace = 0.55
)


# In[1]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata,color=['n_counts', 'percent_mito', 'scrublet_score','TMC1','OTOG','OTOGL','GATA3','TUBB3']
          ,save="_cochlea_human_raw.pdf")


# In[18]:


small_marker_dict={
    'EPITHELIUM':["EPCAM"],
     'PAN-SC':["OTOGL","OTOG","USH1C",],
    'PROSENSORY DOMAIN':['SOX2','SOX9','SOX10',"RORB","ISL1","LGR5",],
    'ROOF CELLS, LATERAL WALL':[  'OTX2','OC90',],
    'MESENCHYMAL':["PRRX1",],
    'FLOOR LATERAL':['GATA3','FGFR3','PROX1'],
    'FLOOR MEDIAL':['TECTA','FGF20'],
'HAIR CELLS':["MYO7A","STRC","OTOF","SLC26A5"],
'NEURONAL':["PRPH","TUBB3","ESRRG"],
'GLIAL':["MPZ",],
'ENDOTHELIAL' :[ "PECAM1"],
'MACROPHAGES' :[ "PTPRC"],
'MELANOCYTES' :[ "MLANA"],
'CHONDROCYTES' :[ "ACAN"],
    'CYCLING' :[ "MKI67","TOP2A","HMGB2",],
        
}

sc.pl.umap(adata,color=['EPCAM', 'OTOGL',"USH1C", 'SOX2','RORB','ISL1','LGR5','OC90','PRRX1','GATA3','FGFR3','TECTA','MYO7A',
                       'STRC','TUBB3','ESRRG','MPZ','PECAM1','PTPRC','MLANA','ACAN','MKI67']
          ,ncols=4,
           save="_cochlea_human_raw2.pdf")


# In[19]:


sc.tl.leiden(adata, resolution=1, key_added = "leiden_res1.0")


# In[20]:


sc.pl.umap(adata, color = ["leiden_res1.0", "gw"], legend_loc = "on data", legend_fontsize =5)


# In[21]:


sc.pl.umap(adata, color = ['sample'], frameon = False,  legend_fontsize = 10, )#legend_loc = "on data",


# In[22]:


sc.pl.umap(adata,color=['scrublet_score','MYO7A'],size=3)


# In[25]:


sc.pl.umap(adata,color=['scrublet_score','leiden_res1.0'],size=3,legend_loc = "on data",legend_fontsize =5)


# In[26]:


adata.obs['is_doublet_scrub'] = adata.obs['scrublet_score']>0.3


# In[27]:


np.mean(adata.obs['scrublet_score']>0.3)


# In[28]:


adata.obs['is_doublet_scrub']


# In[29]:


sc.pl.umap(adata,color=['leiden_res1.0','scrublet_score'],size=4,legend_loc = "on data",legend_fontsize =8)


# In[61]:


#subdivide cluster 8 
adata_db_rm=adata.copy()
sc.tl.leiden(adata_db_rm, resolution=2, restrict_to=('leiden_res1.0', ['8']))
sc.pl.umap(adata_db_rm, color=['leiden_R','scrublet_score'],size=3,legend_fontsize =4,)#legend_loc='on data'


# In[62]:


adata_db_rm.obs['leiden_res1.0']=adata_db_rm.obs['leiden_R']


# In[63]:


exclude_clusters = ['8,4','8,10','8,11']
adata_db_rm= adata_db_rm[~adata_db_rm.obs['leiden_res1.0'].isin(exclude_clusters), :]


# In[64]:


sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=4,legend_fontsize =4)#legend_loc='on data'


# In[65]:


sc.tl.leiden(adata_db_rm, resolution=1, key_added = "leiden_res1.0")


# In[66]:


sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=4,legend_fontsize =8,legend_loc='on data')#legend_loc='on data'


# In[72]:


#subdivide cluster 9 
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.leiden(adata_db_rm, resolution=2, restrict_to=('leiden_res1.0', ['9']))
sc.pl.umap(adata_db_rm, color=['leiden_R','scrublet_score'],size=3,legend_fontsize =4,)#legend_loc='on data'


# In[73]:


sc.set_figure_params(figsize=(8,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_db_rm, ['scrublet_score'], 
             
             groupby='leiden_R',rotation=90)


# In[74]:


adata_db_rm.obs['leiden_res1.0']=adata_db_rm.obs['leiden_R']


# In[75]:


exclude_clusters = ['9,9','9,12']
adata_db_rm= adata_db_rm[~adata_db_rm.obs['leiden_res1.0'].isin(exclude_clusters), :]


# In[76]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=4,legend_fontsize =4)#legend_loc='on data'


# In[77]:


sc.tl.leiden(adata_db_rm, resolution=1, key_added = "leiden_res1.0")


# In[81]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=3,legend_fontsize =10,legend_loc='on data')#legend_loc='on data'


# In[82]:


#subdivide cluster 13 
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.leiden(adata_db_rm, resolution=2, restrict_to=('leiden_res1.0', ['13']))
sc.pl.umap(adata_db_rm, color=['leiden_R','scrublet_score'],size=3,legend_fontsize =4,)#legend_loc='on data'


# In[83]:


sc.set_figure_params(figsize=(8,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_db_rm, ['scrublet_score'], 
             
             groupby='leiden_R',rotation=90)


# In[84]:


adata_db_rm.obs['leiden_res1.0']=adata_db_rm.obs['leiden_R']


# In[85]:


exclude_clusters = ['13,0',]
adata_db_rm= adata_db_rm[~adata_db_rm.obs['leiden_res1.0'].isin(exclude_clusters), :]


# In[86]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=4,legend_fontsize =4)#legend_loc='on data'


# In[87]:


sc.tl.leiden(adata_db_rm, resolution=1, key_added = "leiden_res1.0")


# In[88]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=3,legend_fontsize =10,legend_loc='on data')#legend_loc='on data'


# In[89]:


#subdivide cluster 15 
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.leiden(adata_db_rm, resolution=2, restrict_to=('leiden_res1.0', ['15']))
sc.pl.umap(adata_db_rm, color=['leiden_R','scrublet_score'],size=3,legend_fontsize =4,)#legend_loc='on data'


# In[90]:


sc.set_figure_params(figsize=(8,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_db_rm, ['scrublet_score'], 
             
             groupby='leiden_R',rotation=90)


# In[91]:


adata_db_rm.obs['leiden_res1.0']=adata_db_rm.obs['leiden_R']


# In[92]:


exclude_clusters = ['15,0',]
adata_db_rm= adata_db_rm[~adata_db_rm.obs['leiden_res1.0'].isin(exclude_clusters), :]


# In[93]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=4,legend_fontsize =4)#legend_loc='on data'


# In[94]:


sc.tl.leiden(adata_db_rm, resolution=1, key_added = "leiden_res1.0")


# In[95]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=3,legend_fontsize =10,legend_loc='on data')#legend_loc='on data'


# In[96]:


#subdivide cluster 2 
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.leiden(adata_db_rm, resolution=3, restrict_to=('leiden_res1.0', ['2']))
sc.pl.umap(adata_db_rm, color=['leiden_R','scrublet_score'],size=3,legend_fontsize =4,)#legend_loc='on data'


# In[98]:


sc.set_figure_params(figsize=(10,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_db_rm, ['scrublet_score'], 
             
             groupby='leiden_R',rotation=90)


# In[99]:


adata_db_rm.obs['leiden_res1.0']=adata_db_rm.obs['leiden_R']
exclude_clusters = ['2,2','2,5','2,25','2,30','2,31','2,32','2,33']
adata_db_rm= adata_db_rm[~adata_db_rm.obs['leiden_res1.0'].isin(exclude_clusters), :]
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=4,legend_fontsize =4)#legend_loc='on data'


# In[100]:


sc.tl.leiden(adata_db_rm, resolution=1, key_added = "leiden_res1.0")
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=3,legend_fontsize =10,legend_loc='on data')#legend_loc='on data'


# In[101]:


#subdivide cluster 2 
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.leiden(adata_db_rm, resolution=3, restrict_to=('leiden_res1.0', ['2']))
sc.pl.umap(adata_db_rm, color=['leiden_R','scrublet_score'],size=3,legend_fontsize =4,)#legend_loc='on data'


# In[102]:


sc.set_figure_params(figsize=(10,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_db_rm, ['scrublet_score'],    
             groupby='leiden_R',rotation=90)


# In[103]:


adata_db_rm.obs['leiden_res1.0']=adata_db_rm.obs['leiden_R']
exclude_clusters = ['2,5','2,11','2,9','2,19','2,21','2,23','2,28','2,32','2,33','2,34']
adata_db_rm= adata_db_rm[~adata_db_rm.obs['leiden_res1.0'].isin(exclude_clusters), :]
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=4,legend_fontsize =4)#legend_loc='on data'


# In[104]:


sc.tl.leiden(adata_db_rm, resolution=1, key_added = "leiden_res1.0")
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=3,legend_fontsize =10,legend_loc='on data')#legend_loc='on data'


# In[108]:


#subdivide cluster 15 
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.leiden(adata_db_rm, resolution=2.4, restrict_to=('leiden_res1.0', ['15']))
sc.pl.umap(adata_db_rm, color=['leiden_R','scrublet_score'],size=3,legend_fontsize =4,)#legend_loc='on data'


# In[109]:


sc.set_figure_params(figsize=(10,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_db_rm, ['scrublet_score'],    
             groupby='leiden_R',rotation=90)


# In[110]:


adata_db_rm.obs['leiden_res1.0']=adata_db_rm.obs['leiden_R']
exclude_clusters = ['15,1','15,3','15,8','15,11','15,15','15,17','15,19','15,20','15,23','15,24']
adata_db_rm= adata_db_rm[~adata_db_rm.obs['leiden_res1.0'].isin(exclude_clusters), :]
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=4,legend_fontsize =4)#legend_loc='on data'


# In[111]:


sc.tl.leiden(adata_db_rm, resolution=1, key_added = "leiden_res1.0")
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=3,legend_fontsize =10,legend_loc='on data')#legend_loc='on data'


# In[112]:


#subdivide cluster 8 
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.leiden(adata_db_rm, resolution=2.4, restrict_to=('leiden_res1.0', ['8']))
sc.pl.umap(adata_db_rm, color=['leiden_R','scrublet_score'],size=3,legend_fontsize =4,)#legend_loc='on data'


# In[113]:


sc.set_figure_params(figsize=(10,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_db_rm, ['scrublet_score'],    
             groupby='leiden_R',rotation=90)


# In[114]:


adata_db_rm.obs['leiden_res1.0']=adata_db_rm.obs['leiden_R']
exclude_clusters = ['8,3','8,7','8,12','8,14','8,15','8,23','8,26','8,27',]
adata_db_rm= adata_db_rm[~adata_db_rm.obs['leiden_res1.0'].isin(exclude_clusters), :]
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=4,legend_fontsize =4)#legend_loc='on data'


# In[115]:


sc.tl.leiden(adata_db_rm, resolution=1, key_added = "leiden_res1.0")
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=3,legend_fontsize =10,legend_loc='on data')#legend_loc='on data'


# In[117]:


#subdivide cluster 15 
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.leiden(adata_db_rm, resolution=4, restrict_to=('leiden_res1.0', ['15']))
sc.pl.umap(adata_db_rm, color=['leiden_R','scrublet_score'],size=3,legend_fontsize =4,)#legend_loc='on data'
sc.set_figure_params(figsize=(10,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_db_rm, ['scrublet_score'],    
             groupby='leiden_R',rotation=90)


# In[118]:


sc.tl.leiden(adata_db_rm, resolution=1, key_added = "leiden_res1.0")
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=3,legend_fontsize =10,legend_loc='on data')#legend_loc='on data'


# In[119]:


#subdivide cluster 2 
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.leiden(adata_db_rm, resolution=3, restrict_to=('leiden_res1.0', ['2']))
sc.pl.umap(adata_db_rm, color=['leiden_R','scrublet_score'],size=3,legend_fontsize =4,)#legend_loc='on data'
sc.set_figure_params(figsize=(10,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_db_rm, ['scrublet_score'],    
             groupby='leiden_R',rotation=90)


# In[120]:


adata_db_rm.obs['leiden_res1.0']=adata_db_rm.obs['leiden_R']


# In[121]:


sc.tl.leiden(adata_db_rm, resolution=1, key_added = "leiden_res1.0")
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=3,legend_fontsize =10,legend_loc='on data')#legend_loc='on data'


# In[122]:


#subdivide cluster 3 
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.leiden(adata_db_rm, resolution=3, restrict_to=('leiden_res1.0', ['3']))
sc.pl.umap(adata_db_rm, color=['leiden_R','scrublet_score'],size=3,legend_fontsize =4,)#legend_loc='on data'
sc.set_figure_params(figsize=(10,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_db_rm, ['scrublet_score'],    
             groupby='leiden_R',rotation=90)


# In[123]:


adata_db_rm.obs['leiden_res1.0']=adata_db_rm.obs['leiden_R']
exclude_clusters = ['3,21','3,33','3,34','3,35','3,36']
adata_db_rm= adata_db_rm[~adata_db_rm.obs['leiden_res1.0'].isin(exclude_clusters), :]
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=4,legend_fontsize =4)#legend_loc='on data'


# In[124]:


sc.tl.leiden(adata_db_rm, resolution=1, key_added = "leiden_res1.0")
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=3,legend_fontsize =10,legend_loc='on data')#legend_loc='on data'


# In[125]:


#subdivide cluster 12 
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.leiden(adata_db_rm, resolution=3, restrict_to=('leiden_res1.0', ['12']))
sc.pl.umap(adata_db_rm, color=['leiden_R','scrublet_score'],size=3,legend_fontsize =4,)#legend_loc='on data'
sc.set_figure_params(figsize=(10,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_db_rm, ['scrublet_score'],    
             groupby='leiden_R',rotation=90)


# In[126]:


adata_db_rm.obs['leiden_res1.0']=adata_db_rm.obs['leiden_R']
exclude_clusters = ['12,18','12,27',]
adata_db_rm= adata_db_rm[~adata_db_rm.obs['leiden_res1.0'].isin(exclude_clusters), :]
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=4,legend_fontsize =4)#legend_loc='on data'


# In[127]:


sc.tl.leiden(adata_db_rm, resolution=1, key_added = "leiden_res1.0")
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=3,legend_fontsize =10,legend_loc='on data')#legend_loc='on data'


# In[129]:


sc.tl.leiden(adata_db_rm, resolution=2, key_added = "leiden_res2.0")
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res2.0','scrublet_score'],size=3,legend_fontsize =10,legend_loc='on data')#legend_loc='on data'


# In[130]:


sc.pl.umap(adata_db_rm, color=['leiden_res2.0','scrublet_score'],size=3,legend_fontsize =10,)#legend_loc='on data'


# In[131]:


exclude_clusters = ['4',]
adata_db_rm= adata_db_rm[~adata_db_rm.obs['leiden_res2.0'].isin(exclude_clusters), :]
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_db_rm, color=['leiden_res1.0','scrublet_score'],size=4,legend_fontsize =4)#legend_loc='on data'


# In[132]:


import scvi


# In[133]:


def run_scvi(adata_db_rm, batch_hv = "sample", batch_scvi = "gw", \
             cat_cov_scvi = ["sample",  "sequencing", "gender","dataset"], cont_cov_scvi = ["percent_mito"], \
             include_genes = [], exclude_cc_genes=True, vae_name = "", **kwargs):
    #adata_scvi = sc.AnnData(X = adata.layers['counts'].copy(), obs = adata.obs.copy(), var = adata.var.copy())
    adata_scvi = adata_db_rm.copy()
    adata_scvi.layers["counts"] = adata_db_rm.X.copy()
    sc.pp.normalize_total(adata_scvi, target_sum=1e4)
    sc.pp.log1p(adata_scvi)
     # keep full dimension safe
    sc.pp.highly_variable_genes(adata_scvi, flavor="seurat_v3", n_top_genes=10000, layer="counts",\
                                batch_key=batch_hv)
    selected_genes = list(set(adata_scvi.var.loc[adata_scvi.var['highly_variable']].index.tolist()+ include_genes))
    print(len(selected_genes))

    adata_scvi = adata_scvi[:, selected_genes].copy()
    scvi.model.SCVI.setup_anndata(adata_scvi, layer="counts", batch_key=batch_scvi, \
                             categorical_covariate_keys=cat_cov_scvi, \
                             continuous_covariate_keys=cont_cov_scvi)
    scvi_kwargs = {k: v for k,v in kwargs.items() if k in scvi.model.SCVI.__init__.__code__.co_varnames}
    vae = scvi.model.SCVI(adata_scvi, **scvi_kwargs)
    train_kwargs = {k: v for k,v in kwargs.items() if k in vae.train.__code__.co_varnames}
    vae.train(**train_kwargs)
    adata_scvi.obsm["X_scVI"] = vae.get_latent_representation()
    sc.pp.neighbors(adata_scvi, use_rep="X_scVI")
    sc.tl.leiden(adata_scvi)
    sc.tl.umap(adata_scvi)
    sc.pl.umap(
        adata_scvi,
        color=["sample", "gw"],
        frameon=False,
        ncols=2,
    )
    sc.pl.umap(
        adata_scvi,
        color=['n_counts', 'percent_mito', 
               'scrublet_score'])
    
    return(adata_scvi)


# In[134]:


adata_db_rm.layers["counts"] = adata_db_rm.X.copy() # preserve counts


# In[135]:


adata_scvi = run_scvi(adata_db_rm, batch_hv = "sample", batch_scvi = "gw", \
             cat_cov_scvi = ["sample",  "sequencing", "gender","dataset"], cont_cov_scvi = ["percent_mito"], \
             #include_genes = all_mrkrs,  
                      #vae_name = "cochlea_snuclei", 
   n_hidden=512, n_layers=2, n_latent=50, gene_likelihood='nb', dispersion='gene-batch',use_observed_lib_size=False,
                           train_size=0.99, max_epochs=1000, batch_size=1024, early_stopping = True)


# In[151]:


adata_scvi


# In[153]:


vae = scvi.model.SCVI(adata_scvi, n_layers=3, n_latent=30)
vae


# In[154]:


vae.view_anndata_setup()


# In[156]:


max_epochs_scvi = np.min([round((100000 / adata_scvi.n_obs) * 400), 400])
max_epochs_scvi


# In[157]:


vae.train()


# In[158]:


# Normally we would need to run scVI first but we have already done that here
# model_scvi = scvi.model.SCVI(adata_scvi) etc.
model_scanvi = scvi.model.SCANVI.from_scvi_model(
    vae, labels_key='cell_type', unlabeled_category="unlabelled"
)#unlabelled
print(model_scanvi)
model_scanvi.view_anndata_setup()


# In[159]:


max_epochs_scanvi = int(np.min([10, np.max([2, round(max_epochs_scvi / 3.0)])]))
model_scanvi.train(max_epochs=max_epochs_scanvi)


# In[183]:


adata_scanvi = adata_scvi.copy()
adata_scanvi.obsm["X_scANVI"] = model_scanvi.get_latent_representation()
sc.pp.neighbors(adata_scanvi, use_rep="X_scANVI",n_neighbors=30)#n_neighbors=50
sc.tl.leiden(adata_scanvi,resolution=2,key_added = "leiden_scanvi_2.0")
sc.tl.umap(adata_scanvi,min_dist=0.8)#0.4
sc.pl.umap(adata_scanvi, color=['cell_type','gw','leiden_scanvi_2.0'], wspace=0.5,ncols=1)


# In[184]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.tl.umap(adata_scanvi,min_dist=0.8)#0.4
sc.pl.umap(adata_scanvi, color=['cell_type'], wspace=0.5,ncols=1,legend_loc='on data', legend_fontsize=3)


# In[207]:


# remove clusters that are not combined into one cell type.
sc.settings.figdir = './figures/annotation/'


# In[186]:


adata_scvi.write('human_cochlea_scVI.h5ad')


# In[187]:


adata_scanvi.write('human_cochlea_scanVI.h5ad')


# In[202]:


adata_scanvi2=adata_scanvi.copy()
adata_scanvi2=adata_scanvi2.raw.to_adata()
adata_scanvi2


# In[209]:


adata_scanvi2.layers["counts"] = adata_db_rm.X.copy() # preserve counts
adata_scanvi2


# In[217]:


adata_scanvi2.write('human_cochlea_scanVI.h5ad')


# In[210]:


adata_db_rm


# In[211]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_scanvi2, color=['cell_type'], wspace=0.5,ncols=1,legend_loc='on data', legend_fontsize=3)


# In[189]:


#Save latent space
pd.DataFrame(adata_scvi.obsm["X_scVI"]).to_csv(data_dir + 'all_XscVI_latent_space_db_rm.csv')


# In[190]:


#Save latent space
pd.DataFrame(adata_scanvi.obsm["X_scANVI"]).to_csv(data_dir + 'all_XscANVI_latent_space_db_rm.csv')


# In[193]:


#Load latent space scVI
X_scVI = pd.read_csv(data_dir+'all_XscVI_latent_space_db_rm.csv', index_col=0)
adata_db_rm.obsm["X_scVI"] = X_scVI.to_numpy()


# In[194]:


#Load latent space scVI
X_scANVI = pd.read_csv(data_dir+'all_XscANVI_latent_space_db_rm.csv', index_col=0)
adata_db_rm.obsm["X_scANVI"] = X_scVI.to_numpy()


# In[212]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_scanvi2,color=['n_counts', 'percent_mito', 'scrublet_score','TMC1','OTOG','OTOGL','GATA3','TUBB3']
          ,save="_cochlea_human_scanvi.pdf")


# In[213]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_scanvi2,color=['EPCAM', 'OTOGL',"USH1C", 'SOX2','RORB','ISL1','LGR5','OC90','PRRX1','GATA3','FGFR3','TECTA','MYO7A',
                       'STRC','TUBB3','ESRRG','MPZ','PECAM1','PTPRC','MLANA','ACAN','MKI67']
          ,ncols=4,
           save="_cochlea_human_scanvi2.pdf")


# In[214]:


adata_db_rm.write('all_db_rm.h5ad')


# In[215]:


sc.set_figure_params(figsize=(8,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_scanvi2, ['scrublet_score'], 
             
             groupby='cell_type',rotation=90)


# In[216]:


# remove clusters that are not combined into one cell type.
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_scanvi2, color=['cell_type'], wspace=0.5,ncols=1,
           legend_fontsize=8,
          save='_human_cochlea_scanvi.pdf') #legend_loc='on data',


# In[ ]:




