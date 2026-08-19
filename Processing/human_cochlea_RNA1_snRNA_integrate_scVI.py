#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

import scanpy as sc
import scvi

from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all" # to show output from all the lines in a cells
pd.set_option('display.max_column',None) # display all the columns in pandas
pd.options.display.max_rows = 100

from datetime import date
today = str(date.today())

import matplotlib
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42
sc.settings.set_figure_params(dpi = 150, color_map = 'RdPu', dpi_save = 600, vector_friendly = True, format = 'pdf')

#%cd /...
import scrublet as scr
from doubletdetection import *

#folder = '/...


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


palette = ["#E31A1C", "#1F78B4", "#A6CEE3",  "#B2DF8A", "#33A02C", "#FB9A99",  "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#66C2A5",
               "#FC8D62", "#8DA0CB", "#B3B3B3", "#A6D854", "#FFD92F", "#E5C494", "#E78AC3"]


# # set up scVI Environment

# In[4]:


import scvi
import scanpy as sc
import pandas as pd 
import numpy as np
sc.set_figure_params(figsize=(4, 4))


# In[5]:


cd F:\PROJECTS\PROJECT_HUMAN_FETAL_COCHLEAE\WORKPLACE\R\sc_pipeline\


# In[6]:


data_dir = '/PROJECTS/PROJECT_HUMAN_FETAL_COCHLEAE/WORKPLACE/R/sc_pipeline/'
meta = pd.read_csv(data_dir+'metadata.csv',index_col=0)
meta


# In[7]:


adata = sc.read(data_dir + 'all_rawcounts.h5ad')


# In[8]:


adata


# In[9]:


adata.layers["counts"] = adata.X.copy() # preserve counts


# In[11]:


adata.obs['sample'].value_counts()


# In[16]:


def run_scvi(adata, batch_hv = "sample", batch_scvi = "gw", \
             cat_cov_scvi = ["sample",  "sequencing", "gender","dataset"], cont_cov_scvi = ["percent_mito"], \
             include_genes = [], exclude_cc_genes=True, vae_name = "", **kwargs):
    #adata_scvi = sc.AnnData(X = adata.layers['counts'].copy(), obs = adata.obs.copy(), var = adata.var.copy())
    adata_scvi = adata.copy()
    adata_scvi.layers["counts"] = adata_scvi.X.copy()
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
    
    
    return(adata_scvi)


# In[18]:


adata_scvi = run_scvi(adata, batch_hv = "sample", batch_scvi = "gw", \
             cat_cov_scvi = ["sample",  "sequencing", "gender","dataset"], cont_cov_scvi = ["percent_mito"], \
             #include_genes = all_mrkrs,  
                      #vae_name = "cochlea_snuclei", 
   n_hidden=512, n_layers=2, n_latent=50, gene_likelihood='nb', dispersion='gene-batch',use_observed_lib_size=False,
                           train_size=0.99, max_epochs=1000, batch_size=1024, early_stopping = True)


# In[22]:


adata_scvi


# In[27]:


sc.pp.neighbors(adata_scvi, use_rep="X_scVI")
sc.tl.leiden(adata_scvi)
sc.tl.umap(adata_scvi)    


# In[29]:


adata_raw_scvi = adata.copy()
adata_raw_scvi.obsm['X_scVI'] = adata_scvi.obsm['X_scVI'].copy()
adata_raw_scvi.obsm['X_umap'] = adata_scvi.obsm['X_umap'].copy()
adata_raw_scvi.obsp = adata_scvi.obsp.copy()
adata_raw_scvi.uns = adata_scvi.uns.copy()
adata_raw_scvi.obs['leiden'] = adata_scvi.obs['leiden'].copy()
    
sc.pl.umap(
        adata_raw_scvi,
        color=[ "gw",  "sample",  "sequencing", "gender","dataset"],
        frameon=False,
        ncols=2,
    )
sc.pl.umap(
        adata_raw_scvi,
        color=['n_counts',  'percent_mito',  \
               'scrublet_score', 'is_doublet'])


# In[30]:


sc.pl.umap(
        adata_raw_scvi,
        color=[ "cell_type", ],
        frameon=False,
        ncols=2,
    )


# In[31]:


sc.settings.figdir = './figures/processing/'


# In[32]:


#Save latent space
pd.DataFrame(adata_raw_scvi.obsm["X_scVI"]).to_csv(data_dir + 'all_XscVI_latent_space.csv')


# In[33]:


adata_raw_scvi.write('scvi_integrated_raw.h5ad')


# In[ ]:




