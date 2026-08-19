#!/usr/bin/env python
# coding: utf-8

# pip install git+https://github.com/gillislab/pyMN.git --user

# In[8]:


# it works finally!
pip install git+https://github.com/gillislab/pyMN#egg=pymetaneighbor


# In[1]:


import pymn


# In[2]:


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
sc.settings.figdir = './metaneighbor_analysis/'
sc.settings.set_figure_params(dpi = 150, color_map = 'RdPu', dpi_save = 600, vector_friendly = True, format = 'pdf')


# In[3]:


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


# In[4]:


cd F:\PROJECTS\PROJECT_HUMAN_FETAL_COCHLEAE\WORKPLACE\R\sc_pipeline\


# In[5]:


data_dir = '/PROJECTS/PROJECT_HUMAN_FETAL_COCHLEAE/WORKPLACE/R/sc_pipeline/'
meta = pd.read_csv(data_dir+'metadata.csv',index_col=0)
meta


# In[6]:


adata=sc.read('human_cochlea_annotation.h5ad')
adata


# In[7]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(
    adata,legend_fontsize=8,
    color="cell_type_final",
    palette={"Erythrocytes": "#9a9ae8", #Circulating cells
             "Macrophages": "#5757dd",
             "CCs": "#45457f",

              "Tympanic_border_cells": "#806114",# Surrounding structures
            "Reissner_membrane": "#e0a31b",
             "pre-Osteoblasts": "#b59a5d",
             "pre-Osteocytes": "#726750",
             "Osteoblasts": "#c6b691",
             "Osteocytes": "#edc264",
             "SS": "#9e6e04",

             "Endothelial": "#ef7b75",# Lateral wall
             "Pericytes": "#992b28",
             "Smooth_Muscle_Cells": "#cc9997",
             "Marginal_stria": "#f23833",
             "Intermediate_stria": "#c61612",
             "Basal_stria": "#ef7048",
             "Spindle/Root_cells": "#930300",
             "Fibrocytes": "#efb4a8",

             "Interdental_cells": "#049eaa",# Supporting cells
              "Pillar_cells": "#32c3c6",
              "Deiters_cells": "#90dbdd",
              "Inner_border-phalangeal/Hensen_cells": "#0d9982",
              "Claudius/Inner-Outer_sulcus_cells": "#63b7a9",
              "SCs": "#50efd4",

              "Cochlear_HCs": "#ffde0d",# Hair cells
               "Vestibular_HCs": "#8e7c20",
             
             "Vestibular_SCs": "#efd6bb", # Vestibular
               "Vestibular_Epithelial_cells": "#f2902f",
             "Vestibular_Roof_cells": "#c4b470",
               "Vestibular_Dark_cells": "#abb21b",

             "CoE_Lateral": "#b85ff2", # early cochlear domain
               "CoE_Medial": "#8760a0",
             "CoE_Roof_cells": "#b597c9",
               "CoE_prosensory": "#ed6eff",

             "SGNs": "#89d3a4",# Neurons
               "VGNs": "#025921",
             "GCs": "#8ce264",# Glial cells

              "Mesenchymal": "#eaea78",
             "Chondrocytes": "#bcbc86",
               "Melanocytes": "#c68d8d",
       
            
            },

    size=3,


    
)


# In[8]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[11]:


#These save characters as text in PDFs
import matplotlib
import seaborn as sns
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

#These change plot aesthetics

sns.set(style='white', font_scale=1.25)
plt.rc("axes.spines", top=False, right=False)
plt.rc('xtick', bottom=True)
plt.rc('ytick', left=True)


# In[43]:


cut_labels_4 = ['< 12', '< 16','< 20','< 24','< 28', ]
cut_bins = [0, 12,16,20,24,28 ]

adata.obs["age_bins"] = pd.cut(adata.obs["gw"], bins=cut_bins, labels=cut_labels_4).astype("str")
adata.obs["age_bins"] = pd.Categorical(adata.obs["age_bins"], categories=cut_labels_4, ordered=True)
adata.obs["age_bins"]


# In[44]:


# Convert covariate to continous
adata.obs['age_bins'] = adata.obs['age_bins'].cat.codes
adata.obs['age_bins']
adata.obs.groupby(['gw','sample_id','age_bins']).apply(len)


# In[45]:


adata.obs.columns = adata.obs.columns.astype(str)


# In[46]:


adata


# In[47]:


adata.obs['cell_type_final'] = adata.obs['cell_type_final'].astype(str)
adata.obs['study_id'] = adata.obs['age_bins'].astype(str)


# In[48]:


pymn.variableGenes(adata, study_col='study_id')


# In[49]:


pymn.MetaNeighborUS(adata,
                    study_col='study_id',
                    ct_col='cell_type_final',
                    fast_version=True)


# In[50]:


adata.uns['MetaNeighborUS']


# In[51]:


adata.uns['MetaNeighborUS_params']


# In[52]:


pd.DataFrame(adata.uns['MetaNeighborUS']).to_csv("Metaneughbor_similarity_results.csv")


# In[53]:


pymn.plotMetaNeighborUS(adata, figsize=(10, 10), cmap='coolwarm', fontsize=10,)


# In[54]:


pymn.topHits(adata, threshold=0.9)

adata.uns['MetaNeighborUS_topHits']


# In[55]:


level1_split = pymn.splitClusters(
    adata, k=2,
    save_uns=False)  #Returning the splitClusters instead of saving them
level1_split


# In[56]:


first_split = level1_split[1]


# In[57]:


to_keep = np.in1d(
    pymn.join_labels(adata.obs['study_id'].values,
                     adata.obs['cell_type_final'].values), first_split)
subdata = adata[to_keep, :]
subdata.shape


# In[58]:


pymn.variableGenes(subdata, study_col='study_id')


# In[59]:


res = pymn.MetaNeighborUS(subdata,
                          study_col='study_id',
                          ct_col='cell_type_final',
                          fast_version=True)


# In[60]:


res


# In[61]:


pymn.plotMetaNeighborUS(subdata,cmap='coolwarm',figsize=(10, 10),  fontsize=5)


# In[62]:


pymn.MetaNeighborUS(adata,
                    study_col='study_id',
                    ct_col='cell_type_final',
                    fast_version=True,
                    symmetric_output=False,
                    one_vs_best=True)


# In[64]:


adata


# In[65]:


pd.DataFrame(adata.uns['MetaNeighborUS_1v1']).to_csv("Metaneughbor_similarity_results_one_vs_best.csv")


# In[68]:


pymn.plotMetaNeighborUS(adata,
                        cmap='coolwarm',
                        figsize=(18, 18),
                        mn_key='MetaNeighborUS_1v1',
                        xticklabels=True,
                        yticklabels=True,
                        fontsize=7)


# In[33]:


pymn.extractMetaClusters(adata, threshold=.7)
pymn.score_meta_clusters(adata)
mcsummary = adata.uns['MetaNeighborUS_metacluster_scores']
mcsummary[mcsummary.index != 'outliers']


# In[36]:


pd.DataFrame(mcsummary[mcsummary.index != 'outliers']).to_csv("Metaneughbor_metacluster_scores_cell_type_age_bins.csv")


# In[37]:


pymn.plotUpset(adata)


# In[42]:


pymn.makeClusterGraph(adata, low_threshold=.3)
pymn.plotClusterGraph(adata, font_size=6, figsize=(10, 10))


# In[69]:


# Pathway analysis with AUCell,https://omicverse.readthedocs.io/en/latest/Tutorials-single/t_aucell/#part3-pathway-anaylsis
import omicverse as ov
import scanpy as sc
import scvelo as scv

ov.utils.ov_plot_set()


# In[70]:


adata.X.max()# to see whether it has benn log-transformed.


# In[73]:


pathway_dict=ov.utils.geneset_prepare('GO_Biological_Process_2021.txt',organism='Human')
pathway_dict


# In[74]:


##Assest one geneset
geneset_name='actin cytoskeleton reorganization (GO:0031532)'
ov.single.geneset_aucell(adata,
                            geneset_name=geneset_name,
                            geneset=pathway_dict[geneset_name])
sc.pl.embedding(adata,
                basis='umap',
          color=["{}_aucell".format(geneset_name)])


# In[75]:


##Assest all pathways
adata_aucs=ov.single.pathway_aucell_enrichment(adata,
                                                  pathways_dict=pathway_dict,
                                                  num_workers=8)


# In[76]:


adata_aucs.obs=adata[adata_aucs.obs.index].obs
adata_aucs.obsm=adata[adata_aucs.obs.index].obsm
adata_aucs.obsp=adata[adata_aucs.obs.index].obsp
adata_aucs


# In[77]:


adata_aucs.write_h5ad('hu_coch_auce.h5ad',compression='gzip')


# In[78]:


adata_aucs=sc.read('hu_coch_auce.h5ad')


# In[84]:


sc.pl.embedding(adata_aucs,
                basis='umap',
          color=['inner ear receptor cell stereocilium organization (GO:0060122)'])


# In[85]:


#adata_aucs.uns['log1p']['base']=None
sc.tl.rank_genes_groups(adata_aucs, 'cell_type_final', method='t-test',n_genes=100)
sc.pl.rank_genes_groups_dotplot(adata_aucs,groupby='cell_type_final',
                                cmap='Spectral_r',
                                standard_scale='var',n_genes=2)


# In[88]:


sc.set_figure_params(figsize=(4,6),frameon=False,dpi=150,dpi_save=600)
sc.pl.rank_genes_groups_dotplot(adata_aucs,groupby='cell_type_final',
                                cmap='Spectral_r',save='_AUCell_GO_enrichment.pdf',
                                standard_scale='var',n_genes=3)


# In[89]:


degs = sc.get.rank_genes_groups_df(adata_aucs, group='Cochlear_HCs', key='rank_genes_groups', log2fc_min=2, 
                                    pval_cutoff=0.05)['names'].squeeze()
degs


# In[90]:


adata.uns['log1p']['base']=None
sc.tl.rank_genes_groups(adata, 'cell_type_final', method='t-test',n_genes=100)


# In[91]:


res=ov.single.pathway_enrichment(adata,pathways_dict=pathway_dict,organism='Human',
                                     group_by='cell_type_final',plot=True)


# In[98]:


sc.set_figure_params(figsize=(30,30),frameon=False,dpi=150,dpi_save=600)
ax=ov.single.pathway_enrichment_plot(res,plot_title='Enrichment',cmap='Reds',
                                         xticklabels=True,cbar=False,square=True,vmax=10,
                                         yticklabels=True,cbar_kws={'label': '-log10(qvalue)','shrink': 0.5,})


# In[101]:


res


# In[103]:


pd.DataFrame(res).to_csv("Enrichment of geneset in human cochlea.csv")


# In[ ]:





# In[ ]:




