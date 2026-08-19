#!/usr/bin/env python
# coding: utf-8

# In[71]:


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
sc.settings.figdir = './figures/annotation/'
sc.settings.set_figure_params(dpi = 150, color_map = 'RdPu', dpi_save = 600, vector_friendly = True, format = 'pdf')


# In[72]:


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


# In[73]:


cd F:\PROJECTS\PROJECT_HUMAN_FETAL_COCHLEAE\WORKPLACE\R\sc_pipeline\


# In[74]:


data_dir = '/PROJECTS/PROJECT_HUMAN_FETAL_COCHLEAE/WORKPLACE/R/sc_pipeline/'
meta = pd.read_csv(data_dir+'metadata.csv',index_col=0)
meta


# In[22]:


adata = sc.read(data_dir + 'human_cochlea_scanVI.h5ad')


# In[23]:


adata


# In[17]:


# remove clusters that are not combined into one cell type.
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata, color=['cell_type'], wspace=0.5,ncols=1,
           legend_fontsize=8,
          save='_human_cochlea_scanvi.pdf') #legend_loc='on data',


# In[24]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(
    adata,
    color=['gw', 'scrublet_score',  'percent_mito', 'gender'], color_map = "PiYG",legend_loc="on data",
    frameon=True, ncols = 2, wspace = 0.3
)


# In[16]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata,color=['EPCAM', 'OTOGL',"USH1C", 'SOX2','RORB','ISL1','LGR5','OC90','PRRX1','GATA3','FGFR3','TECTA','MYO7A',
                       'STRC','TUBB3','ESRRG','MPZ','PECAM1','PTPRC','MLANA','ACAN','MKI67','gw', 'scrublet_score',]
          ,ncols=4,color_map = "inferno",size=3,
           save="_human_cochlea_scanvi2.pdf")


# In[77]:


# remove clusters that are not combined into one cell type.
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata, color=['leiden_scanvi_4.0'], wspace=0.5,ncols=1,
           legend_fontsize=8,legend_loc='on data',size=2,
          save='_human_cochlea_scanvi3.pdf') #legend_loc='on data',


# In[69]:


# remove clusters that are not combined into one cell type.
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata, color=['leiden_scanvi_4.0'], wspace=0.5,ncols=1,
           legend_fontsize=8,size=3,
          save='_human_cochlea_scanvi4.pdf') #legend_loc='on data',


# In[75]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata, color=['cell_type'], wspace=0.5,ncols=1,
           legend_fontsize=2,legend_loc='on data',size=2,
          save='_human_cochlea_scanvi5.pdf') #legend_loc='on data',


# In[73]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata, color=['cell_type'], wspace=0.5,ncols=1,
           legend_fontsize=8,size=2,
          save='_human_cochlea_scanvi6.pdf') #legend_loc='on data',


# In[42]:


sc.tl.leiden(adata,resolution=3,key_added = "leiden_scanvi_3.0")
sc.pl.umap(adata, color=['cell_type','gw','leiden_scanvi_3.0'], wspace=0.5,ncols=1)


# In[43]:


# remove clusters that are not combined into one cell type.
df = (
    adata.obs.groupby(["cell_type", 'leiden_scanvi_3.0'])
    .size()
    .unstack(fill_value=0)
)
df
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]
conf_mat


# In[44]:


#Save latent space
pd.DataFrame(df).to_csv(data_dir + 'cell_porprotion.csv')


# In[45]:


import torch
import seaborn as sns
sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")
#save_dir = tempfile.TemporaryDirectory()

get_ipython().run_line_magic('config', 'InlineBackend.print_figure_kwargs={"facecolor": "w"}')
get_ipython().run_line_magic('config', 'InlineBackend.figure_format="retina"')
plt.figure(figsize=(16, 10))
_ = plt.pcolor(conf_mat)
_ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
_ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xlabel("leiden clusters")
plt.ylabel("cell_type")
plt.savefig('cochlea_scanvi.jpeg', dpi=300)


# In[46]:


#leiden cluster (3.0) that to be removed
exclude_clusters = ['21','45','46','49','58','71']
adata= adata[~adata.obs['leiden_scanvi_3.0'].isin(exclude_clusters), :]
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata, color=['leiden_scanvi_3.0','scrublet_score'],size=4,legend_fontsize =4)#legend_loc='on data'


# In[47]:


import scvi


# In[48]:


adata


# In[49]:


sc.pp.neighbors(adata, use_rep="X_scANVI",n_neighbors=30)#n_neighbors=50
sc.tl.leiden(adata,resolution=4,key_added = "leiden_scanvi_4.0")
sc.tl.umap(adata,min_dist=0.8)#0.4
sc.pl.umap(adata, color=['cell_type','gw','leiden_scanvi_4.0'], wspace=0.5,ncols=1)


# In[50]:


# remove clusters that are not combined into one cell type.
df = (
    adata.obs.groupby(["cell_type", 'leiden_scanvi_4.0'])
    .size()
    .unstack(fill_value=0)
)
df
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]
conf_mat


# In[51]:


#Save latent space
pd.DataFrame(df).to_csv(data_dir + 'cell_porprotion2.csv')


# In[55]:


#leiden cluster (4.0) that to be removed
exclude_clusters = ['47','59','86','66','82','61']
adata= adata[~adata.obs['leiden_scanvi_4.0'].isin(exclude_clusters), :]
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata, color=['leiden_scanvi_4.0','scrublet_score'],size=4,legend_fontsize =4)#legend_loc='on data'


# In[56]:


sc.pp.neighbors(adata, use_rep="X_scANVI",n_neighbors=30)#n_neighbors=50
sc.tl.leiden(adata,resolution=4,key_added = "leiden_scanvi_4.0")
sc.tl.umap(adata,min_dist=0.8)#0.4
sc.pl.umap(adata, color=['cell_type','gw','leiden_scanvi_4.0'], wspace=0.5,ncols=1)


# In[64]:


sc.tl.umap(adata,min_dist=0.9)#0.4
sc.pl.umap(adata, color=['cell_type','gw','leiden_scanvi_4.0'], wspace=0.5,ncols=1)


# In[59]:


# remove clusters that are not combined into one cell type.
df = (
    adata.obs.groupby(["cell_type", 'leiden_scanvi_4.0'])
    .size()
    .unstack(fill_value=0)
)
df
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]
conf_mat
#Save latent space
pd.DataFrame(df).to_csv(data_dir + 'cell_porprotion3.csv')


# # Anotation manually

# In[81]:


cluster_annotation = {
 'Basal_stria': 'Basal_stria',
    'Basal_stria/Spindle/Root_cells': 'Basal_stria',
    'CCs': 'CCs',
    'Chondrocytes': 'Chondrocytes',
    'Claudius/Inner-Outer_sulcus_cells': 'Claudius/Inner-Outer_sulcus_cells',
    'CoE_Lateral': 'CoE_Lateral',
    'CoE_Medial': 'CoE_Medial',
    'CoE_Roof_cells': 'CoE_Roof_cells',
    'CoE_prosensory': 'CoE_prosensory',
    'Deiters_cells': 'Deiters_cells',
    'Deiters_cells/Pillar_cells': 'Deiters_cells',
    'Endothelial': 'Endothelial',
    'Erythrocytes': 'Erythrocytes',
    'Fibrocytes': 'Fibrocytes',
    'GCs': 'GCs',
    'HCs': 'Cochlear_HCs',
    'Inner_border-phalangeal/Hensen_cells': 'Inner_border-phalangeal/Hensen_cells',
    'Interdental_cells': 'Interdental_cells',
    'Intermediate_stria': 'Intermediate_stria',
    'Macrophages': 'Macrophages',
    'Marginal_stria': 'Marginal_stria',
    'Melanocytes': 'Melanocytes',
    'Mesenchymal': 'Mesenchymal',
    'Osteoblasts': 'Osteoblasts',
    'Osteocytes': 'Osteocytes',
    'Osteocytes/Osteoblasts': 'Osteocytes',
    'Pillar_cells': 'Pillar_cells',
    'Reissner_membrane': 'Reissner_membrane',
    'SCs': 'SCs',
    'SGNs': 'SGNs',
    'SS': 'SS',
    'Smooth_Muscle_Cells': 'Smooth_Muscle_Cells',
    'Spindle/Root_cells': 'Spindle/Root_cells',
    'Tympanic_border_cells': 'Tympanic_border_cells',
    'VGNs': 'VGNs',
    'Vestibular_Dark_cells': 'Vestibular_Dark_cells',
    'Vestibular_Epithelial_cells': 'Vestibular_Epithelial_cells',
    'Vestibular_HCs': 'Vestibular_HCs',
    'Vestibular_Roof_cells': 'Vestibular_Roof_cells',
    'Vestibular_SCs': 'Vestibular_SCs',
    'pericytes': 'Pericytes',
    'pre-Osteoblasts': 'pre-Osteoblasts',
    'pre-Osteocytes': 'pre-Osteocytes',
    

    


}
adata.obs['cell_type_final'] = adata.obs['cell_type'].map(cluster_annotation).astype('category')


# In[82]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata, color=['cell_type_final'], wspace=0.5,ncols=1,
           legend_fontsize=8,
          save='_human_cochlea_scanvi_cell_type_final.pdf') #legend_loc='on data',


# In[92]:


adata.uns['cell_type_final_colors']=['#FFA500','#8A2BE2','#008080','#008B8B','#483D8B','#D2691E','#006400','#800000','#1E90FF','#F5F5DC',
                                  '#4169E1','#8FBC8F','#AFEEEE',
'#B0C4DE','#87CEEB','#7B68EE','#2E8B57','#3CB371','#CD5C5C','#556B2F','#468274','#20B2AA','#808000','#E6E6FA','#BC8F8F','#d62728',
'#E9967A','#00BEFE','#32CD32','#FFA07A','#DCDCDC','#BDB76B','#5F9EA0','#A52A2A','#A0522D','#F0F8F8',
                                 '#e2a7cc','#5858d6','#42a5e8','#6bddbc',]


# In[93]:


adata.uns['cell_type_final_colors']


# In[95]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata, color=['cell_type_final'], wspace=0.5,ncols=1,
           legend_fontsize=8,
          save='_human_cochlea_scanvi_cell_type_final.pdf') #legend_loc='on data',


# In[97]:


adata.obs['cell_type_final'].value_counts()


# In[107]:


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

save='_human_cochlea_scanvi_cell_type_final.pdf',
    size=3,


    
)
sc.pl.umap(
    adata,legend_fontsize=3,
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

save='_human_cochlea_scanvi_cell_type_final2.pdf',
    size=3,
legend_loc='on data',

    
)


# In[110]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata,color=['EPCAM', 'OTOGL',"USH1C", 'SOX2','RORB','ISL1','LGR5','OC90','PRRX1','GATA3','FGFR3','TECTA','MYO7A',
                       'STRC','TUBB3','ESRRG','MPZ','PECAM1','PTPRC','MLANA','ACAN','MKI67',]
          ,ncols=4,color_map = "inferno",size=4,
           save="_human_cochlea_scanvi_subclasses.pdf")


# In[26]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(
    adata,
    color=['gw', 'gender',  'scrublet_score',  'percent_mito', ], color_map = "PiYG",
    frameon=True, ncols = 2, wspace = 0.3, save="_human_cochlea_scanvi_gw_gender.pdf"
)#legend_loc="on data",


# In[ ]:


sc.set_figure_params(figsize=(8,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata, ['scrublet_score'], 
             
             groupby='sample',rotation=90)


# In[123]:


import omicverse as ov


# In[130]:


import scvelo as scv
import pertpy as pt


# In[131]:


adata


# In[5]:


adata.write('human_cochlea_annotation.h5ad')
adata


# In[5]:


adata=sc.read('human_cochlea_annotation.h5ad')
adata


# In[6]:


adata.raw=adata


# In[9]:


sc.pl.umap(adata, color=['cell_type_final','gw','LGR5',], wspace=0.65,ncols=2)


# In[10]:


adata.obs.index


# In[7]:


print('Total number of cells: {:d}'.format(adata.n_obs))
print('Total number of genes: {:d}'.format(adata.n_vars))
adata.obs['sample'].values.describe()
pd.DataFrame(adata.obs).to_csv(str(data_dir)+'/metadata_annotation.csv')


# In[8]:


adata.obs['cell_type_final'].values.describe()


# In[36]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(
    adata,
    color=['gw', 'gender',  'scrublet_score',  'percent_mito','sample' ], color_map = "PiYG",
    frameon=True, ncols = 2, wspace = 0.3, save="_human_cochlea_scanvi_gw_gender.pdf"
)#legend_loc="on data",


# In[38]:


sc.set_figure_params(figsize=(8,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata, ['scrublet_score'], save="_human_cochlea_scanvi_scrublet_score.pdf",
             
             groupby='sample',rotation=90)


# In[39]:


sc.set_figure_params(figsize=(8,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata, ['n_genes_by_counts'], save="_human_cochlea_scanvi_genes.pdf",
             
             groupby='sample',rotation=90)


# In[40]:


sc.set_figure_params(figsize=(8,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata, ['total_counts'], save="_human_cochlea_scanvi_counts.pdf",
             
             groupby='sample',rotation=90)


# In[11]:


adata_integration=adata
adata_integration


# In[12]:


# Benjamini-Hochberg and Bonferroni FDR helper functions.

def bh(pvalues):
    """
    Computes the Benjamini-Hochberg FDR correction.
    
    Input:
        * pvals - vector of p-values to correct
    """
    pvalues = np.array(pvalues)
    n = int(pvalues.shape[0])
    new_pvalues = np.empty(n)
    values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
    values.sort()
    values.reverse()
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n/rank) * pvalue)
    for i in range(0, int(n)-1):
        if new_values[i] < new_values[i+1]:
            new_values[i+1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return new_pvalues

def bonf(pvalues):
    """
    Computes the Bonferroni FDR correction.
    
    Input:
        * pvals - vector of p-values to correct
    """
    new_pvalues = np.array(pvalues) * len(pvalues)
    new_pvalues[new_pvalues>1] = 1
    return new_pvalues


# In[13]:


meta['sample'] = meta['sample_id'].astype('str')
plotmeta = list(meta.columns)
plotmeta.append('sample')
print('Number of sample: ', meta.index.size)
meta


# In[14]:


holder = []
for sample in meta.index:
    print(sample)
    # Load 10x data as AnnData
    #adata_sample = sc.read_10x_mtx(data_dir+'rawdata_pbx3_neo/'+sample+'/outs/filtered_feature_bc_matrix/',cache=True)
    holder.append(sc.read(data_dir+sample+'/human_gw_raw.h5ad',cache=True)) 
    # Set names of observation as sample + _ + barcode/probe
    holder[-1].obs_names = [sample+'_'+i.split('-')[0] for i in holder[-1].obs_names]
    # Filer genes expressed in less than 3 cells
    sc.pp.filter_genes(holder[-1], min_cells=3)
    # Filer cells with less than 10 genes expressed
    sc.pp.filter_cells(holder[-1], min_genes=10)
    # add in metadata
    holder[-1].obs['sample'] = sample
    for val in meta.columns:
        holder[-1].obs[val] = meta[val][sample]
    # Extract mitochondial genes
    mito_genes = [name for name in holder[-1].var_names if name.startswith('MT-')]
    #for each cell compute fraction of counts in mito genes vs. all genes
    #the `.A1` is only necessary, as X is sparse - it transform to a dense array after summing
    holder[-1].obs['percent_mito'] = np.sum(
        holder[-1][:, mito_genes].X, axis=1).A1 / np.sum(holder[-1].X, axis=1).A1
    #add the total counts per cell as observations-annotation to adata
    holder[-1].obs['n_counts'] = holder[-1].X.sum(axis=1).A1
    print('Total number of cells: {:d}'.format(holder[-1].n_obs))
    print('Total number of genes: {:d}'.format(holder[-1].n_vars))


# In[15]:


# confirm N samples
print(len(holder))
# merge datasets
adata = holder[0].concatenate(holder[1:],join='outer',index_unique=None)
# copy of this matrix in Compressed Sparse Row format
adata.X = adata.X.tocsr()
adata


# In[16]:


adata.obs.index


# In[17]:


adata_integration.obs.index


# In[18]:


adata_new=adata[adata_integration.obs.index]


# In[19]:


adata_new


# In[22]:


adata_integration


# In[24]:


#adata_new.layers['counts']=adata_integration.layers['counts']
adata_new.obsp['connectivities']=adata_integration.obsp['connectivities']
adata_new.obsp['distances']=adata_integration.obsp['distances']

adata_new.obsm['X_pca']=adata_integration.obsm['X_pca']
adata_new.obsm['X_scANVI']=adata_integration.obsm['X_scANVI']
adata_new.obsm['X_scVI']=adata_integration.obsm['X_scVI']
adata_new.obsm['X_umap']=adata_integration.obsm['X_umap']
adata_new.obsm['_scvi_extra_categorical_covs']=adata_integration.obsm['_scvi_extra_categorical_covs']
adata_new.obsm['_scvi_extra_continuous_covs']=adata_integration.obsm['_scvi_extra_continuous_covs']

adata_new.uns['_scvi_manager_uuid']=adata_integration.uns['_scvi_manager_uuid']
adata_new.uns['_scvi_uuid']=adata_integration.uns['_scvi_uuid']
adata_new.uns['annotation_colors']=adata_integration.uns['annotation_colors']
adata_new.uns['ccgs']=adata_integration.uns['ccgs']
adata_new.uns['cell_type_final_colors']=adata_integration.uns['cell_type_final_colors']
adata_new.uns['hvg']=adata_integration.uns['hvg']
adata_new.uns['leiden']=adata_integration.uns['leiden']
adata_new.uns['leiden_R']=adata_integration.uns['leiden_R']
adata_new.uns['leiden_R_colors']=adata_integration.uns['leiden_R_colors']
adata_new.uns['leiden_res2.0_colors']=adata_integration.uns['leiden_res2.0_colors']
adata_new.uns['leiden_res2.0']=adata_integration.uns['leiden_res2.0']
adata_new.uns['leiden_scanvi_2.0_colors']=adata_integration.uns['leiden_scanvi_2.0_colors']
adata_new.uns['log1p']=adata_integration.uns['log1p']
adata_new.uns['neighbors']=adata_integration.uns['neighbors']
adata_new.uns['sample_colors']=adata_integration.uns['sample_colors']
adata_new.uns['umap']=adata_integration.uns['umap']
adata_new.uns['leiden_scanvi_2.0']=adata_integration.uns['leiden_scanvi_2.0']
adata_new.var['highly_variable']=adata_integration.var['highly_variable']
adata_new.var['highly_variable_rank']=adata_integration.var['highly_variable_rank']
adata_new.var['means']=adata_integration.var['means']
adata_new.var['variances']=adata_integration.var['variances']
adata_new.var['variances_norm']=adata_integration.var['variances_norm']

adata_new.obs['cell_type_final']=adata_integration.obs['cell_type_final']
adata_new.obs['low_ncounts']=adata_integration.obs['low_ncounts']
adata_new.obs['low_ncounts_high_mito']=adata_integration.obs['low_ncounts_high_mito']
adata_new.obs['S_score']=adata_integration.obs['S_score']
adata_new.obs['G2M_score']=adata_integration.obs['G2M_score']
adata_new.obs['phase']=adata_integration.obs['phase']
adata_new.obs['scrublet_score']=adata_integration.obs['scrublet_score']
adata_new.obs['scrublet_cluster_score']=adata_integration.obs['scrublet_cluster_score']
adata_new.obs['zscore']=adata_integration.obs['zscore']
adata_new.obs['bh_pval']=adata_integration.obs['bh_pval']
adata_new.obs['bonf_pval']=adata_integration.obs['bonf_pval']
adata_new.obs['is_doublet']=adata_integration.obs['is_doublet']
adata_new.obs['is_doublet_scrub']=adata_integration.obs['is_doublet_scrub']
adata_new.obs['_scvi_batch']=adata_integration.obs['_scvi_batch']
adata_new.obs['_scvi_labels']=adata_integration.obs['_scvi_labels']
adata_new.obs['leiden_scanvi_2.0']=adata_integration.obs['leiden_scanvi_2.0']
adata_new.obs['leiden_res2.0']=adata_integration.obs['leiden_res2.0']


# In[26]:


adata_new.layers["counts"] = adata_new.X.copy() # preserve counts


# In[28]:


sc.pl.umap(adata_new, color=['cell_type_final','gw','LGR5','RORB'], wspace=0.65,ncols=2)


# In[30]:


small_marker_dict={
      'Epithelium':[  "EPCAM",],
 'Vestibular supporting cells/ epithelial cells':[ "MEIS2","ADAMTSL1","OTOGL","USH1C",],
  'Cochlear duct floor medial':["TECTA","FGF10","JAG1",],
    'Cochlear duct floor lateral':["GATA3","FGFR3","PROX1","BMP4",],
 'Cochlear duct floor prosensory':[ "RORB","ISL1","LGR5","SOX2","FGF20",],
    'Vestibular roof cells':[ "NTN1","SMOC2","WNT3",],
    'Cochlear roof cells':[ "OTX2",	"FGF9",	"WNT4",	"GSC",],
   'Vestibular Dark cells':["KCNE1","ATP1B2","SPP1",],
    'Vestibular hair cells':[ "STRC","OTOF","USH2A","MYO15A",],
 'Cycling cells' :[ "MKI67","TOP2A","HMGB2",],
 'Mesenchymal':["PRRX1",],
'Chondrocytes' :["ACAN"],
 'Fibrocytes':['OTOS','CAR3','COL9A2','COL9A3'],
    'Mesenchymal':["PRRX1",], 
    'Tympanic border cells':['EMILIN2','NOTUM','RARRES1'],
 'pre-Osteoblasts':["DLX5","RUNX2",],
    'Osteoblasts':["BGLAP", "IBSP","IFITM5"],
      'Osteocytes':["DMP1","PHEX",],
     'Marginal stria':['KCNE1','ESRRB','DCLK1'],  
     'Basal stria':['CLDN11','ATP6V0A4','TJP1'],
'Spindle/Root_cells':['SLC26A4','ANXA1'],
     'Reissner membrane':["SLC26A7"],
    'Deiters_cells':['LGR5','FGFR3','PROX1','CEACAM16'],
     'Pillar_cells' :['SMPX','LGR6','ENAH'],
  'Interdental_cells':['OTOA'],
     'Pan_Supporting_cells':['OTOG','OTOGl','USH1C','GATA3'],
    'Claudius/Inner-Outer_sulcus_cells':['APOE','EPYC'],
    'Inner_border-phalangeal/Hensen_cells':['S100A1','SLC1A3'],
     'Scala vestibuli border cells':['FXYD2'],
  
     'Endothelial cells':["PECAM1",'VWF','ESAM'],
     'Smooth Muscle Cells':["TAGLN"],
     'pericytes':['RGS5'],
   
    'macrophages':["PTPRC",'AIF1','CD163'],
     'Erythrocytes':["RHD"],
     'Glial Precursor Cells':["OLIG1","OLIG2"],
    'Glial.cells':['MAG','MOG','MOBP'],
    'Schwann.Cells':['MPZ','PRX',"PLP1","PMP22",],
 'Intermediate stria':['TYR','DCT'],
     'Melanocytes' :[ "MLANA"],
     'Neutrophils':["LY6G"],
    'NK Cells':["NKG7"],
    'Monocytes':["CD68"],
'Mast Cells':["PRSS34"],
    'B Cells':["CD19"],
    'T Cells':["CD3G"],
     'Surrounding structures':["OSR2","BMP6","SLC7A11","CHRDL1","COL1A1","ALDH1A2"],
    'VGN':[ "TLX3",],
     'SGN' :["EPHA5",'SNAP25','NEFL','PVALB','CALB2','PRPH','ANO2'],
    'Cochlear hair cell':['MYO7A','SLC26A5','CALB1','POU4F3','TMC1','OTOF','SLC17A8',],

}
# check if the markers are in the data
smarker_genes_in_data = dict()
for ct, markers in small_marker_dict.items():
    markers_found = list()
    for marker in markers:
        if marker in adata_new.var.index:
            markers_found.append(marker)
    smarker_genes_in_data[ct] = markers_found
#del [] # remove the last marker
del_markers = list()
for ct, markers in smarker_genes_in_data.items():
    if markers==[]:
        del_markers.append(ct)
for ct in del_markers:
    del smarker_genes_in_data[ct]


# In[31]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
#sc.tl.dendrogram(adata, groupby='cell_type_final')
sc.pl.dotplot(
    adata_new,
    groupby="cell_type_final",
    var_names=smarker_genes_in_data,
  # dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
     save="human_fetal_cochlea_annotation_markers.pdf"
)


# In[32]:


adata_new.obs['cell_type_final'].value_counts()


# In[68]:


small_marker_dict={
 'Cochlear duct floor lateral':["GATA3",],
     'Cochlear duct floor medial':["TECTA",],
    'Cochlear roof cells':[ "OTX2","OC90",	"FGF9",	"WNT4",	"GSC",],
    'Cochlear duct floor prosensory':[ "RORB","ISL1","LGR5",],
'Hair cell':['MYO7A',"STRC",],
    'Pan_Supporting_cells':['OTOG','OTOGL'],
'Pillar_cells' :['LGR6',],
   'Deiters_cells':['FGFR3','PROX1'],
 'Interdental_cells':['OTOA'],
     'Claudius/Inner-Outer_sulcus_cells':['EPYC'],
     'Inner_border-phalangeal/Hensen_cells':['S100A1','SLC1A3'],
   'Neurons' :["EPHA5",'SNAP25','NEFL',], 
     'Glial.Cells':['MPZ','PRX',"PLP1","PMP22",],
 'Vestibular Dark cells':["SPP1",],
 'Vestibular epithelial cells':[ "MEIS2","ADAMTSL1",],
    
'Vestibular roof cells':[ "NTN1","SMOC2","WNT3",],
     'Vestibular supporting cells':[ "USH1C",],
'Cycling cells' :[ "MKI67","TOP2A","HMGB2",],
'Erythrocytes':["HBG2",],
    'macrophages':["PTPRC",'AIF1','CD163'], 
    'Tympanic border cells':['EMILIN2','NOTUM','RARRES1'],
'Reissner membrane':["SLC26A7"],
 'pre-Osteoblasts':["DLX5","RUNX2",],
    'Osteoblasts':[ "IBSP","IFITM5"],
      'Osteocytes':["DMP1","PHEX",],
'Surrounding structures':["OSR2",],
'Endothelial cells':["PECAM1",'VWF','ESAM'],
   'Basal stria':['CLDN11','ATP6V0A4','TJP1'], 
    'Marginal stria':['KCNE1','ESRRB',], 
    'Intermediate stria':['TYR','DCT'],
'Fibrocytes':['OTOS','CAR3','COL9A2','COL9A3'],
 'pericytes':['RGS5'],
    'Spindle/Root_cells':['SLC26A4',],
     'Smooth Muscle Cells':["TAGLN"],
   'Chondrocytes' :["ACAN"],
    'Melanocytes' :[ "MLANA"],
     'Mesenchymal':["PRRX1",], 
    
  


}
# check if the markers are in the data
smarker_genes_in_data = dict()
for ct, markers in small_marker_dict.items():
    markers_found = list()
    for marker in markers:
        if marker in adata_new.var.index:
            markers_found.append(marker)
    smarker_genes_in_data[ct] = markers_found
#del [] # remove the last marker
del_markers = list()
for ct, markers in smarker_genes_in_data.items():
    if markers==[]:
        del_markers.append(ct)
for ct in del_markers:
    del smarker_genes_in_data[ct]


# In[34]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
#sc.tl.dendrogram(adata, groupby='cell_type_final')
sc.pl.dotplot(
    adata_new,
    groupby="cell_type_final",
    var_names=smarker_genes_in_data,
  # dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
     save="human_fetal_cochlea_annotation_markers.pdf"
)


# In[69]:


adata_new.obs['cell_type_final'].cat.reorder_categories(['CoE_Lateral','CoE_Medial','CoE_Roof_cells','CoE_prosensory','Cochlear_HCs',
                                                     'Vestibular_HCs',
                                                     'SCs','Pillar_cells','Deiters_cells','Interdental_cells','Claudius/Inner-Outer_sulcus_cells',
                                                     'Inner_border-phalangeal/Hensen_cells','SGNs','VGNs','GCs',
                                                     'Vestibular_Dark_cells','Vestibular_Epithelial_cells','Vestibular_Roof_cells','Vestibular_SCs',
                                                     'CCs','Erythrocytes','Macrophages',
                                                     'Tympanic_border_cells','Reissner_membrane','pre-Osteoblasts','pre-Osteocytes',
                                                     'Osteoblasts','Osteocytes',
                                                     'SS',
                                                     'Endothelial','Basal_stria',
                                                     'Marginal_stria', 'Intermediate_stria',
                                                     'Fibrocytes','Pericytes','Spindle/Root_cells','Smooth_Muscle_Cells',
                                                     'Chondrocytes', 'Melanocytes', 'Mesenchymal',
           
                                                         ], 
                                                    inplace=True,)


# In[70]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
#sc.tl.dendrogram(adata, groupby='cell_type_final')
sc.pl.dotplot(
    adata_new,
    groupby="cell_type_final",
    var_names=smarker_genes_in_data,
  # dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
     save="human_fetal_cochlea_annotation_markers.pdf"
)


# In[39]:


adata_new.write('human_cochlea_annotation_add_40000genes.h5ad',compression='gzip')


# In[40]:


adata_new.obs.head()


# In[67]:


adata_new


# In[62]:


#we first want to determine whether opposing wnt/shh signaling pathway exist in early first trimester stage.
#as support for recent CSC paper of IEO.
#subset cochlear epithelium cell-type cells and vestibular epithelium cell-type cells
adata_coch_vesb=adata_new[adata_new.obs['cell_type_final'].isin(['CoE_Roof_cells','CoE_Lateral','CoE_Medial','CoE_prosensory','Cochlear_HCs',
                                                                 'SCs','Deiters_cells','Pillar_cells','Interdental_cells','Inner_border-phalangeal/Hensen_cells',
                                                                 'Claudius/Inner-Outer_sulcus_cells',
                                                                 'Vestibular_Dark_cells','Vestibular_Epithelial_cells','Vestibular_HCs',
                                                                 'Vestibular_Roof_cells','Vestibular_SCs',
                                                                
                                                                
                                                                ]
    
)]
adata_coch_vesb


# In[63]:


sc.pp.normalize_total(adata_coch_vesb, target_sum=1e4)
sc.pp.log1p(adata_coch_vesb)


# In[61]:


sc.pp.regress_out(adata_coch_vesb, ['total_counts'])


# In[64]:


sc.pp.scale(adata_coch_vesb, max_value=10)


# In[65]:


adata_coch_vesb.obs['cell_type_final'].cat.reorder_categories(['CoE_prosensory','CoE_Lateral','CoE_Medial','CoE_Roof_cells','Cochlear_HCs',
                                                                 'SCs','Deiters_cells','Pillar_cells','Interdental_cells',
                                                               'Inner_border-phalangeal/Hensen_cells',
                                                                 'Claudius/Inner-Outer_sulcus_cells',
                                                                 'Vestibular_HCs','Vestibular_Epithelial_cells','Vestibular_SCs',
                                                                 'Vestibular_Roof_cells','Vestibular_Dark_cells',
                                                    
           
                                                         ], 
                                                   ) # inplace=True,


# In[66]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
#sc.tl.dendrogram(adata, groupby='cell_type_final')
sc.pl.dotplot(
    adata_coch_vesb,
    groupby="cell_type_final",
    var_names=[ 
                               'DLX5','MSX1','GPR155','UBE2C','PCDH20','NEUROD6','WNT3A',  'HMX3',   #VESTIBULAR/DORSAL MARKERS,
                                'NR2F2', 'NR2F1','GATA3','INSM1','HES6','TMPRSS3','FGFR3','LGR5',                # COCHLEAR/VENTRAL MARKERS,
                                'SULF1',  'LRP2', 'GAS1', 'PTCH1',              # SHH SIGNALING
                                 'ATOH1','CCER2','KCNH6','GRXCR2','MYO7A','LHX3','POU4F3',  # HAIR CELLS
                                'EPCAM','FBXO2', # OTIC MARKERS,
                                
                            ],
  # dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    save='_markers_cochlear_vestibular_epithelium.pdf'
)


# In[ ]:




