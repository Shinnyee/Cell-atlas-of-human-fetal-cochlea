#!/usr/bin/env python
# coding: utf-8

# In[2]:


cd F:\PROJECTS\PROJECT_HUMAN_FETAL_COCHLEAE\WORKPLACE\R\


# this dataset is from van der valk et al. Cell Report,2023.

# In[3]:


import scanpy as sc
sc.set_figure_params()
adata = sc.read_10x_mtx('./FW7_5')
print(adata)


# In[3]:


adata.obs['batch']='gw9.5'


# In[4]:


adata.obs_names_make_unique()
adata.var_names_make_unique()


# In[5]:


adata=sc.AnnData(adata.X,obs=adata.obs,var=adata.var)
adata.var["Gene"]=adata.var_names
adata.obs["CellID"]=adata.obs_names
adata


# In[6]:


adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# In[7]:


sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.3, multi_panel=True)


# In[8]:


sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


# In[9]:


adata_filtered = adata[adata.obs.n_genes_by_counts < 5000, :]
adata_filtered = adata_filtered[adata_filtered.obs.total_counts > 500, :]
adata_filtered = adata_filtered[adata_filtered.obs.pct_counts_mt < 5, :]
adata_filtered


# In[10]:


sc.pl.scatter(adata_filtered, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata_filtered, x='total_counts', y='n_genes_by_counts')


# In[11]:


sc.pp.normalize_total(adata_filtered, target_sum=1e4)
sc.pp.log1p(adata_filtered)
sc.pp.highly_variable_genes(adata_filtered, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata_filtered)


# In[12]:


adata_filtered.raw = adata_filtered


# In[13]:


adata_filtered = adata_filtered[:, adata_filtered.var.highly_variable]
adata_filtered


# In[14]:


sc.pp.regress_out(adata_filtered, ['total_counts'])
sc.pp.scale(adata_filtered, max_value=10)
sc.tl.pca(adata_filtered, svd_solver='arpack',n_comps=100)
sc.pl.pca(adata_filtered, color=['MYO7A','RORB','ATOH1','ISL1', 'SOX2','STRC','OTOF','USH2A' ] )


# In[2]:


import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys
import scipy
import anndata


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams["font.family"] = "Arial"

sc.set_figure_params(figsize=(4, 4))
sc.settings.set_figure_params(dpi = 150, color_map = 'RdPu', dpi_save = 600, vector_friendly = True, format = 'pdf')
palette = ['#fbbc04', '#199de5', '#cfe6d5']


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


# In[17]:


adata=adata_filtered


# In[18]:


adata


# In[19]:


adata=adata.raw.to_adata()
adata


# In[20]:


sc.pp.neighbors(adata, n_neighbors=50, n_pcs=100)
sc.tl.leiden(adata,resolution=2)
sc.tl.louvain(adata)
sc.tl.umap(adata,min_dist=0.4)


# In[22]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata, color=[
              "RORB","ISL1","LGR5","SOX2","FGF20","ATOH1","CDKN1B",   # COCHLEAR DUCT FLOOR PROSENSORY
                   "TECTA","FGF10","JAG1",#FLOOR MEDIAL
                   "GATA3","FGFR3","PROX1","BMP4",   #FLOOR LATERAL
                   "OTX2",	"FGF9",	"WNT4",	"GSC",#COCHLEAR ROOF CELLS
                   "MEIS2",	"ADAMTSL1",	"OTOGL","USH1C",#VESTIBULAR SIPPORTING CELLS/EPITHELIAL CELLS
                   "NTN1","SMOC2","WNT3",#VESTIBULAR ROOF CELLS
                   "KCNE1","ATP1B2","SPP1",#DARK CELLS
                   "STRC","OTOF","USH2A","MYO15A","SLC26A5","SLC17A8","SLC7A14",  #HAIR CELLS
                   "SNAP25","TUBB3","PRPH","CALB1","PBX3","ESRRG", # SGN
                   "MBP","MPZ","PLP1","PMP22", # GC
                   "EPCAM", #EPITHELIUM
                   "PRRX1",#MESENCHYMAL
                   "ACAN", #CHONDROCYTES
                   "PECAM1", #ENDOTHELIAL
                   "PTPRC", #MACROPHAGES
                   "MLANA",  #MELANOCYTES
                   "MKI67","TOP2A","HMGB2", #CYCLING               
        
                          'leiden'])


# In[27]:


small_marker_dict={
    'Cochlear epithelium':[  "EPCAM",],
    'Cochlear duct floor prosensory':[ "RORB","ISL1","LGR5","SOX2","FGF20",],
     'Cochlear duct floor medial':["TECTA","FGF10","JAG1",],
    'Cochlear duct floor lateral':["GATA3","FGFR3","PROX1","BMP4",],
    'Cochlear roof cells':[ "OTX2",	"FGF9",	"WNT4",	"GSC",],
 
 'Vestibular supporting cells/ epithelial cells':[ "MEIS2","ADAMTSL1","OTOGL","USH1C",],
     'Vestibular roof cells':[ "NTN1","SMOC2","WNT3",],
    'Vestibular Dark cells':["KCNE1","ATP1B2","SPP1",],
    'Vestibular hair cells':[ "STRC","OTOF","USH2A","MYO15A",],
    
'NEURONAL':[ "SNAP25","TUBB3","PRPH","CALB1","PBX3","ESRRG",],
    'SGN':[ "EPHA5",],
    'VGN':[ "TLX3","SALL3"],
'GLIAL':[ "MBP","MPZ","PLP1","PMP22",],
 'MESENCHYMAL':["PRRX1",],
'ENDOTHELIAL' :[ "PECAM1"],
'MACROPHAGES' :[ "PTPRC"],
'MELANOCYTES' :[ "MLANA"],
'CHONDROCYTES' :[ "ACAN"],
 'CYCLING' :[ "MKI67","TOP2A","HMGB2",],
       
}
# check if the markers are in the data
smarker_genes_in_data = dict()
for ct, markers in small_marker_dict.items():
    markers_found = list()
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    smarker_genes_in_data[ct] = markers_found
#del [] # remove the last marker
del_markers = list()
for ct, markers in smarker_genes_in_data.items():
    if markers==[]:
        del_markers.append(ct)
for ct in del_markers:
    del smarker_genes_in_data[ct]


# In[28]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.dendrogram(adata,
    groupby="leiden",)
sc.pl.dotplot(
    adata,
    groupby="leiden",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    save="_cochlea-cell_type_gw9_5_FIRST_ROUND_SCREENING_1.pdf"
)


# In[26]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata,color=["leiden",],legend_fontsize =8,legend_loc="on data",)


# In[42]:


# 1st remove ambigous clusrer then subdivide particular clusters
exclude_clusters = ['0','1','4','5','7','8','11','13','14','15','19','22','23','26','28','31',]
adata_rm= adata[~adata.obs['leiden'].isin(exclude_clusters), :]
adata_rm


# In[43]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_rm,color=["leiden",],legend_fontsize =8,legend_loc="on data",)


# In[44]:


cluster_annotation = {
    '16': 'Mesenchymal',
    '10': 'Mesenchymal',
 '25': 'CoE_prosensory',
 '24': 'CoE_Medial',
 '6': 'CoE_Lateral',
    '12': 'GCs', 
     '21': 'SGNs',
    '9': 'SGNs',
'29': 'VGNs',
 '3': 'Chondrocytes',
  '32': 'Melanocytes',
 '27': 'Macrophages',
'20': 'Endothelial',
 '2': 'CCs',
     '30': 'Vestibular_HCs',
   '17': 'Vestibular_SCs',
  '18': 'Vestibular_Roof_cells',
   
}
adata_rm.obs['cell_type'] = adata_rm.obs['leiden'].map(cluster_annotation).astype('category')


# In[45]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_rm,color=["leiden","cell_type"],legend_loc="on data",legend_fontsize =10)


# In[64]:


sc.pp.neighbors(adata_rm, n_neighbors=30, n_pcs=100)
sc.tl.umap(adata_rm,min_dist=0.6,)


# In[65]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_rm,color=["leiden","cell_type"],legend_fontsize =5,legend_loc="on data",)


# In[66]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.correlation_matrix(adata_rm, 'cell_type',
                         save="cell-type-correlation_HUMAN_GW9_5.pdf")


# In[67]:


adata_rm


# In[68]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_rm,color='cell_type',legend_loc='on data',frameon=False, legend_fontsize=5, legend_fontoutline=False,
           title="Human_fetal_cochlea_GW9_5 n=4,284 nuclei",
          save="_HUMAN_GW9_5_CELL_TYPE_ANNO.pdf")
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_rm,color='cell_type',frameon=False, legend_fontsize=8, legend_fontoutline=False,
           title="Human_fetal_cochlea_GW9_5 n=4,284 nuclei",
          save="_HUMAN_GW9_5_CELL_TYPE_ANNO.pdf")


# In[69]:


adata_rm.write("human_gw9_5_raw.h5ad")


# In[70]:


small_marker_dict={
    'Cochlear epithelium':[  "EPCAM",],
    'Cochlear duct floor prosensory':[ "RORB","ISL1","LGR5","SOX2","FGF20",],
     'Cochlear duct floor medial':["TECTA","FGF10","JAG1",],
    'Cochlear duct floor lateral':["GATA3","FGFR3","PROX1","BMP4",],
    'Cochlear roof cells':[ "OTX2",	"FGF9",	"WNT4",	"GSC",],
 
 'Vestibular supporting cells/ epithelial cells':[ "MEIS2","ADAMTSL1","OTOGL","USH1C",],
     'Vestibular roof cells':[ "NTN1","SMOC2","WNT3",],
    'Vestibular Dark cells':["KCNE1","ATP1B2","SPP1",],
    'Vestibular hair cells':[ "STRC","OTOF","USH2A","MYO15A",],
    
'NEURONAL':[ "SNAP25","TUBB3","PRPH","CALB1","PBX3","ESRRG",],
    'SGN':[ "EPHA5",],
    'VGN':[ "TLX3","SALL3"],
'GLIAL':[ "MBP","MPZ","PLP1","PMP22",],
 'MESENCHYMAL':["PRRX1",],
'ENDOTHELIAL' :[ "PECAM1"],
'MACROPHAGES' :[ "PTPRC"],
'MELANOCYTES' :[ "MLANA"],
'CHONDROCYTES' :[ "ACAN"],
 'CYCLING' :[ "MKI67","TOP2A","HMGB2",],
       
}
# check if the markers are in the data
smarker_genes_in_data = dict()
for ct, markers in small_marker_dict.items():
    markers_found = list()
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    smarker_genes_in_data[ct] = markers_found
#del [] # remove the last marker
del_markers = list()
for ct, markers in smarker_genes_in_data.items():
    if markers==[]:
        del_markers.append(ct)
for ct in del_markers:
    del smarker_genes_in_data[ct]


# In[71]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.tl.dendrogram(adata_rm, groupby='cell_type')
sc.pl.dotplot(
    adata_rm,
    groupby="cell_type",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
     save="cochlea-cell_type_refined_human_gw9_5-order.pdf"
)


# In[74]:


small_marker_dict={
'NEURONAL':[ "SNAP25","TUBB3","PRPH","CALB1","PBX3","ESRRG",],
    'SGN':[ "EPHA5",],
    'VGN':[ "TLX3",],
    'Epithelium':[  "EPCAM",],
    'Vestibular hair cells':[ "STRC","OTOF","USH2A","MYO15A",],
     'Vestibular roof cells':[ "NTN1","SMOC2","WNT3",],
 'Vestibular supporting cells':[ "MEIS2","ADAMTSL1","OTOGL","USH1C",],
     'Cochlear duct floor prosensory':[ "RORB","ISL1","LGR5","SOX2","FGF20",],
    'Cochlear duct floor lateral':["GATA3","FGFR3","PROX1","BMP4",],
    'Cochlear duct floor medial':["TECTA","FGF10","JAG1",],
    'CHONDROCYTES' :[ "ACAN"],
    'CYCLING' :[ "MKI67","TOP2A","HMGB2",],
 'MESENCHYMAL':["PRRX1",],
'ENDOTHELIAL' :[ "PECAM1"],
    'MACROPHAGES' :[ "PTPRC"],
'GLIAL':[ "MBP","MPZ","PLP1","PMP22",],
'MELANOCYTES' :[ "MLANA"],

 
       
}
# check if the markers are in the data
smarker_genes_in_data = dict()
for ct, markers in small_marker_dict.items():
    markers_found = list()
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    smarker_genes_in_data[ct] = markers_found
#del [] # remove the last marker
del_markers = list()
for ct, markers in smarker_genes_in_data.items():
    if markers==[]:
        del_markers.append(ct)
for ct in del_markers:
    del smarker_genes_in_data[ct]


# In[75]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.tl.dendrogram(adata_rm, groupby='cell_type')
sc.pl.dotplot(
    adata_rm,
    groupby="cell_type",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
     save="cochlea-cell_type_refined_human_gw9_5-order.pdf"
)


# In[76]:


sc.set_figure_params(figsize=(6,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_rm, ['n_genes_by_counts'], 
             save="cochlea-cell_type_refined_human_gw9_5_gene_detection.pdf",
             groupby='cell_type',rotation=90)


# In[77]:


sc.set_figure_params(figsize=(6,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_rm, ['total_counts'], 
             save="cochlea-cell_type_refined_human_gw9_5_total_counts.pdf",
             groupby='cell_type',rotation=90)


# In[78]:


sc.set_figure_params(figsize=(6,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_rm, ['pct_counts_mt'], 
             save="cochlea-cell_type_refined_human_gw9_5_mt_pct.pdf",
             groupby='cell_type',rotation=90)


# In[79]:


adata_rm.write("human_gw9_5_raw.h5ad")


# In[80]:


adata_rm.obs['cell_type'].value_counts()


# In[4]:


adata_gw9=sc.read("human_gw9_5_raw.h5ad")
adata_gw9


# In[7]:


adata = sc.read_10x_mtx('./FW7_5')
print(adata)


# In[8]:


adata.X.max()


# In[11]:


adata_new=adata[adata_gw9.obs.index]
adata_new


# In[12]:


adata_new.X.max()


# In[14]:


adata_new.obsp['connectivities']=adata_gw9.obsp['connectivities']
adata_new.obsp['distances']=adata_gw9.obsp['distances']

adata_new.obsm['X_pca']=adata_gw9.obsm['X_pca']
adata_new.obsm['X_umap']=adata_gw9.obsm['X_umap']

adata_new.uns['cell_type_colors']=adata_gw9.uns['cell_type_colors']
adata_new.uns['dendrogram_cell_type']=adata_gw9.uns['dendrogram_cell_type']
adata_new.uns['dendrogram_leiden']=adata_gw9.uns['dendrogram_leiden']
adata_new.uns['hvg']=adata_gw9.uns['hvg']
adata_new.uns['leiden']=adata_gw9.uns['leiden']
adata_new.uns['leiden_colors']=adata_gw9.uns['leiden_colors']
adata_new.uns['log1p']=adata_gw9.uns['log1p']
adata_new.uns['louvain']=adata_gw9.uns['louvain']
adata_new.uns['neighbors']=adata_gw9.uns['neighbors']
adata_new.uns['pca']=adata_gw9.uns['pca']
adata_new.uns['umap']=adata_gw9.uns['umap']

adata_new.var['gene_ids']=adata_gw9.var['gene_ids']
adata_new.var['feature_types']=adata_gw9.var['feature_types']
adata_new.var['Gene']=adata_gw9.var['Gene']
adata_new.var['n_cells_by_counts']=adata_gw9.var['n_cells_by_counts']
adata_new.var['mean_counts']=adata_gw9.var['mean_counts']
adata_new.var['pct_dropout_by_counts']=adata_gw9.var['pct_dropout_by_counts']
adata_new.var['total_counts']=adata_gw9.var['total_counts']
adata_new.var['highly_variable']=adata_gw9.var['highly_variable']
adata_new.var['means']=adata_gw9.var['means']
adata_new.var['dispersions']=adata_gw9.var['dispersions']
adata_new.var['dispersions_norm']=adata_gw9.var['dispersions_norm']

adata_new.obs['batch']=adata_gw9.obs['batch']
adata_new.obs['CellID']=adata_gw9.obs['CellID']
adata_new.obs['n_genes_by_counts']=adata_gw9.obs['n_genes_by_counts']
adata_new.obs['total_counts']=adata_gw9.obs['total_counts']
adata_new.obs['total_counts_mt']=adata_gw9.obs['total_counts_mt']
adata_new.obs['pct_counts_mt']=adata_gw9.obs['pct_counts_mt']
adata_new.obs['leiden']=adata_gw9.obs['leiden']
adata_new.obs['louvain']=adata_gw9.obs['louvain']
adata_new.obs['cell_type']=adata_gw9.obs['cell_type']


# In[16]:


# Saving count data
adata_new.layers["counts"] = adata_new.X.copy()


# In[17]:


# Normalizing to median total counts
sc.pp.normalize_total(adata_new)
# Logarithmize the data
sc.pp.log1p(adata_new)


# In[18]:


adata_new.layers["logcounts"] = adata_new.X.copy()


# In[20]:


adata_new.layers["logcounts"].max()


# In[21]:


sc.pl.umap(adata_gw9, color=['cell_type', 
                                
                            ],
           size=20,
           legend_fontsize=6,ncols = 2,wspace = 0.5,
           
           
          )
sc.pl.umap(adata_new, color=['cell_type', 
                                
                            ],
           size=20,
           legend_fontsize=6,ncols = 2,wspace = 0.5,
           
           
          )


# In[22]:


adata_new.write("human_gw9_5_rawcounts.h5ad")


# In[5]:


adata_gw9.obs['cell_type'].value_counts()


# In[17]:


adata_gw9_coe=adata_gw9[adata_gw9.obs['cell_type'].isin(['CoE_Lateral','CoE_Medial',
                                                   'CoE_prosensory',]
    
)]
adata_gw9_coe


# In[18]:


adata_gw9_coe.write("human_gw9_5_raw_coe.h5ad")


# In[9]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_gw9, color=['cell_type', 
                               'DLX5','MSX1','GPR155','UBE2C','PCDH20','NEUROD6','WNT3A',  'HMX3',   #VESTIBULAR/DORSAL MARKERS,
                                'NR2F2', 'NR2F1','GATA3','INSM1','HES6','TMPRSS3','FGFR3','LGR5',                # COCHLEAR/VENTRAL MARKERS,
                                'SULF1',  'LRP2', 'GAS1', 'PTCH1',              # SHH SIGNALING
                                 'ATOH1','CCER2','KCNH6','GRXCR2','MYO7A','LHX3','POU4F3',  # HAIR CELLS
                                'EPCAM','FBXO2', # OTIC MARKERS,
                                
                            ],
           size=20,
           legend_fontsize=6,ncols = 2,wspace = 0.5,
           
           
          )


# In[10]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
#sc.tl.dendrogram(adata, groupby='cell_type_final')
sc.pl.dotplot(
    adata_gw9,
    groupby="cell_type",
    var_names=[ 
                               'DLX5','MSX1','GPR155','UBE2C','PCDH20','NEUROD6','WNT3A',  'HMX3',   #VESTIBULAR/DORSAL MARKERS,
                                'NR2F2', 'NR2F1','GATA3','INSM1','HES6','TMPRSS3','FGFR3','LGR5',                # COCHLEAR/VENTRAL MARKERS,
                                'SULF1',  'LRP2', 'GAS1', 'PTCH1',              # SHH SIGNALING
                                 'ATOH1','CCER2','KCNH6','GRXCR2','MYO7A','LHX3','POU4F3',  # HAIR CELLS
                                'EPCAM','FBXO2', # OTIC MARKERS,
                                
                            ],
  # dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
)


# In[ ]:




