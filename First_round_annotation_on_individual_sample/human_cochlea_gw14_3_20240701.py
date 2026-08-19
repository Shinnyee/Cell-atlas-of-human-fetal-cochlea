#!/usr/bin/env python
# coding: utf-8

# In[1]:


cd F:\PROJECTS\PROJECT_HUMAN_FETAL_COCHLEAE\WORKPLACE\R\


# In[2]:


ls


# In[2]:


import scanpy as sc
sc.set_figure_params()
adata=sc.read_10x_h5("filtered_feature_bc_matrix_gw14.h5")


# In[21]:


adata


# In[22]:


adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`


# In[23]:


adata


# In[25]:


adata=sc.AnnData(adata.X,obs=adata.obs,var=adata.var)
adata.var["Gene"]=adata.var_names
adata.obs["CellID"]=adata.obs_names


# In[26]:


adata


# In[27]:


adata.obs['batch']='GW14'


# In[28]:


adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# In[32]:


sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.3, multi_panel=True)


# In[33]:


sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


# In[34]:


adata


# In[35]:


adata_filtered = adata[adata.obs.n_genes_by_counts < 5000, :]
adata_filtered = adata_filtered[adata_filtered.obs.total_counts > 500, :]
adata_filtered = adata_filtered[adata_filtered.obs.pct_counts_mt < 5, :]
adata_filtered


# In[36]:


sc.pl.scatter(adata_filtered, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata_filtered, x='total_counts', y='n_genes_by_counts')


# In[37]:


sc.pp.normalize_total(adata_filtered, target_sum=1e4)
sc.pp.log1p(adata_filtered)
sc.pp.highly_variable_genes(adata_filtered, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata_filtered)


# In[38]:


adata_filtered.raw = adata_filtered


# In[39]:


adata_filtered


# In[40]:


adata_filtered = adata_filtered[:, adata_filtered.var.highly_variable]
adata_filtered


# In[41]:


sc.pp.regress_out(adata_filtered, ['total_counts'])
sc.pp.scale(adata_filtered, max_value=10)
sc.tl.pca(adata_filtered, svd_solver='arpack',n_comps=100)
sc.pl.pca(adata_filtered, color=['MYO7A','TMC1','OTOG','OTOGL','SNAP25','MBP' ] )


# In[47]:


sc.pl.pca(adata_filtered, color=['MYO7A','TMC1','ISL1',       'OTOG','OTOGL', 'GATA3','SNAP25','MBP' ] )


# In[3]:


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


# In[4]:


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


# In[91]:


data_dir = '/PROJECTS/PROJECT_HUMAN_FETAL_COCHLEAE/WORKPLACE/R/'


# In[93]:


adata_filtered


# In[94]:


adata_filtered.write("HUMAN_cochlea_GW14.h5ad")


# In[95]:


adata = sc.read(data_dir + 'HUMAN_cochlea_GW14.h5ad')
adata


# In[98]:


sc.pp.neighbors(adata, n_neighbors=50, n_pcs=100)
sc.tl.leiden(adata,resolution=2)
sc.tl.louvain(adata)
sc.tl.umap(adata,min_dist=0.4)


# FOR THE FIRST ROUD OF SCREENING, WE USED WELL-ESTABLISHED MARKERS FROM VAB DE VALK, ER AL. 2023,CELL REPORT AND LOCHER, ET AL. 2013 NEURAL DEVELOPMENT.
# WE REASONED THAT HUAM FETAL COCHLEA AT GW14 STAGE STILL MANIFEST AS EARLY AS THE COUNTERPART IN 2023 CELL REPORT PAPER, AS FIRST APPERACE OF HAIR CELL IS IN THE BASAL TURN AS EARLY AS GW14 THAT HAVE BEEN DISCRIBED IN 2013 NEURAL DEVELOPMENT.

# In[103]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata, color=[
                           "EPCAM",   # EPITHELIUM
                        "OTOGL","OTOG","USH1C", # PAN-SC
                      'SOX2','SOX9','SOX10', # PROSENSORY DOMAIN
                        'OTX2','OC90', # ROOF CELLS, LATERAL WALL
                        'TECTA', # FLOOR MEDIAL
                        'GATA3',  # FLOOR LATERAL
                             "PRRX1",#MESENCHYMAL
                               "PRPH","TUBB3",   #NEURONAL
                            "ACAN",#CHONDROCYTES
                           "MPZ",    #GLIAL
                          "PECAM1",#ENDOTHELIAL  
                           "PTPRC",#MACROPHAGES
                         "MLANA", #MELANOCYTES
                          "MYO7A","STRC", # HAIR CELLS  
                      "MKI67","TOP2A","HMGB2", #CYCLING
                          'leiden'])


# In[126]:


small_marker_dict={
    'EPITHELIUM':["EPCAM"],
     'PAN-SC':["OTOGL","OTOG","USH1C",],
    'PROSENSORY DOMAIN':['SOX2','SOX9','SOX10',"RORB","ISL1","LGR5",],
    'ROOF CELLS, LATERAL WALL':[  'OTX2','OC90',],
    'MESENCHYMAL':["PRRX1",],
    'FLOOR LATERAL':['GATA3','FGFR3','PROX1'],
    'FLOOR MEDIAL':['TECTA','FGF20'],
'HAIR CELLS':["MYO7A","STRC","OTOF","SLC26A5"],
'NEURONAL':["PRPH","TUBB3",],
'GLIAL':["MPZ",],
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


# In[127]:


smarker_genes_in_data


# In[128]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.dendrogram(adata,
    groupby="leiden",)
sc.pl.dotplot(
    adata,
    groupby="leiden",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    save="_cochlea-cell_type_gw14_FIRST_ROUND_SCREENING.pdf"
)


# In[110]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata, color=[
                           "EPCAM",   # EPITHELIUM
                        "OTOGL","OTOG","USH1C", # PAN-SC
                      'SOX2','SOX9','SOX10', # PROSENSORY DOMAIN
                        'OTX2','OC90', # ROOF CELLS, LATERAL WALL
                        'TECTA', # FLOOR MEDIAL
                        'GATA3',  # FLOOR LATERAL
                             "PRRX1",#MESENCHYMAL
                               "PRPH","TUBB3",   #NEURONAL
                            "ACAN",#CHONDROCYTES
                           "MPZ",    #GLIAL
                          "PECAM1",#ENDOTHELIAL  
                           "PTPRC",#MACROPHAGES
                         "MLANA", #MELANOCYTES
                          "MYO7A","STRC", # HAIR CELLS  
                      "MKI67","TOP2A","HMGB2", #CYCLING
                          'leiden'],
          save="_cochlea-cell_type_gw14_FIRST_ROUND_SCREENING.pdf")


# In[111]:


sc.pl.umap(adata, color=['leiden'],legend_loc="on data")


# In[2]:


# AFTER MANUAL CELL-TYPE ANNOTATION, THERE IS STILL TWO CLUSTER THAT ARE AMBIGOUS, DIFFICULT TO DEFINED, CLUSTER 18 AND CLUSTER 31
# REMOVE CLUSTER 18 & 31 PRIOR TO DOWNSTREAM ANALYSIS
exclude_clusters = ['18','31',]
adata_rm= adata[~adata.obs['leiden'].isin(exclude_clusters), :]
adata_rm


# In[124]:


adata 


# In[125]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata, color=['leiden'],legend_loc="on data")
sc.pl.umap(adata_rm, color=['leiden'],legend_loc="on data")


# In[3]:


import scanpy as sc
adata_rm=sc.read("human_gw14_raw.h5ad")


# In[4]:


adata_rm


# In[5]:


cluster_annotation = {
    '0': 'Mesenchymal',
    '1': 'Mesenchymal',
    '2': 'Mesenchymal',
    '3': 'Mesenchymal',
    '4': 'Mesenchymal',
    '5': 'Mesenchymal',
    '6': 'Mesenchymal',
    '7': 'CoE_Medial',
    '8': 'CoE_Lateral',
    '9': 'Mesenchymal',
    '10': 'CoE_Medial',
    '11': 'Mesenchymal',
    '12': 'CoE_Lateral',
    '13': 'CoE_Roof_cells',
    '14': 'Mesenchymal',
    '15': 'CoE_Medial',
    '16': 'CCs',
    '17': 'Chondrocytes',
   
    '19': 'Mesenchymal',
    '20': 'Macrophages',
    '21': 'Mesenchymal',
    '22': 'Mesenchymal',
    '23': 'Endothelial',
    '24': 'GCs',
    '25': 'CoE_Lateral',
    '26': 'Melanocytes',
    '27': 'HCs',
    '28': 'Mesenchymal',
    '29': 'CoE_prosensory',
    '30': 'SGNs',
   
}
adata_rm.obs['cell_type'] = adata_rm.obs['leiden'].map(cluster_annotation).astype('category')


# In[6]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_rm,color=["leiden","cell_type"],legend_loc="on data",legend_fontsize =10)


# In[7]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_rm,color=["leiden","cell_type"],legend_fontsize =10)


# In[8]:


sc.pp.neighbors(adata_rm, n_neighbors=30, n_pcs=80)
sc.tl.umap(adata_rm,min_dist=0.4,)


# In[9]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_rm,color=["leiden","cell_type"],legend_fontsize =8,legend_loc="on data",)


# In[10]:


sc.set_figure_params(figsize=(6,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_rm, ['n_genes_by_counts'], groupby='cell_type',rotation=90)


# In[11]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.correlation_matrix(adata_rm, 'cell_type',
                         save="cell-type-correlation_HUMAN_GW14.pdf")


# In[12]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_rm,color='cell_type',legend_loc='on data',frameon=False, legend_fontsize=10, legend_fontoutline=False,
           title="Human_fetal_cochlea_GW14 n=12,081 nuclei",
          save="_HUMAN_GW14_CELL_TYPE_ANNO.pdf")
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_rm,color='cell_type',frameon=False, legend_fontsize=10, legend_fontoutline=False,
           title="Human_fetal_cochlea_GW14 n=12,081 nuclei",
          save="_HUMAN_GW14_CELL_TYPE_ANNO_2.pdf")


# In[13]:


adata_rm


# In[14]:


adata_rm.write("human_gw14_raw.h5ad")


# In[6]:


adata_rm=sc.read("human_gw14_raw.h5ad")


# In[7]:


adata_rm


# In[8]:


small_marker_dict={
    'EPITHELIUM':["EPCAM"],
     #'PAN-SC':["OTOGL","OTOG","USH1C",],
    'PROSENSORY DOMAIN':['SOX2','SOX9','SOX10',"RORB","ISL1","LGR5",],
    'ROOF CELLS, LATERAL WALL':[  'OTX2','OC90',],
    'MESENCHYMAL':["PRRX1",],
    'FLOOR LATERAL':['GATA3','FGFR3','PROX1'],
    'FLOOR MEDIAL':['TECTA','FGF20'],
'HAIR CELLS':["MYO7A","STRC","OTOF","SLC26A5","ATOH1"],
'NEURONAL':["PRPH","TUBB3",],
'GLIAL':["MPZ",],
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
        if marker in adata_rm.var.index:
            markers_found.append(marker)
    smarker_genes_in_data[ct] = markers_found
#del [] # remove the last marker
del_markers = list()
for ct, markers in smarker_genes_in_data.items():
    if markers==[]:
        del_markers.append(ct)
for ct in del_markers:
    del smarker_genes_in_data[ct]


# In[9]:


smarker_genes_in_data


# In[13]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.tl.dendrogram(adata_rm, groupby='cell_type')
sc.pl.dotplot(
    adata_rm,
    groupby="cell_type",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
     save="cochlea-cell_type_refined_human_gw14-order.pdf"
)


# In[16]:


small_marker_dict={
 'CYCLING' :[ "MKI67","TOP2A","HMGB2",],
   'MESENCHYMAL':["PRRX1",],
     'EPITHELIUM':["EPCAM"],
    'FLOOR MEDIAL':['TECTA','FGF20'],
     'FLOOR LATERAL':['GATA3','FGFR3','PROX1'],
    'PROSENSORY DOMAIN':['SOX2','SOX9','SOX10',"RORB","ISL1","LGR5",],
    'ENDOTHELIAL' :[ "PECAM1"],
    'MACROPHAGES' :[ "PTPRC"],
    'GLIAL':["MPZ",],
    'MELANOCYTES' :[ "MLANA"],
      'ROOF CELLS, LATERAL WALL':[  'OTX2','OC90',],
    'CHONDROCYTES' :[ "ACAN"],
    'HAIR CELLS':["MYO7A","STRC","OTOF","SLC26A5","ATOH1","INSM1","IKZF2","TBX2"],
'NEURONAL':["PRPH","TUBB3",],
     #'PAN-SC':["OTOGL","OTOG","USH1C",],
        
}
# check if the markers are in the data
smarker_genes_in_data = dict()
for ct, markers in small_marker_dict.items():
    markers_found = list()
    for marker in markers:
        if marker in adata_rm.var.index:
            markers_found.append(marker)
    smarker_genes_in_data[ct] = markers_found
#del [] # remove the last marker
del_markers = list()
for ct, markers in smarker_genes_in_data.items():
    if markers==[]:
        del_markers.append(ct)
for ct in del_markers:
    del smarker_genes_in_data[ct]


# In[17]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.tl.dendrogram(adata_rm, groupby='cell_type')
sc.pl.dotplot(
    adata_rm,
    groupby="cell_type",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
     save="cochlea-cell_type_refined_human_gw14-order.pdf"
)


# In[29]:


sc.set_figure_params(figsize=(6,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_rm, ['n_genes_by_counts'], 
             save="cochlea-cell_type_refined_human_gw14_gene_detection.pdf",
             groupby='cell_type',rotation=90)


# In[30]:


sc.set_figure_params(figsize=(6,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_rm, ['total_counts'], 
             save="cochlea-cell_type_refined_human_gw14_total_counts.pdf",
             groupby='cell_type',rotation=90)


# In[31]:


sc.set_figure_params(figsize=(6,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_rm, ['pct_counts_mt'], 
             save="cochlea-cell_type_refined_human_gw14_mt_pct.pdf",
             groupby='cell_type',rotation=90)


# In[8]:


small_marker_dict={
    'Surrounding structures':["OSR2","BMP6","SLC7A11"],
     'Osteocytes':['CCL11'],
    'Reissner membrane':["SLC26A7"],
    'Osteoblasts':["DLX5","BGLAP","PHEX", "DMP1","IBSP"],
    'Smooth Muscle Cells':["TAGLN"],
    'Neutrophils':["LY6G"],
'Monocytes':["CD68"],
'macrophages':['AIF1','CD163','PTPRC',],
'T Cells':["CD3G"],
    'Mast Cells':["PRSS34"],
 'B Cells':["CD19"],
      'NK Cells':["NKG7"],
  'Erythrocytes':["RHD"],
'Endo.Cells':['VWF','ESAM','PECAM1'],
'pericytes':['RGS5',],
 'macrophages':['AIF1','CD163'],
    'OHC':['TMC1','MYO7A','SLC26A5','CALB1'],
    'IHC':['OTOF','SLC17A8'],
    'TypeII_SGN':['PRPH','ANO2'],
    'TypeI_SGN' :['SNAP25','NEFL','PVALB','CALB2','PROX1','ESRRG','TUBB3'],
   'Schwann.Cells':['MPZ','PRX'],
'Satellite.glial.cells':['MAG','MOG','MOBP'],
    'Glial Precursor Cells':["OLIG1","OLIG2"],
     'Is':['MLANA','TYR','DCT'],
 'Fb':['OTOS','CAR3','COL9A2','COL9A3',],
  'TBC':['EMILIN2','NOTUM','RARRES1'],
 'Sp/Rt':['SLC26A4','ANXA1'],
     'Bs':['CLDN11','ATP6V0A4','TJP1'],
 'Ms':['KCNE1','ESRRB','DCLK1'],
 'Pan_SC':['OTOG','OTOGL','USH1C','GATA3'],
   'IBC_IPh_HeC':['S100A1','SLC1A3'],
      'CC_ISC_OSC':['APOE','EPYC'],
    'IDC':['OTOA'],
    'DC':['LGR5','FGFR3','PROX1','CEACAM16'],
    'PC' :['SMPX','LGR6','ENAH'],
   
}
# check if the markers are in the data
smarker_genes_in_data = dict()
for ct, markers in small_marker_dict.items():
    markers_found = list()
    for marker in markers:
        if marker in adata_gw14.var.index:
            markers_found.append(marker)
    smarker_genes_in_data[ct] = markers_found
#del [] # remove the last marker
del_markers = list()
for ct, markers in smarker_genes_in_data.items():
    if markers==[]:
        del_markers.append(ct)
for ct in del_markers:
    del smarker_genes_in_data[ct]


# In[9]:


smarker_genes_in_data


# In[10]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.dotplot(
    adata_gw14,
    groupby="cell_type",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    save="cochlea-cell_type_gw14.pdf"
)


# In[201]:


small_marker_dict={'Type I':["SYT2","EPHA4","TSC22D3","KCNA1","RYR2"],
                  'Type II':["GATA3","SMAD6","ATP2B4", "ANO2","TH","PRPH"], 
'Type IA':["CALB2","MDGA1","B3GAT1","RXRG"],
'Type IB':["SEMA3E","NTNG1", "CALB1","LRRC52"],
  'Type IC':["RUNX1","GRM8","POU4F1","LYPD1" ],   
    
}

# check if the markers are in the data
smarker_genes_in_data = dict()
for ct, markers in small_marker_dict.items():
    markers_found = list()
    for marker in markers:
        if marker in adata_rm.var.index:
            markers_found.append(marker)
    smarker_genes_in_data[ct] = markers_found
#del [] # remove the last marker
del_markers = list()
for ct, markers in smarker_genes_in_data.items():
    if markers==[]:
        del_markers.append(ct)
for ct in del_markers:
    del smarker_genes_in_data[ct]


# In[202]:


smarker_genes_in_data


# In[204]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.dotplot(
    adata_rm,
    groupby="cell_type",
    var_names=smarker_genes_in_data,
   #dendrogram=True,
    cmap='Blues',
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
     save="_SGN_human_gw14_selective_markers.pdf"
)


# In[35]:


adata_rm.obs['cell_type'].value_counts()


# In[37]:


adata_rm.write("human_gw14_raw.h5ad")


# In[16]:


adata_gw14=sc.read("human_gw14_raw.h5ad")
adata_gw14


# In[17]:


adata=sc.read_10x_h5("filtered_feature_bc_matrix_gw14.h5")
adata


# In[18]:


adata.X.max()


# In[19]:


adata_new=adata[adata_gw14.obs.index]
adata_new


# In[20]:


adata_new.X.max()


# In[21]:


adata_gw11=adata_gw14.copy()
adata_new.obsp['connectivities']=adata_gw11.obsp['connectivities']
adata_new.obsp['distances']=adata_gw11.obsp['distances']

adata_new.obsm['X_pca']=adata_gw11.obsm['X_pca']
adata_new.obsm['X_umap']=adata_gw11.obsm['X_umap']

adata_new.uns['cell_type_colors']=adata_gw11.uns['cell_type_colors']
adata_new.uns['dendrogram_cell_type']=adata_gw11.uns['dendrogram_cell_type']
adata_new.uns['dendrogram_leiden']=adata_gw11.uns['dendrogram_leiden']
adata_new.uns['hvg']=adata_gw11.uns['hvg']
adata_new.uns['leiden']=adata_gw11.uns['leiden']
adata_new.uns['leiden_colors']=adata_gw11.uns['leiden_colors']
adata_new.uns['log1p']=adata_gw11.uns['log1p']
adata_new.uns['louvain']=adata_gw11.uns['louvain']
adata_new.uns['neighbors']=adata_gw11.uns['neighbors']
adata_new.uns['pca']=adata_gw11.uns['pca']
adata_new.uns['umap']=adata_gw11.uns['umap']

adata_new.var['gene_ids']=adata_gw11.var['gene_ids']
adata_new.var['feature_types']=adata_gw11.var['feature_types']
adata_new.var['Gene']=adata_gw11.var['Gene']
adata_new.var['n_cells_by_counts']=adata_gw11.var['n_cells_by_counts']
adata_new.var['mean_counts']=adata_gw11.var['mean_counts']
adata_new.var['pct_dropout_by_counts']=adata_gw11.var['pct_dropout_by_counts']
adata_new.var['total_counts']=adata_gw11.var['total_counts']
adata_new.var['highly_variable']=adata_gw11.var['highly_variable']
adata_new.var['means']=adata_gw11.var['means']
adata_new.var['dispersions']=adata_gw11.var['dispersions']
adata_new.var['dispersions_norm']=adata_gw11.var['dispersions_norm']

adata_new.obs['batch']=adata_gw11.obs['batch']
adata_new.obs['CellID']=adata_gw11.obs['CellID']
adata_new.obs['n_genes_by_counts']=adata_gw11.obs['n_genes_by_counts']
adata_new.obs['total_counts']=adata_gw11.obs['total_counts']
adata_new.obs['total_counts_mt']=adata_gw11.obs['total_counts_mt']
adata_new.obs['pct_counts_mt']=adata_gw11.obs['pct_counts_mt']
adata_new.obs['leiden']=adata_gw11.obs['leiden']
adata_new.obs['louvain']=adata_gw11.obs['louvain']
adata_new.obs['cell_type']=adata_gw11.obs['cell_type']


# In[22]:


# Saving count data
adata_new.layers["counts"] = adata_new.X.copy()


# In[23]:


# Normalizing to median total counts
sc.pp.normalize_total(adata_new)
# Logarithmize the data
sc.pp.log1p(adata_new)


# In[24]:


adata_new.layers["logcounts"] = adata_new.X.copy()


# In[25]:


adata_new.layers["logcounts"].max()


# In[26]:


sc.pl.umap(adata_gw14, color=['cell_type', 
                                
                            ],
           size=20,
           legend_fontsize=6,ncols = 2,wspace = 0.5,
           
           
          )
sc.pl.umap(adata_new, color=['cell_type', 
                                
                            ],
           size=20,
           legend_fontsize=6,ncols = 2,wspace = 0.5,
           
           
          )


# In[27]:


adata_new.write("human_gw14_rawcounts.h5ad")


# In[5]:


adata_gw14.obs['cell_type'].value_counts()


# In[6]:


adata_gw14_coe=adata_gw14[adata_gw14.obs['cell_type'].isin(['CoE_Medial','CoE_Lateral',
                                                   'CoE_Roof_cells','CoE_prosensory',
                                                           'HCs',]
    
)]
adata_gw14_coe


# In[7]:


adata_gw14_coe.write("human_gw14_3_raw_coe.h5ad")


# In[ ]:




