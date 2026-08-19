#!/usr/bin/env python
# coding: utf-8

# In[1]:


cd F:\PROJECTS\PROJECT_HUMAN_FETAL_COCHLEAE\WORKPLACE\R\


# In[2]:


import scanpy as sc
sc.set_figure_params()
adata = sc.read_10x_mtx('./filtered_feature_bc_matrix_23W')
print(adata)


# In[3]:


adata


# In[4]:


adata.obs['batch']='gw23'


# In[5]:


adata.obs_names_make_unique()
adata.var_names_make_unique()


# In[6]:


adata=sc.AnnData(adata.X,obs=adata.obs,var=adata.var)
adata.var["Gene"]=adata.var_names
adata.obs["CellID"]=adata.obs_names


# In[7]:


adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# In[8]:


sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.3, multi_panel=True)


# In[9]:


sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


# In[10]:


adata_filtered = adata[adata.obs.n_genes_by_counts < 5000, :]
adata_filtered = adata_filtered[adata_filtered.obs.total_counts > 500, :]
adata_filtered = adata_filtered[adata_filtered.obs.pct_counts_mt < 5, :]
adata_filtered


# In[11]:


sc.pl.scatter(adata_filtered, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata_filtered, x='total_counts', y='n_genes_by_counts')


# In[12]:


sc.pp.normalize_total(adata_filtered, target_sum=1e4)
sc.pp.log1p(adata_filtered)
sc.pp.highly_variable_genes(adata_filtered, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata_filtered)


# In[13]:


adata_filtered.raw = adata_filtered


# In[14]:


adata_filtered = adata_filtered[:, adata_filtered.var.highly_variable]
adata_filtered


# In[15]:


sc.pp.regress_out(adata_filtered, ['total_counts'])
sc.pp.scale(adata_filtered, max_value=10)
sc.tl.pca(adata_filtered, svd_solver='arpack',n_comps=100)
sc.pl.pca(adata_filtered, color=['MYO7A','TMC1','OTOG','OTOGL','SNAP25','MBP' ] )


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


# In[18]:


adata=adata_filtered
adata


# In[19]:


adata=adata.raw.to_adata()
adata


# In[20]:


sc.pp.neighbors(adata, n_neighbors=50, n_pcs=100)
sc.tl.leiden(adata,resolution=2)
sc.tl.louvain(adata)
sc.tl.umap(adata,min_dist=0.4)


# In[21]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata, color=[
                              "MYO7A",'TMC1','PCP4','OTOF','STRC','SLC26A5','SLC17A8', 'CALB1',    # HAIR CELLS
                           "EPCAM",   # EPITHELIUM
                        "OTOGL","OTOG","USH1C", # PAN-SC
                      'SOX2','SOX9','SOX10', # PROSENSORY DOMAIN
                        'OTX2','OC90', # ROOF CELLS, LATERAL WALL
                        'TECTA', # FLOOR MEDIAL
                        'GATA3',  # FLOOR LATERAL
                             "PRRX1",#MESENCHYMAL
                               "PRPH","TUBB3","SNAP25",   #NEURONAL
                            "ACAN",#CHONDROCYTES
                           "MPZ",    #GLIAL
                          "PECAM1",#ENDOTHELIAL  
                           "PTPRC",#MACROPHAGES
                         "MLANA", #MELANOCYTES
                         
                      "MKI67","TOP2A","HMGB2", #CYCLING
                          'leiden'])


# In[22]:


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


# In[23]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.dendrogram(adata,
    groupby="leiden",)
sc.pl.dotplot(
    adata,
    groupby="leiden",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    save="_cochlea-cell_type_gw23_FIRST_ROUND_SCREENING_1.pdf"
)


# In[24]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata,color=["leiden",],legend_fontsize =8,legend_loc="on data",)


# In[25]:


small_marker_dict={
     'Neutrophils':["LY6G"],
    'NK Cells':["NKG7"],
    'Monocytes':["CD68"],
'Mast Cells':["PRSS34"],
    'B Cells':["CD19"],
    'T Cells':["CD3G"],
    'Erythrocytes':["RHD"],
   'Reissner membrane':["SLC26A7"],
    'Smooth Muscle Cells':["TAGLN"],
    'Lateral wall':['MLANA',],
     'Intermediate stria':['TYR','DCT'],
     'Scala vestibuli border cells':['FXYD2'],
    'Marginal stria':['KCNE1','ESRRB','DCLK1'],
    'Glial Precursor Cells':["OLIG1","OLIG2"],
    'Glial.cells':['MAG','MOG','MOBP'],
    'Schwann.Cells':['MPZ','PRX'],
     'Fibrocytes':['OTOS','CAR3','COL9A2','COL9A3'],
 'Surrounding structures':["OSR2","BMP6","SLC7A11","CHRDL1","COL1A1","ALDH1A2"],
     'Osteocytes':["DMP1","PHEX",],
      'Osteoblasts':["DLX5","BGLAP", "IBSP","RUNX2","IFITM5"],
     'Basal stria':['CLDN11','ATP6V0A4','TJP1'],
    'pericytes':['RGS5'],
    'Endothelial cells':['VWF','ESAM'],
     'Tympanic border cells':['EMILIN2','NOTUM','RARRES1'],
    'macrophages':['AIF1','CD163'],
    'Spindle/Root_cells':['SLC26A4','ANXA1'],
    'Pan_Supporting_cells':['OTOG','OTOGl','USH1C','GATA3'],
    'Inner_border-phalangeal/Hensen_cells':['S100A1','SLC1A3'],
     'Deiters_cells':['LGR5','FGFR3','PROX1','CEACAM16'],
     'Pillar_cells' :['SMPX','LGR6','ENAH'],
      'Interdental_cells':['OTOA'],
     'Claudius/Inner-Outer_sulcus_cells':['APOE','EPYC'],
     'TypeII_SGN':['PRPH','ANO2'],
     'TypeI_SGN' :['SNAP25','NEFL','PVALB','CALB2','LYPD1','RUNX1'],
    'HC':['TMC1','MYO7A','ATOH1','CALB1','POU4F3','STRC'],
    'OHC':['IKZF2','INSM1','SLC26A5','KCNQ4'],
       'IHC':['OTOF','SLC17A8','TBX2'],
   
    
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


# In[26]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.dendrogram(adata,
    groupby="leiden",)
sc.pl.dotplot(
    adata,
    groupby="leiden",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    save="_cochlea-cell_type_gw23_FIRST_ROUND_SCREENING_2.pdf"
)


# In[73]:


# 1st remove ambigous clusrer then subdivide particular clusters
exclude_clusters = ['0','11','7','8','10','12','16','21','22','25']
adata_rm= adata[~adata.obs['leiden'].isin(exclude_clusters), :]
adata_rm


# In[74]:


sc.pl.umap(adata_rm, color='leiden',legend_loc='on data')


# In[75]:


#subdivide cluster 4 for rm,bs,root cells
sc.tl.leiden(adata_rm, resolution=0.3, restrict_to=('leiden', ['4']))
sc.pl.umap(adata_rm, color='leiden_R',legend_loc='on data')


# In[76]:


adata_rm.obs['leiden']=adata_rm.obs['leiden_R']


# In[77]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.dendrogram(adata_rm,
    groupby="leiden",)
sc.pl.dotplot(
    adata_rm,
    groupby="leiden",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    save="_cochlea-cell_type_gw23_FIRST_ROUND_SCREENING_3.pdf"
)


# In[78]:


#subdivide cluster 27 for SMC AND PERICYTES
sc.tl.leiden(adata_rm, resolution=0.4, restrict_to=('leiden', ['27']))
sc.pl.umap(adata_rm, color='leiden_R',legend_loc='on data')


# In[79]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.dendrogram(adata_rm,
    groupby="leiden_R",)
sc.pl.dotplot(
    adata_rm,
    groupby="leiden_R",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    save="_cochlea-cell_type_gw23_FIRST_ROUND_SCREENING_4.pdf"
)


# In[80]:


exclude_clusters = ['27,2',]
adata_rm= adata_rm[~adata_rm.obs['leiden_R'].isin(exclude_clusters), :]
adata_rm


# In[81]:


adata_rm.obs['leiden']=adata_rm.obs['leiden_R']


# In[83]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.dendrogram(adata_rm,
    groupby="leiden",)
sc.pl.dotplot(
    adata_rm,
    groupby="leiden",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    save="_cochlea-cell_type_gw23_FIRST_ROUND_SCREENING_4.pdf"
)


# In[84]:


#subdivide cluster 2,6 for PRE-/OSTEOBLAST OSTEOCYTES
sc.tl.leiden(adata_rm, resolution=0.4, restrict_to=('leiden', ['2','6']))
sc.pl.umap(adata_rm, color='leiden_R',legend_loc='on data')


# In[86]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.dendrogram(adata_rm,
    groupby="leiden_R",)
sc.pl.dotplot(
    adata_rm,
    groupby="leiden_R",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    save="_cochlea-cell_type_gw23_FIRST_ROUND_SCREENING_5.pdf"
)


# In[87]:


adata_rm.obs['leiden']=adata_rm.obs['leiden_R']


# In[88]:


sc.pl.umap(adata_rm, color='leiden',legend_loc='on data')


# In[93]:


cluster_annotation = {
    '31': 'HCs',
    '33': 'SGN',
    '29': 'GC',
    '13': 'Endothelial',
    '19': 'Macrophages',
    '3': 'Intermediate_stria',
    '30': 'Intermediate_stria',
    '9': 'Chondrocytes',
    '23': 'CCs',
    '17': 'Marginal stria',
    '18': 'Fibrocytes',
     '4,0': 'Reissner_membrane',   
   '4,1': 'Basal_stria',
 '27,0': 'pericytes',
    '27,1': 'Smooth_Muscle_Cells',
    '24': 'Tympanic_border_cells',   
    '1': 'Tympanic_border_cells',
   '5': 'Tympanic_border_cells',
    '32': 'Spindle/Root_cells',
 '2-6,0': 'pre-Osteoblasts',
  '2-6,1': 'Osteoblasts', 
'2-6,2': 'Osteocytes',
    '15': 'Inner_border-phalangeal/Hensen_cells',
 '28': 'Deiters_cells', 
  '20': 'Pillar_cells',
     '26': 'Interdental_cells',
    
    '14': 'Claudius/Inner-Outer_sulcus_cells',
     
    
   
}
adata_rm.obs['cell_type'] = adata_rm.obs['leiden'].map(cluster_annotation).astype('category')


# In[94]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_rm,color=["leiden","cell_type"],legend_loc="on data",legend_fontsize =5)


# In[95]:


sc.pp.neighbors(adata_rm, n_neighbors=50, n_pcs=100)
sc.tl.umap(adata_rm,min_dist=0.4,)


# In[96]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_rm,color=["leiden","cell_type"],legend_fontsize =5,legend_loc="on data",)


# In[97]:


exclude_clusters = ['1',]
adata_rm= adata_rm[~adata_rm.obs['leiden'].isin(exclude_clusters), :]
adata_rm


# In[130]:


sc.pp.neighbors(adata_rm, n_neighbors=60, n_pcs=80)
sc.tl.umap(adata_rm,min_dist=0.5,)


# In[131]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_rm,color=["leiden","cell_type"],legend_fontsize =5,legend_loc="on data",)


# In[133]:


cluster_annotation = {
    '31': 'HCs',
    '33': 'SGNs',
    '29': 'GCs',
    '13': 'Endothelial',
    '19': 'Macrophages',
    '3': 'Intermediate_stria',
    '30': 'Intermediate_stria',
    '9': 'Chondrocytes',
    '23': 'CCs',
    '17': 'Marginal_stria',
    '18': 'Fibrocytes',
     '4,0': 'Reissner_membrane',   
   '4,1': 'Basal_stria',
 '27,0': 'pericytes',
    '27,1': 'Smooth_Muscle_Cells',
    '24': 'Tympanic_border_cells',   
   '5': 'Tympanic_border_cells',
    '32': 'Spindle/Root_cells',
 '2-6,0': 'pre-Osteoblasts',
  '2-6,1': 'Osteoblasts', 
'2-6,2': 'Osteocytes',
    '15': 'Inner_border-phalangeal/Hensen_cells',
 '28': 'Deiters_cells', 
  '20': 'Pillar_cells',
     '26': 'Interdental_cells',
    
    '14': 'Claudius/Inner-Outer_sulcus_cells',
     
    
   
}
adata_rm.obs['cell_type'] = adata_rm.obs['leiden'].map(cluster_annotation).astype('category')


# In[134]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_rm,color=["leiden","cell_type"],legend_fontsize =5,legend_loc="on data",)


# In[135]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.correlation_matrix(adata_rm, 'cell_type',
                         save="cell-type-correlation_HUMAN_GW23.pdf")


# In[136]:


adata_rm


# In[137]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_rm,color='cell_type',legend_loc='on data',frameon=False, legend_fontsize=5, legend_fontoutline=False,
           title="Human_fetal_cochlea_GW23 n=3,597 nuclei",
          save="_HUMAN_GW23_CELL_TYPE_ANNO.pdf")
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_rm,color='cell_type',frameon=False, legend_fontsize=8, legend_fontoutline=False,
           title="Human_fetal_cochlea_GW23 n=3,597 nuclei",
          save="_HUMAN_GW23_CELL_TYPE_ANNO_2.pdf")


# In[138]:


adata_rm.write("human_gw23_raw.h5ad")


# In[139]:


small_marker_dict={
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
     'Endothelial cells':['VWF','ESAM'],
     'Smooth Muscle Cells':["TAGLN"],
     'pericytes':['RGS5'],
    'macrophages':['AIF1','CD163'],
     'Cycling cells' :[ "MKI67","TOP2A","HMGB2",],
     'Erythrocytes':["RHD"],
     'Glial Precursor Cells':["OLIG1","OLIG2"],
    'Glial.cells':['MAG','MOG','MOBP'],
    'Schwann.Cells':['MPZ','PRX'],
 'Intermediate stria':['TYR','DCT'],
     'Melanocytes' :[ "MLANA"],
   
     'Neutrophils':["LY6G"],
    'NK Cells':["NKG7"],
    'Monocytes':["CD68"],
'Mast Cells':["PRSS34"],
    'B Cells':["CD19"],
    'T Cells':["CD3G"],
     'Surrounding structures':["OSR2","BMP6","SLC7A11","CHRDL1","COL1A1","ALDH1A2"],
     'TypeII_SGN':['PRPH','ANO2'],
     'TypeI_SGN' :['SNAP25','NEFL','PVALB','CALB2','LYPD1','RUNX1'],
    'OHC':['TMC1','MYO7A','SLC26A5','CALB1','POU4F3'],
       'IHC':['OTOF','SLC17A8','STRC'],
       
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


# In[140]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.tl.dendrogram(adata_rm, groupby='cell_type')
sc.pl.dotplot(
    adata_rm,
    groupby="cell_type",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
     save="cochlea-cell_type_refined_human_gw23-order.pdf"
)


# In[143]:


small_marker_dict={

 'Osteoblasts':["BGLAP", "IBSP","IFITM5"],
      'Osteocytes':["DMP1","PHEX",],   
'Chondrocytes' :["ACAN"],
 'pre-Osteoblasts':["DLX5","RUNX2",],
'Deiters_cells':['LGR5','FGFR3','PROX1','CEACAM16'],
     'Pillar_cells' :['SMPX','LGR6','ENAH'],
 'Claudius/Inner-Outer_sulcus_cells':['APOE','EPYC'],
    'Inner_border-phalangeal/Hensen_cells':['S100A1','SLC1A3'],
    'Interdental_cells':['OTOA'],
     'Pan_Supporting_cells':['OTOG','OTOGl','USH1C','GATA3'],
 'Marginal stria':['KCNE1','ESRRB','DCLK1'],  
     'Basal stria':['CLDN11','ATP6V0A4','TJP1'],
 'Reissner membrane':["SLC26A7"],
     'HC':['TMC1','MYO7A','CALB1','POU4F3','STRC'],
    'OHC':['SLC26A5','OCM'],
       'IHC':['OTOF','SLC17A8'], 
     'TypeII_SGN':['PRPH','ANO2'],
     'TypeI_SGN' :['SNAP25','NEFL','PVALB','CALB2','LYPD1','RUNX1'],

 'Glial Precursor Cells':["OLIG1","OLIG2"],
    'Glial.cells':['MAG','MOG','MOBP'],
    'Schwann.Cells':['MPZ','PRX'],
    'Melanocytes' :[ "MLANA"],
 'Intermediate stria':['TYR','DCT'],
     'macrophages':['AIF1','CD163'],
'Spindle/Root_cells':['SLC26A4','ANXA1'],
     'Cycling cells' :[ "MKI67","TOP2A","HMGB2",],
 'Fibrocytes':['OTOS','CAR3','COL9A2','COL9A3'],
   
    'Tympanic border cells':['EMILIN2','NOTUM','RARRES1'],
 'Endothelial cells':['VWF','ESAM'],
 'Smooth Muscle Cells':["TAGLN"],
     'pericytes':['RGS5'],
    
   
       
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


# In[144]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.tl.dendrogram(adata_rm, groupby='cell_type')
sc.pl.dotplot(
    adata_rm,
    groupby="cell_type",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
     save="cochlea-cell_type_refined_human_gw23-order.pdf"
)


# In[145]:


sc.set_figure_params(figsize=(6,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_rm, ['n_genes_by_counts'], 
             save="cochlea-cell_type_refined_human_gw23_gene_detection.pdf",
             groupby='cell_type',rotation=90)


# In[146]:


sc.set_figure_params(figsize=(6,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_rm, ['total_counts'], 
             save="cochlea-cell_type_refined_human_gw23_total_counts.pdf",
             groupby='cell_type',rotation=90)


# In[147]:


sc.set_figure_params(figsize=(6,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_rm, ['pct_counts_mt'], 
             save="cochlea-cell_type_refined_human_gw23_mt_pct.pdf",
             groupby='cell_type',rotation=90)


# In[148]:


adata_rm


# In[149]:


adata_rm.obs['cell_type'].value_counts()


# In[150]:


adata_rm.write("human_gw23_raw.h5ad")


# In[4]:


adata=sc.read("human_gw23_raw.h5ad")


# In[5]:


adata


# In[8]:


sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


# In[17]:


sc.set_figure_params(figsize=(6,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata, ['CHRNA10',], 
             save="cochlea-cell_type_refined_human_gw23_CHRNA10.pdf",
             groupby='cell_type',rotation=90)


# In[3]:


adata_rm=sc.read("human_gw23_raw.h5ad")
adata_rm


# In[4]:


adata = sc.read_10x_mtx('./filtered_feature_bc_matrix_23W')
print(adata)


# In[5]:


adata.X.max()


# In[6]:


adata_new=adata[adata_rm.obs.index]
adata_new


# In[7]:


adata_new.X.max()


# In[8]:


adata_gw11=adata_rm.copy()
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


# In[9]:


# Saving count data
adata_new.layers["counts"] = adata_new.X.copy()


# In[10]:


# Normalizing to median total counts
sc.pp.normalize_total(adata_new)
# Logarithmize the data
sc.pp.log1p(adata_new)


# In[11]:


adata_new.layers["logcounts"] = adata_new.X.copy()


# In[12]:


adata_new.layers["logcounts"].max()


# In[13]:


sc.pl.umap(adata_rm, color=['cell_type', 
                                
                            ],
           size=20,
           legend_fontsize=6,ncols = 2,wspace = 0.5,
           
           
          )
sc.pl.umap(adata_new, color=['cell_type', 
                                
                            ],
           size=20,
           legend_fontsize=6,ncols = 2,wspace = 0.5,
           
           
          )


# In[15]:


adata_new.var_names_make_unique()


# In[16]:


adata_new.write("human_gw23_rawcounts.h5ad")


# In[5]:


adata_rm.obs['cell_type'].value_counts()


# In[6]:


adata_gw23_coe=adata_rm[adata_rm.obs['cell_type'].isin(['Claudius/Inner-Outer_sulcus_cells','Interdental_cells',
                                                        'Deiters_cells','Pillar_cells','Inner_border-phalangeal/Hensen_cells',
                                                   'HCs',
                                                           ]
    
)]
adata_gw23_coe


# In[7]:


adata_gw23_coe.write("human_gw23_raw_coe.h5ad")


# In[ ]:




