#!/usr/bin/env python
# coding: utf-8

# In[1]:


cd F:\PROJECTS\PROJECT_HUMAN_FETAL_COCHLEAE\WORKPLACE\R\


# In[2]:


import scanpy as sc
sc.set_figure_params()
adata=sc.read_10x_h5("coch-331_filtered_feature_bc_matrix.h5")


# In[3]:


adata


# In[4]:


adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`


# In[5]:


adata=sc.AnnData(adata.X,obs=adata.obs,var=adata.var)
adata.var["Gene"]=adata.var_names
adata.obs["CellID"]=adata.obs_names


# In[6]:


adata.obs['batch']='GW22'


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


# In[16]:


sc.pl.pca(adata_filtered, color=['MYO7A','TMC1','ISL1','OTOG','OTOGL', 'GATA3','SNAP25','MBP' ] )


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


# In[19]:


data_dir = '/PROJECTS/PROJECT_HUMAN_FETAL_COCHLEAE/WORKPLACE/R/'


# In[20]:


adata_filtered


# In[21]:


adata_filtered.write("HUMAN_cochlea_GW22_2.h5ad")


# In[23]:


adata = sc.read(data_dir + 'HUMAN_cochlea_GW22_2.h5ad')
adata


# In[24]:


adata=adata.raw.to_adata()
adata


# In[25]:


sc.pp.neighbors(adata, n_neighbors=50, n_pcs=100)
sc.tl.leiden(adata,resolution=2)
sc.tl.louvain(adata)
sc.tl.umap(adata,min_dist=0.4)


# In[26]:


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


# In[41]:


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


# In[42]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.dendrogram(adata,
    groupby="leiden",)
sc.pl.dotplot(
    adata,
    groupby="leiden",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    save="_cochlea-cell_type_gw22_2_FIRST_ROUND_SCREENING_1.pdf"
)


# In[82]:


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
    'OHC':['TMC1','MYO7A','SLC26A5','CALB1','POU4F3'],
       'IHC':['OTOF','SLC17A8','STRC'],
   
    
}


# In[83]:


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


# In[84]:


smarker_genes_in_data


# In[86]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.dendrogram(adata,
    groupby="leiden",)
sc.pl.dotplot(
    adata,
    groupby="leiden",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    save="_cochlea-cell_type_gw22_2_FIRST_ROUND_SCREENING_2.pdf"
)


# In[47]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata,color=["leiden",],legend_fontsize =8,legend_loc="on data",)


# In[78]:


adata1=adata.copy()
adata1


# In[79]:


sc.tl.leiden(adata1, resolution=0.3, restrict_to=('leiden', ['17']))
sc.pl.umap(adata1, color='leiden_R',legend_loc='on data')


# In[87]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.dendrogram(adata1,
    groupby="leiden_R",)
sc.pl.dotplot(
    adata1,
    groupby="leiden_R",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    save="_cochlea-cell_type_gw22_2_FIRST_ROUND_SCREENING_3.pdf"
)


# In[88]:


sc.tl.leiden(adata1, resolution=0.3, restrict_to=('leiden', ['21']))
sc.pl.umap(adata1, color='leiden_R',legend_loc='on data')


# In[89]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.dendrogram(adata1,
    groupby="leiden_R",)
sc.pl.dotplot(
    adata1,
    groupby="leiden_R",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    save="_cochlea-cell_type_gw22_2_FIRST_ROUND_SCREENING_4.pdf"
)


# In[92]:


sc.tl.leiden(adata1, resolution=0.8, restrict_to=('leiden', ['21']))
sc.pl.umap(adata1, color='leiden_R',legend_loc='on data')


# In[93]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.dendrogram(adata1,
    groupby="leiden_R",)
sc.pl.dotplot(
    adata1,
    groupby="leiden_R",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    save="_cochlea-cell_type_gw22_2_FIRST_ROUND_SCREENING_5.pdf"
)


# In[96]:


exclude_clusters = ['4','9','22','24','27','21,0','21,1','7']
adata_rm= adata1[~adata1.obs['leiden_R'].isin(exclude_clusters), :]
adata_rm


# In[97]:


sc.pl.umap(adata_rm, color='leiden_R',legend_loc='on data')


# In[98]:


adata_rm.obs['leiden'] = adata_rm.obs['leiden_R']


# In[99]:


sc.pl.umap(adata_rm, color='leiden_R',legend_loc='on data')


# In[103]:


sc.tl.leiden(adata_rm, resolution=0.3, restrict_to=('leiden', ['17']))
sc.pl.umap(adata_rm, color='leiden_R',legend_loc='on data')


# In[104]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.tl.dendrogram(adata_rm,
    groupby="leiden_R",)
sc.pl.dotplot(
    adata_rm,
    groupby="leiden_R",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    save="_cochlea-cell_type_gw22_2_FIRST_ROUND_SCREENING_6.pdf"
)


# In[105]:


adata_rm.obs['leiden'] = adata_rm.obs['leiden_R']


# In[106]:


sc.pl.umap(adata_rm, color='leiden',legend_loc='on data')


# In[107]:


adata_rm.obs['leiden'].value_counts()


# In[157]:


cluster_annotation = {
    '0': 'GCs',
    '1': 'Osteoblasts',
    '2': 'Osteocytes',
    '3': 'Intermediate_stria',
    '5': 'Tympanic_border_cells',
    '6': 'Fibrocytes',
    '8': 'Tympanic_border_cells',
    '10': 'Chondrocytes',
    '11': 'Macrophages',
    '12': 'pre-Osteoblasts',
    '13': 'Endothelial',
    '14': 'Chondrocytes',
    '15': 'Chondrocytes',
    '16': 'Osteocytes',
      '18': 'Mesenchymal',   
    '19': 'Claudius/Inner-Outer_sulcus_cells',
    '20': 'Claudius/Inner-Outer_sulcus_cells',
 '17,0': 'Basal_stria/Spindle/Root_cells',
     '17,1': 'Reissner_membrane',
    '23': 'Chondrocytes',
    '25': 'Chondrocytes',
    '26': 'GCs',
    '28': 'Claudius/Inner-Outer_sulcus_cells',
    '29': 'Marginal_stria',
    '30': 'CCs',
     '31': 'Inner_border-phalangeal/Hensen_cells',
     '32': 'Interdental_cells',
     '33': 'Endothelial',
     '34': 'Melanocytes',
     '35': 'Erythrocytes',
     '36': 'Pillar_cells',
     '37': 'Deiters_cells',
     '21,2': 'Smooth_Muscle_Cells',
     '21,3': 'pericytes',
   
}
adata_rm.obs['cell_type'] = adata_rm.obs['leiden'].map(cluster_annotation).astype('category')


# In[158]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_rm,color=["leiden","cell_type"],legend_loc="on data",legend_fontsize =5)


# In[159]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_rm,color=["leiden","cell_type"],legend_fontsize =5)


# In[160]:


sc.pp.neighbors(adata_rm, n_neighbors=80, n_pcs=100)
sc.tl.umap(adata_rm,min_dist=0.6,)


# In[161]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.pl.umap(adata_rm,color=["leiden","cell_type"],legend_fontsize =5,legend_loc="on data",)


# In[162]:


sc.set_figure_params(figsize=(6,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_rm, ['n_genes_by_counts'], groupby='cell_type',rotation=90)


# In[163]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.correlation_matrix(adata_rm, 'cell_type',
                         save="cell-type-correlation_HUMAN_GW22_2.pdf")


# In[164]:


adata_rm


# In[165]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_rm,color='cell_type',legend_loc='on data',frameon=False, legend_fontsize=5, legend_fontoutline=False,
           title="Human_fetal_cochlea_GW22_2 n=7,861 nuclei",
          save="_HUMAN_GW22_2_CELL_TYPE_ANNO.pdf")
sc.set_figure_params(figsize=(4,4),frameon=False,dpi=150,dpi_save=600)
sc.pl.umap(adata_rm,color='cell_type',frameon=False, legend_fontsize=8, legend_fontoutline=False,
           title="Human_fetal_cochlea_GW22_2 n=7,861 nuclei",
          save="_HUMAN_GW22_2_CELL_TYPE_ANNO_2.pdf")


# In[166]:


adata_rm.write("human_gw22_2_331_raw.h5ad")


# In[167]:


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


# In[168]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.tl.dendrogram(adata_rm, groupby='cell_type')
sc.pl.dotplot(
    adata_rm,
    groupby="cell_type",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
     save="cochlea-cell_type_refined_human_gw22_2-order.pdf"
)


# In[169]:


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


# In[170]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.tl.dendrogram(adata_rm, groupby='cell_type')
sc.pl.dotplot(
    adata_rm,
    groupby="cell_type",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
     save="cochlea-cell_type_refined_human_gw22_2-order.pdf"
)


# In[171]:


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


# In[172]:


sc.set_figure_params(figsize=(4,4),frameon=False,dpi=300,dpi_save=600)
sc.tl.dendrogram(adata_rm, groupby='cell_type')
sc.pl.dotplot(
    adata_rm,
    groupby="cell_type",
    var_names=smarker_genes_in_data,
   dendrogram=True,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
     save="cochlea-cell_type_refined_human_gw22_2-order.pdf"
)


# In[173]:


sc.set_figure_params(figsize=(6,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_rm, ['n_genes_by_counts'], 
             save="cochlea-cell_type_refined_human_gw22_2_gene_detection.pdf",
             groupby='cell_type',rotation=90)


# In[174]:


sc.set_figure_params(figsize=(6,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_rm, ['total_counts'], 
             save="cochlea-cell_type_refined_human_gw22_2_total_counts.pdf",
             groupby='cell_type',rotation=90)


# In[175]:


sc.set_figure_params(figsize=(6,3),frameon=False,dpi=150,dpi_save=600)
sc.pl.violin(adata_rm, ['pct_counts_mt'], 
             save="cochlea-cell_type_refined_human_gw22_2_mt_pct.pdf",
             groupby='cell_type',rotation=90)


# In[176]:


adata_rm


# In[177]:


adata_rm.write("human_gw22_2_331_raw.h5ad")


# In[3]:


adata_rm=sc.read("human_gw22_2_331_raw.h5ad")
adata_rm


# In[4]:


adata=sc.read_10x_h5("coch-331_filtered_feature_bc_matrix.h5")
adata


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


adata_new.write("human_gw22_2_331_rawcounts.h5ad")


# In[5]:


adata_rm.obs['cell_type'].value_counts()


# In[6]:


adata_gw22_coe=adata_rm[adata_rm.obs['cell_type'].isin(['Claudius/Inner-Outer_sulcus_cells','Interdental_cells',
                                                        'Deiters_cells','Pillar_cells','Inner_border-phalangeal/Hensen_cells',
                                                   
                                                           ]
    
)]
adata_gw22_coe


# In[7]:


adata_gw22_coe.write("human_gw22_2_331_raw_coe.h5ad")


# In[ ]:




