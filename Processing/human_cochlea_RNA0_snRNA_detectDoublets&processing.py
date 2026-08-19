#!/usr/bin/env python
# coding: utf-8

# In[1]:


#import cupy as cp
import scrublet as scr
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys
import scipy


def MovePlots(plotpattern, subplotdir):
    os.system('mkdir -p '+str(sc.settings.figdir)+'/'+subplotdir)
    os.system('mv '+str(sc.settings.figdir)+'/*'+plotpattern+'** '+str(sc.settings.figdir)+'/'+subplotdir)


sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.figdir = './figures/preprocessing/'
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80)  # low dpi (dots per inch) yields small inline figures

sys.executable


# # Detect doublet

# In[2]:


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


# Scrumblet
# (Courtesy of K Polansky)
# 
# Two-step doublet score processing, mirroring the approach from Popescu et al. https://www.nature.com/articles/s41586-019-1652-y which was closely based on Pijuan-Sala et al. https://www.nature.com/articles/s41586-019-0933-9
# 
# The first step starts with some sort of doublet score, e.g. Scrublet, and ends up with a per-cell p-value (with significant values marking doublets). For each sample individually:
# 
# run Scrublet to obtain each cell's score
# overcluster the manifold - run a basic Scanpy pipeline up to clustering, then additionally cluster each cluster separately
# compute per-cluster Scrublet scores as the median of the observed values, and use those going forward
# identify p-values:
# compute normal distribution parameters: centered at the median of the scores, with a MAD-derived standard deviation
# the score distribution is zero-truncated, so as per the paper I only use above-median values to compute the MAD
# K deviates from the paper a bit, at least the exact wording captured within it, and multiply the MAD by 1.4826 to obtain a literature-derived normal distribution standard deviation estimate
# FDR-correct the p-values via Benjamini-Hochberg
# write out all this doublet info into CSVs for later use
# NOTE: The second step is performed later, in a multi-sample space

# In[3]:


cd F:\PROJECTS\PROJECT_HUMAN_FETAL_COCHLEAE\WORKPLACE\R\sc_pipeline\


# In[4]:


data_dir = '/PROJECTS/PROJECT_HUMAN_FETAL_COCHLEAE/WORKPLACE/R/sc_pipeline/'
meta = pd.read_csv(data_dir+'metadata.csv',index_col=0)


# In[5]:


meta


# In[6]:


meta['sample'] = meta['sample_id'].astype('str')
plotmeta = list(meta.columns)
plotmeta.append('sample')
print('Number of sample: ', meta.index.size)
meta


# In[7]:


samples = meta.index.to_list()
samples


# In[8]:


meta['sample']


# In[23]:


#there's loads of clustering going on, so set verbosity low unless you enjoy walls of text
sc.settings.verbosity = 0  # verbosity: errors (0), warnings (1), info (2), hints (3)

scorenames = ['scrublet_score','scrublet_cluster_score','zscore','bh_pval','bonf_pval']
if not os.path.exists('scrublet-scores'):
    os.makedirs('scrublet-scores')
    #loop over the subfolders of the rawdata folder

samples = meta.index.to_list()

for sample in reversed(list(samples)):
    print(sample)
    #import data
    adata_sample = sc.read(data_dir+sample+'/human_gw_raw.h5ad',cache=True)
    adata_sample.var_names_make_unique()
    #rename cells to SAMPLE_BARCODE
    adata_sample.obs_names = [sample+'_'+i for i in adata_sample.obs_names]
    #do some early filtering to retain meaningful cells for doublet inspection
    sc.pp.filter_cells(adata_sample, min_genes=200)
    sc.pp.filter_genes(adata_sample, min_cells=3)
    #convert to lower to be species agnostic: human mito start with MT-, mouse with mt-
    mito_genes = [name for name in adata_sample.var_names if name.lower().startswith('MT-')]
    # for each cell compute fraction of counts in mito genes vs. all genes
    # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    adata_sample.obs['percent_mito'] = np.sum(
        adata_sample[:, mito_genes].X, axis=1).A1 / np.sum(adata_sample.X, axis=1).A1
    adata_sample = adata_sample[adata_sample.obs['percent_mito'] < 0.5, :]

    #set up and run Scrublet, seeding for replicability
    np.random.seed(0)
    scrub = scr.Scrublet(adata_sample.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    adata_sample.obs['scrublet_score'] = doublet_scores

    #overcluster prep. run turbo basic scanpy pipeline
    sc.pp.normalize_per_cell(adata_sample, counts_per_cell_after=1e4)
    sc.pp.log1p(adata_sample)
    sc.pp.highly_variable_genes(adata_sample, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_sample = adata_sample[:, adata_sample.var['highly_variable']]
    sc.pp.scale(adata_sample, max_value=10)
    sc.tl.pca(adata_sample, svd_solver='arpack')
    sc.pp.neighbors(adata_sample)
    #overclustering proper - do basic clustering first, then cluster each cluster
    sc.tl.leiden(adata_sample)
    adata_sample.obs['leiden'] = [str(i) for i in adata_sample.obs['leiden']]
    for clus in np.unique(adata_sample.obs['leiden']):
        adata_sub = adata_sample[adata_sample.obs['leiden']==clus].copy()
        sc.tl.leiden(adata_sub)
        adata_sub.obs['leiden'] = [clus+','+i for i in adata_sub.obs['leiden']]
        adata_sample.obs.loc[adata_sub.obs_names,'leiden'] = adata_sub.obs['leiden']

    #compute the cluster scores - the median of Scrublet scores per overclustered cluster
    for clus in np.unique(adata_sample.obs['leiden']):
        adata_sample.obs.loc[adata_sample.obs['leiden']==clus, 'scrublet_cluster_score'] = \
            np.median(adata_sample.obs.loc[adata_sample.obs['leiden']==clus, 'scrublet_score'])
    #now compute doublet p-values. figure out the median and mad (from above-median values) for the distribution
    med = np.median(adata_sample.obs['scrublet_cluster_score'])
    mask = adata_sample.obs['scrublet_cluster_score']>med
    mad = np.median(adata_sample.obs['scrublet_cluster_score'][mask]-med)
    #let's do a one-sided test. the Bertie write-up does not address this but it makes sense
    zscores = (adata_sample.obs['scrublet_cluster_score'].values - med) / (1.4826 * mad)
    adata_sample.obs['zscore'] = zscores
    pvals = 1-scipy.stats.norm.cdf(zscores)
    adata_sample.obs['bh_pval'] = bh(pvals)
    adata_sample.obs['bonf_pval'] = bonf(pvals)

    #create results data frame for single sample and copy stuff over from the adata object
    scrublet_sample = pd.DataFrame(0, index=adata_sample.obs_names, columns=scorenames)
    for score in scorenames:
        scrublet_sample[score] = adata_sample.obs[score]
    #write out complete sample scores
    scrublet_sample.to_csv('scrublet-scores/'+sample+'.csv')


# ### Preprocessing
# 
# 
# Load 10x
# Filter: 1) cells (< 10 genes); 2) genes (< 3 cells)
# 
# Quantify: 1) % mitochondrial genes; 2) total counts

# In[24]:


adata_sample


# In[25]:


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


# In[26]:


# confirm N samples
print(len(holder))
# merge datasets
adata = holder[0].concatenate(holder[1:],join='outer',index_unique=None)
# copy of this matrix in Compressed Sparse Row format
adata.X = adata.X.tocsr()
adata


# In[13]:


adata.write('all_rawcounts.h5ad')


# QC pplots
# Plot distributions of the values n_genes, n_counts and percent_mito

# In[14]:


sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True, save='.pdf', show=True)
sc.pl.scatter(adata, x='n_counts', y='percent_mito', save='_numi_vs_mito.pdf', show=True)
sc.pl.scatter(adata, x='n_counts', y='n_genes', save='_numi_vs_ngenes.pdf', show=True)


# In[15]:


sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], groupby='sample', rotation = 90, jitter=0.4, save='_batch_collection.pdf', show=False)


# In[16]:


print('Total number of cells: {:d}'.format(adata.n_obs))
print('Total number of genes: {:d}'.format(adata.n_vars))
pd.DataFrame(adata.obs).to_csv(str(sc.settings.figdir)+'/metadata_prefilters.csv')


# In[18]:


sc.pl.violin(adata[[ i == 'female' for i in adata.obs.gender]], ['SRY', 'RPS4Y1','DDX3Y'], groupby='sample', rotation = 90,save="_human_cochlea_male3.pdf",) # FEMALE
sc.pl.violin(adata[[ i == 'male' for i in adata.obs.gender]], ['SRY', 'RPS4Y1','DDX3Y'], groupby='sample', rotation = 90,save="_human_cochlea_male_2.pdf",) # MALE
sc.pl.violin(adata, ['SRY', 'RPS4Y1','DDX3Y'], groupby='gender',save="_human_cochlea_male.pdf",) # MALE


# In[19]:


adata


# Filter cells with few genes 
# Check number of genes per cell distribution and filter cells accordingly

# In[28]:


plt.hist(adata.obs['n_genes'], bins = 100,range=(0,7500))
plt.axvline(350, linestyle = '--', color = 'red')


# In[29]:


plt.hist(adata.obs['n_counts'], bins = 100,range=(0,10000))
plt.axvline(1200, linestyle = '--', color = 'red')


# In[22]:


# confirm N samples
print(len(holder))
# merge datasets
adata = holder[0].concatenate(holder[1:],join='inner',index_unique=None)
# copy of this matrix in Compressed Sparse Row format
adata.X = adata.X.tocsr()
adata


# In[30]:


adata_qc = adata.copy()
# Filter raw cells according to identified QC thresholds:
print('Total number of cells: {:d}'.format(adata_qc.n_obs))

np.quantile(adata.obs['n_counts'], [.99])
sc.pp.filter_cells(adata_qc, min_counts = 1200)
print('Number of cells after min count filter: {:d}'.format(adata_qc.n_obs))

sc.pp.filter_cells(adata_qc, max_counts = 10000)
print('Number of cells after max count filter: {:d}'.format(adata_qc.n_obs))

np.quantile(adata.obs['n_genes'], [.99])
sc.pp.filter_cells(adata_qc, min_genes = 350)
print('Number of cells after min genes filter: {:d}'.format(adata_qc.n_obs))

sc.pp.filter_cells(adata_qc, max_genes = 7500)
print('Number of cells after max genes filter: {:d}'.format(adata_qc.n_obs))

#adata_qc = adata_qc[adata_qc.obs['percent_mito'] < 0.1]
#print('Number of cells after mito filter: {:d}'.format(adata_qc.n_obs))


# In[31]:


print('Total number of cells: {:d}'.format(adata_qc.n_obs))
print('Total number of genes: {:d}'.format(adata_qc.n_vars))
adata_qc.obs['sample_id'].values.describe()


# Filter cells with large % mitochondrial genes

# In[32]:


sc.pl.violin(adata_qc, ['percent_mito'], groupby='sample', rotation = 90) 


# In[33]:


plt.hist(adata_qc.obs['percent_mito'], bins = 100, cumulative=True,range=(0,0.1))
plt.axvline(0.02, linestyle = '--', color = 'red')
plt.axvline(0.05, linestyle = '--', color = 'darkred')
plt.axhline(adata_qc.n_obs*0.99, linestyle = '-', color = 'green')


# In[34]:


# >5%
adata_qc = adata_qc[adata_qc.obs['percent_mito'] < 0.05]
print('Number of cells after mito filter: {:d}'.format(adata_qc.n_obs))


# Remove cells with low counts high mito combo

# In[35]:


x = [ int(i) < 1000 for i in adata_qc.obs['n_counts']]
adata_qc.obs['low_ncounts'] = [ str(i) for i in x ]

y = [ i > 0.1 for i in adata_qc.obs['percent_mito']]
adata_qc.obs['high_mito'] = [ str(i) for i in y ]

mask = [all(tup) for tup in zip(x,y)]
adata_qc.obs['low_ncounts_high_mito'] = [ str(i) for i in mask ]

adata_qc = adata_qc[[ 'False' in i for i in adata_qc.obs['low_ncounts_high_mito']   ]] 


# In[36]:


print('Total number of cells: {:d}'.format(adata_qc.n_obs))
print('Total number of genes: {:d}'.format(adata_qc.n_vars))
adata_qc.obs['sample'].values.describe()
pd.DataFrame(adata_qc.obs).to_csv(str(sc.settings.figdir)+'/metadata_filtered.csv')


# In[37]:


adata_qc


# ## Identify cells behaving like cc genes
# #Per genes analysis: identify genes behaving like known cell cycle genes

# In[253]:


bdata = adata_qc.copy()
# Normalize total counts per cell
sc.pp.normalize_per_cell(bdata, counts_per_cell_after=1e4)
# Logarithmize the data matrix
sc.pp.log1p(bdata)


# In[254]:


# Extract highly variable genes
sc.pp.highly_variable_genes(bdata)
highly_variable_genes = bdata.var["highly_variable"]
bdata = bdata[:, highly_variable_genes]


# In[255]:


# Traspose matrix for a GENE-centered analysis
bdata = bdata.copy().T


# In[256]:


bdata.X.shape


# In[257]:


# Scale data to unit variance and zero mean
sc.pp.scale(bdata, max_value=10)

# Scatter plot in PCA coordinates
sc.tl.pca(bdata)
bdata.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
# Plot the variance ratio
sc.pl.pca_variance_ratio(bdata, log=True, save='_ccg_identification.pdf')


# In[258]:


num_pcs = 10
# Compute a neighborhood graph of observations
sc.pp.neighbors(bdata, n_pcs=num_pcs)
# Embed the neighborhood graph using UMAP
sc.tl.umap(bdata)
# Cluster GENES into subgroups using louvain: resolution < 1 to find less clusters
sc.tl.leiden(bdata, resolution=1)


# In[259]:


# Locate ccs cluster
bdata.obs['known_cyclers'] = [i in ['CDK1','MKI67','CCNB2','PCNA'] for i in bdata.obs_names]
bdata.obs['known_cyclers'] = [ str(i) for i in  bdata.obs['known_cyclers']]
sc.pl.umap(bdata, color=['known_cyclers', 'leiden'], color_map='OrRd',save='_ccg_identification.pdf')
print(bdata.obs.loc[[i in ['CDK1','MKI67','CCNB2','PCNA'] for i in bdata.obs_names],'leiden'])


# In[260]:


ccgs_cl = bdata.obs.loc['MKI67',['leiden']][0]
print("Cell cycle genes cluster is "+ccgs_cl)


# # Flag CYCLING GENES

# In[261]:


# Add unstructured dict-like annotation for ccgs
adata_qc.uns['ccgs'] = list(bdata.obs[bdata.obs['leiden']==ccgs_cl].index)


# In[262]:


# Move plots
MovePlots('ccg_identification','ccg_identification')


# # Score cell cycle

# In[263]:


s_genes = [x.strip() for x in open('S_genes.tsv')]
g2m_genes = [x.strip() for x in open('G2M_genes.tsv')]


# In[264]:


s_genes = [x for x in s_genes if x in adata_qc.var_names]
g2m_genes = [x for x in g2m_genes if x in adata_qc.var_names]


# In[265]:


sc.tl.score_genes_cell_cycle(adata_qc, s_genes=s_genes, g2m_genes=g2m_genes)


# # Identify HVGs Flavor='seurat_v3' expects count data.
# 
# ## For more informations, see "Feature selection for individual datasets" at https://www.sciencedirect.com/science/article/pii/S0092867419305598

# In[266]:


sc.pp.highly_variable_genes(adata_qc, flavor='seurat_v3', n_top_genes=2000, subset=False)


# # Load scrublet

# In[267]:


scorenames = ['scrublet_score','scrublet_cluster_score','zscore','bh_pval','bonf_pval']

scrdf = []
for sample in meta.index:
    scrdf.append(pd.read_csv('scrublet-scores/'+sample+'.csv', header=0, index_col=0))
scrdf = pd.concat(scrdf)
scrdf.index = [i.replace('-1', '') for i in scrdf.index]

idx = [ i in adata_qc.obs_names for i in scrdf.index ]
scrdf = scrdf[idx]
for score in scorenames:
    adata_qc.obs[score] = scrdf[score]
adata_qc.obs['is_doublet'] = adata_qc.obs['bonf_pval'] < 0.01


# In[268]:


adata_qc


# In[269]:


adata_qc.uns['ccgs']


# In[270]:


np.mean(adata_qc.obs['is_doublet'])


# In[271]:


np.mean(adata_qc.obs['scrublet_score']>0.3)


# In[272]:


adata_qc.raw = adata_qc.copy()


# In[274]:


adata_qc.var


# In[275]:


adata_qc.obs_names_make_unique()
adata_qc.var_names_make_unique()


# In[276]:


adata_qc


# In[277]:


adata_qc.write('all_rawcounts.h5ad')


# # we do not remove cc genes as one cluster we named as cycling cell-type which highly expresses cell-cycling markers TOP2A,HMGB2 et al.

# Normalize per cell and log transform After removing unwanted cells and genes from the dataset, the next step is to normalize the data. By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.

# In[278]:


sc.pp.normalize_total(adata_qc, target_sum=1e4)
sc.pp.log1p(adata_qc)


# # PCA
# 
# ## 1. Filter HVGs in bdata and do PCA with them

# In[279]:


bdata = adata_qc[:, adata_qc.var['highly_variable']]
print('Total number of cells: {:d}'.format(bdata.n_obs))
print('Total number of genes: {:d}'.format(bdata.n_vars))


# In[280]:


sc.pp.scale(bdata, max_value=10)
sc.tl.pca(bdata, svd_solver='arpack', n_comps=50)


# ## 2. Transfer PCA to the main adata

# In[281]:


adata_qc.var['highly_variable'].fillna(value=False, inplace=True) # fill NaNs with False so that subsetting to HVGs is possible
adata_qc.obsm['X_pca'] = bdata.obsm['X_pca'].copy()
adata_qc.uns['pca'] = bdata.uns['pca'].copy()
adata_qc.varm['PCs'] = np.zeros(shape=(adata_qc.n_vars, 50))
adata_qc.varm['PCs'][adata_qc.var['highly_variable']] = bdata.varm['PCs']
sc.pl.pca_variance_ratio(adata_qc, log=True, save='.pdf')


# In[282]:


n_pcs = 25
sc.pp.neighbors(adata_qc, n_pcs = n_pcs)
sc.tl.umap(adata_qc)


# In[283]:


adata_qc.write('all_processed.h5ad')


# In[284]:


meta


# In[285]:


adata_qc


# In[287]:


sc.pl.umap(adata_qc, color=[ 'gender', 'sample','S_score', 'G2M_score',], save='_predoublet.pdf', ncols = 2,wspace=0.5)


# In[289]:


sc.pl.umap(adata_qc, color=[ 'cell_type',], save='_predoublet2.pdf', ncols = 2,wspace=0.5)


# In[292]:


sc.pl.umap(adata_qc, color=['n_genes', 'percent_mito'], save='_predoublet_stats.pdf', ncols = 2, color_map='OrRd', use_raw=False)


# In[294]:


sc.pl.umap(adata_qc, color=['scrublet_score', 'scrublet_cluster_score'], color_map='OrRd', save='_predoublet_stats2.pdf',)


# In[295]:


adata_qc.obs['cell_type'].value_counts()


# In[296]:


adata_hc=adata_qc[adata_qc.obs['cell_type']=='HCs']
adata_hc


# In[298]:


adata_hc.obs['sample'].value_counts()


# In[ ]:




