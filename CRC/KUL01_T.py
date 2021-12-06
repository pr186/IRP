### TRYING CELLLO HERE

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from anndata import AnnData
import scvelo as scv
import cello
import os

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=120)

# load the unfiltered matrix
results_file = 'KUL01_T.h5ad'  # the file that will store the analysis results
adata = anndata.read("/home/pr186/Documents/IRP/FINAL/CRC_RAW/KUL01_T/adata.h5ad")
adata

adata.var["gene_id"] = adata.var.index.values

adata.var_names_make_unique()


adata.var.head()

#adata.var.index = adata.var["gene_name"]


# Investigating highest expressed genes
#sc.pl.highest_expr_genes(adata, n_top=20, save = '_1.png')


### BASIC FILTERING

# Removes genes with 0 umi counts
adata = adata[:, np.asarray(adata.X.sum(axis=0)).reshape(-1) > 0]

#adata.var_names_make_unique()

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# FILTER OUT HIGH GENES COUNT

# add the total counts per cell as observations-annotation to adata
#adata.obs['n_counts'] = adata.X.sum(axis=1).A1

#sc.pl.scatter(adata, x='n_counts', y='n_genes', save = '_beforeFiltering_2.png')

# Create a mask to filter out cells with more than 6500 genes
adata = adata[adata.obs.n_genes < 6500, :]
adata


############################################## FOR CELLO!!!!!!!!!!!

cellxgenes = adata.to_df()
adata2 = AnnData(cellxgenes)

adata2.obs.head()

cello_resource_loc = "/home/pr186/Documents/IRP/CRC_2/KUL01_T_FINAL/CELLO"
model_prefix = 'KUL01_T'

adata2.var_names_make_unique()


###################################################################

adata.var.index = adata.var["gene_name"]


# Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell,
# so that counts become comparable among cells.

# normalize counts in each cell to be equal
sc.pp.normalize_total(adata, target_sum=10**4)
sc.pp.normalize_total(adata2, target_sum=10**4)

# Logarithmize the data / Replace raw counts with their logarithm
sc.pp.log1p(adata)
sc.pp.log1p(adata2)


# Lets now look again at the highest expressed genes after filtering, normalization, and log
sc.pl.highest_expr_genes(adata, n_top=20, save = '_afterNormalization_5.png')

# Set the .raw attribute of AnnData object to the normalized and logarithmized raw gene expression 
# for later use in differential testing and visualizations of gene expression. 
# This simply freezes the state of the AnnData object.

adata.raw = adata

# Identify highly-variable genes.

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.highly_variable_genes(adata2, min_mean=0.0125, max_mean=3, min_disp=0.5)


sc.pl.highly_variable_genes(adata, save = '_DISPERSION_PLOT_6.png')


# Regress out effects of total counts per cell and the percentage of mitochondrial 
# genes expressed. Scale the data to unit variance.
# sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
# !!!!!!!! We do not regress out as per https://github.com/theislab/scanpy/issues/526

# SCALING THE DATA 

sc.pp.scale(adata, max_value=10)
#sc.pp.scale(adata2, max_value=10)


# PCA - PRINCIPAL COMPONENT ANALYSIS

sc.tl.pca(adata, svd_solver='arpack')
sc.tl.pca(adata2, svd_solver='arpack')


# Let us inspect the contribution of single PCs to the total variance in the data.

sc.pl.pca_variance_ratio(adata, log=True, save = '_7.png')

adata

# Compute the neighborhood graph !!! check values

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.pp.neighbors(adata2, n_neighbors=10, n_pcs=40)


# UMAP & TSNE

sc.tl.umap(adata)
sc.tl.tsne(adata)

## EXPRESSION OF TARGET GENES ON UMAP & TSNE
genes_TARGET = ['DACH1', 'DACH2', 'PAX6', 'SIX1', 'SIX2', 'EYA1', 'EYA2']
genes_TARGET_PR = []
genes_TARGET2 = ["SIX1", "SIX2", "SIX3", "SIX4", "SIX5", "SIX6", "PAX4", "PAX6", 
              "DACH1", "DACH2", "EYA1", "EYA2", "EYA3", "EYA4", "EYA5", "EYA6"]
genes_TARGET2_PR = []
genes_PROLINE = ['PYCR1', 'PYCR2', 'PYCR3', 'PRODH', 'PRODH2', 'OAT', 'ALDH18A1', 'ALDH4A1']
genes_PROLINE_PR = []
genes_GLYCOLYSIS = ['HK1', 'HK2', 'HK3', 'GPI', 'PFKFB1', 'PFKFB2', 'PFKFB3', 
                    'PFKFB4', 'ALDOA', 'ALDOB', 'ALDOC', 'TRI', 'GAPDH', 'PGK1', 
                    'PGK2', 'PGAM1', 'ENO1', 'ENO2', 'ENO3', 'PKM', 
                    'PFK1', 'PFK2', 'LDHA', 'LDHB', 'GLUT1', 'SLC2A1']
genes_GLYCOLYSIS_PR = []

for gene in range(len(genes_TARGET)):
    try:
        sc.pl.umap(adata, color = genes_TARGET[gene], save = '_' + genes_TARGET[gene] + '_TARGET.png')
        sc.pl.tsne(adata, color = genes_TARGET[gene], save = '_' + genes_TARGET[gene] + '_TARGET.png')
        genes_TARGET_PR.append(genes_TARGET[gene])
    except KeyError:
        print(genes_TARGET[gene] + " not present")
        
sc.pl.umap(adata, color = genes_TARGET_PR, save = '_TARGET_ALL.png')
sc.pl.tsne(adata, color = genes_TARGET_PR, save = '_TARGET_ALL.png')

for gene in range(len(genes_TARGET2)):
    try:
        sc.pl.umap(adata, color = genes_TARGET2[gene], save = '_' + genes_TARGET2[gene] + '_TARGET2.png')
        sc.pl.tsne(adata, color = genes_TARGET2[gene], save = '_' + genes_TARGET2[gene] + '_TARGET2.png')
        genes_TARGET2_PR.append(genes_TARGET2[gene])
    except KeyError:
        print(genes_TARGET2[gene] + " not present")
        
sc.pl.umap(adata, color = genes_TARGET2_PR, save = '_TARGET2_ALL.png')
sc.pl.tsne(adata, color = genes_TARGET2_PR, save = '_TARGET2_ALL.png')


for gene in range(len(genes_PROLINE)):
    try:
        sc.pl.umap(adata, color = genes_PROLINE[gene], save = '_' + genes_PROLINE[gene] + '_PROLINE.png')
        sc.pl.tsne(adata, color = genes_PROLINE[gene], save = '_' + genes_PROLINE[gene] + '_PROLINE.png')
        genes_PROLINE_PR.append(genes_PROLINE[gene])
    except KeyError:
        print(genes_PROLINE[gene] + " not present")
        
sc.pl.umap(adata, color = genes_PROLINE_PR, save = '_PROLINE_ALL.png')
sc.pl.tsne(adata, color = genes_PROLINE_PR, save = '_PROLINE_ALL.png')
        
for gene in range(len(genes_GLYCOLYSIS)):
    try:
        sc.pl.umap(adata, color = genes_GLYCOLYSIS[gene], save = '_' + genes_GLYCOLYSIS[gene] + '_GLYCOLYSIS.png')
        sc.pl.tsne(adata, color = genes_GLYCOLYSIS[gene], save = '_' + genes_GLYCOLYSIS[gene] + '_GLYCOLYSIS.png')
        genes_GLYCOLYSIS_PR.append(genes_GLYCOLYSIS[gene])
    except KeyError:
        print(genes_GLYCOLYSIS[gene] + " not present")

sc.pl.umap(adata, color = genes_GLYCOLYSIS_PR, save = '_GLYCOLYSIS_ALL.png')
sc.pl.tsne(adata, color = genes_GLYCOLYSIS_PR, save = '_GLYCOLYSIS_ALL.png')




# CLUSTERING

sc.pl.umap(adata, save = 'EMPTY_UMAP.png')
sc.pl.tsne(adata, save = 'EMPTY_TSNE.png')


sc.tl.louvain(adata,resolution=0.5, random_state=42)
sc.pl.umap(adata, color=['louvain'], save = '_louvain_8.png')
sc.pl.tsne(adata, color=['louvain'], save = '_louvain_9.png')

sc.tl.leiden(adata, resolution = 0.5)
sc.tl.leiden(adata2, resolution = 0.5)

sc.pl.umap(adata ,color='leiden', save = '_leiden_10.png')
sc.pl.tsne(adata, color='leiden', save = '_leiden_11.png')

# FIND MARKERS GENES - LOUVAIN

sc.tl.rank_genes_groups(adata, 'louvain', method='t-test', corr_method="bonferroni")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save = '_tTEST_louvain_12.png')

sc.settings.verbosity = 2  # reduce the verbosity

sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon', corr_method="bonferroni")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save = '_wilcoxon_louvain_13.png')

# SHOW & SAVE TOP 10 genes in table

pd.set_option('display.max_columns', None)
df1 = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(25)
df1
df1.to_csv("louvain_top25.csv")

# FIND MARKERS GENES - LEIDEN

sc.tl.rank_genes_groups(adata, 'leiden', method='t-test', corr_method="bonferroni")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save = '_tTEST_leiden_15.png')

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', corr_method="bonferroni")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save = '_wilcoxon_leiden_16.png')

# SHOW & SAVE TOP 10 genes in table

pd.set_option('display.max_columns', None)
df2 = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(25)
df2
df2.to_csv("leiden_top25.csv")

### GENERATE A CSV FILE FOR ANNOTATION OF CLUSTERS WITH SCSA

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
dat = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals']})
dat.to_csv("SCSA_KUL01_T.csv")


new_cluster_names_SCSA = [
                     "Cancer Cells (Cluster 0)",
                     "Cancer Stem Cells (Cluster 1)",
                     "B Cells",
                     "T Cells",
                     "Cancer Stem Cells (Cluster 4)",
                     "Cancer Stem Cells (Cluster 5)",
                     "Masenchymal Cells",
                     "Cancer Stem Cells (Cluster 7)",
                     "Cancer Stem Cells (Cluster 8)",
                     "Neutrophil",
                     "Cancer Stem Cells (Cluster 10)",
                     "Cancer Stem Cells (Cluster 11)",
                     "Cancer Stem Cells (Cluster 12)",
                     "Cancer Stem Cells (Cluster 13)"
]

new_cluster_names_leiden = [
                     "Malignant Cells (Cluster 0)",
                     "Malignant Cells (Cluster 1)",
                     "Malignant Cells (Cluster 2)",
                     "CD4 & CD8 T Cells",
                     "Dendritic Cells, Macrophages, Monocytes (Cluster 4)",
                     "Fibroblast (Cluster 5)",
                     "Fibroblast, Endothelial Cells",
                     "Endothelial Cells (Cluster 7)",
                     "Fibroblast (Cluster 8)",
                     "Dendritic Cells, Macrophages, Monocytes (Cluster 9)",
                     "Fibroblast (Cluster 10)",
                     "NK Cells",
                     "Malignant Cells (Cluster 12)",
                     "Endothelial Cells (Cluster 13)",
]

adata.obs.head()
adata.var.head()

####### AUTHORS ANN

adata.obs["barcode_id"] = "KUL01-T_" + adata.obs.index.values

metadata = pd.read_csv("/home/pr186/Documents/IRP/CRC_FINAL/annotation.txt", header = 0, sep = "\t")
metadata.index = metadata.Index
adata.obs["celltype"] = adata.obs.barcode_id.map(metadata["Cell_type"])
adata.obs["cellsubtype"] = adata.obs.barcode_id.map(metadata["Cell_subtype"])


# Relabel the clusters
adata.obs['leiden_2'] = adata.obs['leiden'].cat.rename_categories(new_cluster_names_leiden)
adata.obs['SCSA'] = adata.obs['leiden'].cat.rename_categories(new_cluster_names_SCSA)


sc.pl.tsne(adata, color=['leiden_2'], save = '_annotated_leiden_20.png')
sc.pl.umap(adata, color=['leiden_2'], save = '_annotated_leiden_21.png')

sc.pl.tsne(adata, color=['SCSA'], save = '_annotated_SCSA_20.png')
sc.pl.umap(adata, color=['SCSA'], save = '_annotated_SCSA_21.png')

sc.pl.tsne(adata, color=['celltype'], save = '_annotated_AUTHORS1_20.png')
sc.pl.umap(adata, color=['celltype'], save = '_annotated_AUTHORS1_21.png')

sc.pl.tsne(adata, color=['cellsubtype'], save = '_annotated_AUTHORS2_20.png')
sc.pl.umap(adata, color=['cellsubtype'], save = '_annotated_AUTHORS2_21.png')

############# TRY CELLO


#cello.scanpy_cello(
#  adata2,
#  clust_key='leiden',
#  rsrc_loc=cello_resource_loc,
#  out_prefix=model_prefix,
#  log_dir=os.getcwd()
#)

adata.obs.head()
adata.var.head()

cello.scanpy_cello(
  adata2,
  clust_key='leiden',
  rsrc_loc=cello_resource_loc,
  model_file=f'/home/pr186/Documents/IRP/FINAL/CRC_RAW/KUL01_T/KUL01_T.model.dill'
)

results_file = 'KUL01_T_CELLO.h5ad'
adata2.write(results_file, compression='gzip')

adata2.obs.head()
adata2.obs["barcode_id"] = "KUL01-T_" + adata.obs.index.values
adata2.obs.index = adata2.obs["barcode_id"]


adata.obs["Most specific cell type"] = adata.obs.barcode_id.map(adata2.obs["Most specific cell type"])

adata.obs.head()

sc.pl.tsne(adata, color=['Most specific cell type'], save = '_annotated_CELLO_20.png')
sc.pl.umap(adata, color=['Most specific cell type'], save = '_annotated_CELLO_21.png')


# Try genes of interest

adata.var.index = adata.var["gene_id"]


ax = sc.pl.dotplot(adata, genes_TARGET_PR, gene_symbols = "gene_name", groupby='leiden', use_raw = False, save = '_leiden_TARGET_25.png')
ax = sc.pl.dotplot(adata, genes_TARGET2_PR,  gene_symbols = "gene_name", groupby='leiden', use_raw = False, save = '_leiden_TARGET2_25.png')
ax = sc.pl.dotplot(adata, genes_PROLINE_PR, gene_symbols = "gene_name",  groupby='leiden', use_raw = False, save = '_leiden_PROLINE_26.png')
ax = sc.pl.dotplot(adata, genes_GLYCOLYSIS_PR,  gene_symbols = "gene_name", groupby='leiden', use_raw = False, save = '_leiden_GLYCOLYSIS_27.png')

ax = sc.pl.dotplot(adata, genes_TARGET_PR, gene_symbols = "gene_name", groupby='leiden_2', use_raw = False, save = '_leiden2_TARGET_25.png')
ax = sc.pl.dotplot(adata, genes_TARGET2_PR,  gene_symbols = "gene_name", groupby='leiden_2', use_raw = False, save = '_leiden2_TARGET2_25.png')
ax = sc.pl.dotplot(adata, genes_PROLINE_PR, gene_symbols = "gene_name",  groupby='leiden_2', use_raw = False, save = '_leiden2_PROLINE_26.png')
ax = sc.pl.dotplot(adata, genes_GLYCOLYSIS_PR,  gene_symbols = "gene_name", groupby='leiden_2', use_raw = False, save = '_leiden2_GLYCOLYSIS_27.png')

ax = sc.pl.dotplot(adata, genes_TARGET_PR, gene_symbols = "gene_name", groupby='celltype', use_raw = False, save = '_AUTHORS1_TARGET_25.png')
ax = sc.pl.dotplot(adata, genes_TARGET2_PR,  gene_symbols = "gene_name", groupby='celltype', use_raw = False, save = '_AUTHORS12_TARGET2_25.png')
ax = sc.pl.dotplot(adata, genes_PROLINE_PR, gene_symbols = "gene_name",  groupby='celltype', use_raw = False, save = '_AUTHORS1_PROLINE_26.png')
ax = sc.pl.dotplot(adata, genes_GLYCOLYSIS_PR,  gene_symbols = "gene_name", groupby='celltype', use_raw = False, save = '_AUTHORS1_GLYCOLYSIS_27.png')

ax = sc.pl.dotplot(adata, genes_TARGET_PR, gene_symbols = "gene_name", groupby='cellsubtype', use_raw = False, save = '_AUTHORS2_TARGET_25.png')
ax = sc.pl.dotplot(adata, genes_TARGET2_PR,  gene_symbols = "gene_name", groupby='cellsubtype', use_raw = False, save = '_AUTHORS2_TARGET2_25.png')
ax = sc.pl.dotplot(adata, genes_PROLINE_PR, gene_symbols = "gene_name",  groupby='cellsubtype', use_raw = False, save = '_AUTHORS2_PROLINE_26.png')
ax = sc.pl.dotplot(adata, genes_GLYCOLYSIS_PR,  gene_symbols = "gene_name", groupby='cellsubtype', use_raw = False, save = '_AUTHORS2_GLYCOLYSIS_27.png')

ax = sc.pl.dotplot(adata, genes_TARGET_PR, gene_symbols = "gene_name", groupby='Most specific cell type', use_raw = False, save = '_CELLO_TARGET_25.png')
ax = sc.pl.dotplot(adata, genes_TARGET2_PR,  gene_symbols = "gene_name", groupby='Most specific cell type', use_raw = False, save = '_CELLO_TARGET2_25.png')
ax = sc.pl.dotplot(adata, genes_PROLINE_PR, gene_symbols = "gene_name",  groupby='Most specific cell type', use_raw = False, save = '_CELLO_PROLINE_26.png')
ax = sc.pl.dotplot(adata, genes_GLYCOLYSIS_PR,  gene_symbols = "gene_name", groupby='Most specific cell type', use_raw = False, save = '_CELLO_GLYCOLYSIS_27.png')


### VELOCITY
scv.pl.proportions(adata, groupby = 'louvain', figsize = (20, 4), save = 'louvain_28.png')
scv.pl.proportions(adata, groupby = 'leiden', figsize = (20, 4), save = 'leiden_29.png')

scv.pp.filter_and_normalize(adata, min_shared_counts = 20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs = 30, n_neighbors = 30)

scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)


scv.pl.velocity_embedding_stream(adata, color = 'louvain', basis='tsne', frameon=True, legend_loc = 'right', save = '_louvain_tsne_30.png')
scv.pl.velocity_embedding_stream(adata, color = 'louvain', basis='umap', frameon=True, legend_loc = 'right', save = '_louvain_umap_31.png')

scv.pl.velocity_embedding_stream(adata, color = 'leiden', basis='tsne', frameon=True, legend_loc = 'right', save = '_leiden_tsne_32.png')
scv.pl.velocity_embedding_stream(adata, color = 'leiden', basis='umap', frameon=True, legend_loc = 'right', save = '_leiden_umap_33.png')

scv.pl.velocity_embedding_stream(adata, color = 'celltype', basis='umap', frameon=True, legend_loc = 'right', save = '_origANN5.png')
scv.pl.velocity_embedding_stream(adata, color = 'cellsubtype', basis='umap', frameon=True, legend_loc = 'right', save = '_origANN6.png')
scv.pl.velocity_embedding_stream(adata, color = 'Most specific cell type', basis='umap', frameon=True, legend_loc = 'right', save = '_CELLO.png')


adata.var.index = adata.var["gene_name"]


for gene in range(len(genes_TARGET_PR)):
    try:
        scv.pl.velocity(adata, genes_TARGET_PR[gene], basis = "tsne", figsize = (8, 6), save = "_" + genes_TARGET_PR[gene] + '_TSNE_TARGET.png')
        scv.pl.velocity(adata, genes_TARGET_PR[gene], basis = "umap", figsize = (8, 6), save = "_" + genes_TARGET_PR[gene] + '_UMAP_TARGET.png')
    except KeyError:
        print(genes_TARGET_PR[gene] + " not present")


for gene in range(len(genes_TARGET2_PR)):
    try:
        scv.pl.velocity(adata, genes_TARGET2_PR[gene], basis = "tsne", figsize = (8, 6), save = "_" + genes_TARGET2_PR[gene] + '_TSNE_TARGET2.png')
        scv.pl.velocity(adata, genes_TARGET2_PR[gene], basis = "umap", figsize = (8, 6), save = "_" + genes_TARGET2_PR[gene] + '_UMAP_TARGET2.png')
    except KeyError:
        print(genes_TARGET2_PR[gene] + " not present")


for gene in range(len(genes_PROLINE_PR)):
    try:
        scv.pl.velocity(adata, genes_PROLINE_PR[gene], basis = "tsne", figsize = (8, 6), save = "_" + genes_PROLINE_PR[gene] + '_TSNE_PROLINE.png')
        scv.pl.velocity(adata, genes_PROLINE_PR[gene], basis = "umap", figsize = (8, 6), save = "_" + genes_PROLINE_PR[gene] + '_UMAP_PROLINE.png')
    except KeyError:
        print(genes_PROLINE_PR[gene] + " not present")

for gene in range(len(genes_GLYCOLYSIS_PR)):
    try:
        scv.pl.velocity(adata, genes_GLYCOLYSIS_PR[gene], basis = "tsne", figsize = (8, 6), save = "_" + genes_GLYCOLYSIS_PR[gene] + '_TSNE_GLYCOLYSIS.png')
        scv.pl.velocity(adata, genes_GLYCOLYSIS_PR[gene], basis = "umap", figsize = (8, 6), save = "_" + genes_GLYCOLYSIS_PR[gene] + '_UMAP_GLYCOLYSIS.png')
    except KeyError:
        print(genes_GLYCOLYSIS_PR[gene] + " not present")
        
        
        
scv.pl.velocity_embedding(adata, basis = "umap", frameon=True, arrow_length=3, arrow_size=2, dpi=300, save = '_34.png')
scv.pl.velocity_embedding(adata, basis = "tsne", frameon=True, arrow_length=3, arrow_size=2, dpi=300, save = '_35.png')



scv.tl.rank_velocity_genes(adata, groupby='louvain', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()
df.to_csv("velocity_genes_KUL01.csv")

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], frameon=True, save = '_speedANDcoherence_36.png')

scv.tl.score_genes_cell_cycle(adata, use_raw = False)
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95], save = '_37.png', frameon=True) #, add_outline=True)

scv.pl.velocity_graph(adata, threshold=.1, legend_loc = 'right', save = '_38.png', frameon=True, add_outline=True)


x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell=70)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax, save = '_39.png')

scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', frameon=True, cmap='gnuplot', save = 'pseudotime_40.png')































