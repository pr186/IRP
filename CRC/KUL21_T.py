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
adata = anndata.read("/home/pr186/Documents/IRP/FINAL/CRC_RAW/KUL21_T/adata.h5ad")
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
model_prefix = 'KUL21_T'

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


# UMAP

sc.tl.umap(adata)

## EXPRESSION OF TARGET GENES ON UMAP
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
        genes_TARGET_PR.append(genes_TARGET[gene])
    except KeyError:
        print(genes_TARGET[gene] + " not present")
        
sc.pl.umap(adata, color = genes_TARGET_PR, save = '_TARGET_ALL.png')

for gene in range(len(genes_TARGET2)):
    try:
        sc.pl.umap(adata, color = genes_TARGET2[gene], save = '_' + genes_TARGET2[gene] + '_TARGET2.png')
        genes_TARGET2_PR.append(genes_TARGET2[gene])
    except KeyError:
        print(genes_TARGET2[gene] + " not present")
        
sc.pl.umap(adata, color = genes_TARGET2_PR, save = '_TARGET2_ALL.png')


for gene in range(len(genes_PROLINE)):
    try:
        sc.pl.umap(adata, color = genes_PROLINE[gene], save = '_' + genes_PROLINE[gene] + '_PROLINE.png')
        genes_PROLINE_PR.append(genes_PROLINE[gene])
    except KeyError:
        print(genes_PROLINE[gene] + " not present")
        
sc.pl.umap(adata, color = genes_PROLINE_PR, save = '_PROLINE_ALL.png')
        
for gene in range(len(genes_GLYCOLYSIS)):
    try:
        sc.pl.umap(adata, color = genes_GLYCOLYSIS[gene], save = '_' + genes_GLYCOLYSIS[gene] + '_GLYCOLYSIS.png')
        genes_GLYCOLYSIS_PR.append(genes_GLYCOLYSIS[gene])
    except KeyError:
        print(genes_GLYCOLYSIS[gene] + " not present")

sc.pl.umap(adata, color = genes_GLYCOLYSIS_PR, save = '_GLYCOLYSIS_ALL.png')




# CLUSTERING

sc.pl.umap(adata, save = 'EMPTY_UMAP.png')


sc.tl.louvain(adata,resolution=0.5, random_state=42)
sc.pl.umap(adata, color=['louvain'], save = '_louvain_8.png')

sc.tl.leiden(adata, resolution = 0.5)
sc.tl.leiden(adata2, resolution = 0.5)

sc.pl.umap(adata ,color='leiden', save = '_leiden_10.png')

####### AUTHORS ANN

adata.obs["barcode_id"] = "KUL21-T_" + adata.obs.index.values

metadata = pd.read_csv("/home/pr186/Documents/IRP/CRC_FINAL/annotation.txt", header = 0, sep = "\t")
metadata.index = metadata.Index
adata.obs["celltype"] = adata.obs.barcode_id.map(metadata["Cell_type"])
adata.obs["cellsubtype"] = adata.obs.barcode_id.map(metadata["Cell_subtype"])


sc.pl.umap(adata, color=['celltype'], save = '_annotated_AUTHORS1_21.png')

sc.pl.umap(adata, color=['cellsubtype'], save = '_annotated_AUTHORS2_21.png')

############# TRY CELLO


cello.scanpy_cello(
  adata2,
  clust_key='leiden',
  rsrc_loc=cello_resource_loc,
  out_prefix=model_prefix,
  log_dir=os.getcwd()
)


# File to store results

results_file = 'KUL21_T_CELLO.h5ad'
adata2.write(results_file, compression='gzip')

adata2.obs["barcode_id"] = "KUL21-T_" + adata.obs.index.values
adata2.obs.index = adata2.obs["barcode_id"]


adata.obs["Most specific cell type"] = adata.obs.barcode_id.map(adata2.obs["Most specific cell type"])

adata.obs.head()

sc.pl.umap(adata, color=['Most specific cell type'], save = '_annotated_CELLO_21.png')





# Try genes of interest

adata.var.index = adata.var["gene_id"]


ax = sc.pl.dotplot(adata, genes_TARGET_PR, gene_symbols = "gene_name", groupby='leiden', use_raw = False, save = '_leiden_TARGET_25.png')
ax = sc.pl.dotplot(adata, genes_TARGET2_PR,  gene_symbols = "gene_name", groupby='leiden', use_raw = False, save = '_leiden_TARGET2_25.png')
ax = sc.pl.dotplot(adata, genes_PROLINE_PR, gene_symbols = "gene_name",  groupby='leiden', use_raw = False, save = '_leiden_PROLINE_26.png')
ax = sc.pl.dotplot(adata, genes_GLYCOLYSIS_PR,  gene_symbols = "gene_name", groupby='leiden', use_raw = False, save = '_leiden_GLYCOLYSIS_27.png')

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

scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, color = 'louvain', basis='umap', frameon=True, legend_loc = 'right', save = '_louvain_umap_31.png')

scv.pl.velocity_embedding_stream(adata, color = 'leiden', basis='umap', frameon=True, legend_loc = 'right', save = '_leiden_umap_33.png')

scv.pl.velocity_embedding_stream(adata, color = 'celltype', basis='umap', frameon=True, legend_loc = 'right', save = '_origANN5.png')
scv.pl.velocity_embedding_stream(adata, color = 'cellsubtype', basis='umap', frameon=True, legend_loc = 'right', save = '_origANN6.png')
scv.pl.velocity_embedding_stream(adata, color = 'Most specific cell type', basis='umap', frameon=True, legend_loc = 'right', save = '_CELLO.png')


adata.var.index = adata.var["gene_name"]


for gene in range(len(genes_TARGET_PR)):
    scv.pl.velocity(adata, genes_TARGET_PR[gene], basis = "umap", figsize = (8, 6), save = "_" + genes_TARGET_PR[gene] + '_UMAP_TARGET.png')

for gene in range(len(genes_TARGET2_PR)):
    scv.pl.velocity(adata, genes_TARGET2_PR[gene], basis = "umap", figsize = (8, 6), save = "_" + genes_TARGET2_PR[gene] + '_UMAP_TARGET2.png')

for gene in range(len(genes_PROLINE_PR)):
    scv.pl.velocity(adata, genes_PROLINE_PR[gene], basis = "umap", figsize = (8, 6), save = "_" + genes_PROLINE_PR[gene] + '_UMAP_PROLINE.png')

for gene in range(len(genes_GLYCOLYSIS_PR)):
    scv.pl.velocity(adata, genes_GLYCOLYSIS_PR[gene], basis = "umap", figsize = (8, 6), save = "_" + genes_GLYCOLYSIS_PR[gene] + '_UMAP_GLYCOLYSIS.png')

scv.pl.velocity_embedding(adata, basis = "umap", frameon=True, arrow_length=3, arrow_size=2, dpi=300, save = '_34.png')



scv.tl.rank_velocity_genes(adata, groupby='louvain', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()
df.to_csv("velocity_genes_KUL21.csv")

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































