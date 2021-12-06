
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from anndata import AnnData
from sklearn.decomposition import TruncatedSVD
import matplotlib
import matplotlib.pyplot as plt
import scvelo as scv
import cello
import os
import scanorama

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=120)


# load the unfiltered matrix
adata_T = anndata.read("/home/pr186/Documents/IRP/FINAL/CRC_RAW/KUL28_T/adata.h5ad")
adata_T

adata_T.var.head()

adata_T.var["gene_id"] = adata_T.var.index.values

#adata_T.var.index = adata_T.var["gene_name"]


adata_T.var_names_make_unique()

#adata_T.var.index.is_unique

adata_T.var.index

adata_T.obs["Type"] = "CRC"

adata_T.obs.head()



###############################

adata_N = anndata.read("/home/pr186/Documents/IRP/FINAL/CRC_RAW/KUL28_N/adata.h5ad")

adata_N.var["gene_id"] = adata_N.var.index.values

#adata_N.var.index = adata_N.var["gene_name"]

adata_N.var_names_make_unique()

adata_N.obs["Type"] = "CONTOL"
adata_N


### BASIC FILTERING
    
# Removes genes with 0 umi counts
adata_T = adata_T[:, np.asarray(adata_T.X.sum(axis=0)).reshape(-1) > 0]
adata_N = adata_N[:, np.asarray(adata_N.X.sum(axis=0)).reshape(-1) > 0]
  
#adata.var_names_make_unique()
    
sc.pp.filter_cells(adata_T, min_genes=200)
sc.pp.filter_genes(adata_T, min_cells=3)

sc.pp.filter_cells(adata_N, min_genes=200)
sc.pp.filter_genes(adata_N, min_cells=3)
 
    
# Create a mask to filter out cells with more than 6500 genes
adata_T = adata_T[adata_T.obs.n_genes < 6500, :]
adata_N = adata_N[adata_N.obs.n_genes < 6500, :]

adata_T
adata_N


sc.pp.normalize_total(adata_T, inplace=True)
sc.pp.log1p(adata_T)
sc.pp.highly_variable_genes(adata_T, flavor="seurat", n_top_genes=2000, inplace=True)

sc.pp.normalize_total(adata_N, inplace=True)
sc.pp.log1p(adata_N)
sc.pp.highly_variable_genes(adata_N, flavor="seurat", n_top_genes=2000, inplace=True)


adata_T.obs["barcode_id"] = "KUL28-T_" + adata_T.obs.index.values
adata_N.obs["barcode_id"] = "KUL28-N_" + adata_N.obs.index.values


metadata = pd.read_csv("/home/pr186/Documents/IRP/CRC_FINAL/annotation.txt", header = 0, sep = "\t")
metadata.index = metadata.Index
adata_T.obs["celltype"] = adata_T.obs.barcode_id.map(metadata["Cell_type"])
adata_T.obs["cellsubtype"] = adata_T.obs.barcode_id.map(metadata["Cell_subtype"])

adata_N.obs["celltype"] = adata_N.obs.barcode_id.map(metadata["Cell_type"])
adata_N.obs["cellsubtype"] = adata_N.obs.barcode_id.map(metadata["Cell_subtype"])

adata_N.obs.head()

###########################INTEGRATION
adatas = [adata_T, adata_N]

adatas_cor = scanorama.correct_scanpy(adatas, return_dimred=True)

adatas_cor

adata_spatial = adatas_cor[0].concatenate(
    adatas_cor[1],
    batch_key="Type",
    uns_merge="unique",
    batch_categories = ["CRC", "CONTROL"]
)

adata_spatial

adata_spatial.obs.head()
adata_spatial.var.head()


adata_spatial.obs.index = adata_spatial.obs["barcode_id"]


sc.pp.neighbors(adata_spatial, use_rep="X_scanorama")
sc.tl.umap(adata_spatial)
sc.tl.leiden(adata_spatial)

ann_T = anndata.read("/home/pr186/Documents/IRP/FINAL/CRC_RAW/KUL28_T/KUL28_T_CELLO.h5ad")
ann_N = anndata.read("/home/pr186/Documents/IRP/FINAL/CRC_RAW/KUL28_N/KUL28_N_CELLO.h5ad")

ann_T.obs["barcode_id"] = "KUL28-T_" + ann_T.obs.index.values
ann_N.obs["barcode_id"] = "KUL28-N_" + ann_N.obs.index.values


## ADDING PREFIX TO DISTINGUISH CELL TYPES
ann_T.obs['CELLO'] = "CRC " + ann_T.obs["Most specific cell type"].astype(str)
ann_N.obs['CELLO'] = "CONTROL " + ann_N.obs["Most specific cell type"].astype(str)

ann_T.obs.index = ann_T.obs["barcode_id"]
ann_N.obs.index = ann_N.obs["barcode_id"]

## CREATING METADATAFILE FOR ANNOTATION
#ann_T.obs.to_csv("ann_T.csv")
#ann_N.obs.to_csv("ann_N.csv")

ann_N
ann_T

df_T = ann_T.obs
df_N = ann_N.obs

frames = [df_T, df_N]

annotation = pd.concat(frames)
annotation

annotation.index = annotation.barcode_id

## ADDING CELLO ANNOTATION TO INTEGRATED DATA
adata_spatial.obs["CELLO"] =  adata_spatial.obs.barcode_id.map(annotation["CELLO"])


adata_spatial


sc.pl.umap(adata_spatial, color=["Type"], palette=sc.pl.palettes.default_20, save = '_1.png')
sc.pl.umap(adata_spatial, color=["leiden"], palette=sc.pl.palettes.default_20, save = '_2.png')
sc.pl.umap(adata_spatial, color=["celltype"], palette=sc.pl.palettes.default_20, save = '_3.png')
sc.pl.umap(adata_spatial, color=["cellsubtype"], palette=sc.pl.palettes.godsnot_102, save = '_4.png')
sc.pl.umap(adata_spatial, color=["CELLO"], palette=sc.pl.palettes.godsnot_102, save = '_5.png')


pd.set_option('display.max_columns', None)



###########gene visualisation
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
        sc.pl.umap(adata_spatial, gene_symbols = "gene_name", color = genes_TARGET[gene], save = '_' + genes_TARGET[gene] + '_TARGET.png')
        genes_TARGET_PR.append(genes_TARGET[gene])
    except IndexError:
        print(genes_TARGET[gene] + " not present")
        
sc.pl.umap(adata_spatial, gene_symbols = "gene_name", color = genes_TARGET_PR, save = '_TARGET_ALL.png')

for gene in range(len(genes_TARGET2)):
    try:
        sc.pl.umap(adata_spatial, gene_symbols = "gene_name", color = genes_TARGET2[gene], save = '_' + genes_TARGET2[gene] + '_TARGET2.png')
        genes_TARGET2_PR.append(genes_TARGET2[gene])
    except IndexError:
        print(genes_TARGET2[gene] + " not present")
        
sc.pl.umap(adata_spatial, gene_symbols = "gene_name", color = genes_TARGET2_PR, save = '_TARGET2_ALL.png')


for gene in range(len(genes_PROLINE)):
    try:
        sc.pl.umap(adata_spatial, gene_symbols = "gene_name", color = genes_PROLINE[gene], save = '_' + genes_PROLINE[gene] + '_PROLINE.png')
        genes_PROLINE_PR.append(genes_PROLINE[gene])
    except IndexError:
        print(genes_PROLINE[gene] + " not present")
        
sc.pl.umap(adata_spatial, gene_symbols = "gene_name", color = genes_PROLINE_PR, save = '_PROLINE_ALL.png')
        
for gene in range(len(genes_GLYCOLYSIS)):
    try:
        sc.pl.umap(adata_spatial, gene_symbols = "gene_name", color = genes_GLYCOLYSIS[gene], save = '_' + genes_GLYCOLYSIS[gene] + '_GLYCOLYSIS.png')
        genes_GLYCOLYSIS_PR.append(genes_GLYCOLYSIS[gene])
    except IndexError:
        print(genes_GLYCOLYSIS[gene] + " not present")

sc.pl.umap(adata_spatial, gene_symbols = "gene_name", color = genes_GLYCOLYSIS_PR, save = '_GLYCOLYSIS_ALL.png')





ax = sc.pl.dotplot(adata_spatial, genes_TARGET_PR, gene_symbols = "gene_name", groupby='Type', use_raw = False, save = '_leiden_TARGET_25.png')
ax = sc.pl.dotplot(adata_spatial, genes_TARGET2_PR, gene_symbols = "gene_name", groupby='Type', use_raw = False, save = '_leiden_TARGET2_25.png')
ax = sc.pl.dotplot(adata_spatial, genes_PROLINE_PR, gene_symbols = "gene_name", groupby='Type', use_raw = False, save = '_leiden_PROLINE_26.png')
ax = sc.pl.dotplot(adata_spatial, genes_GLYCOLYSIS_PR, gene_symbols = "gene_name", groupby='Type', use_raw = False, save = '_leiden_GLYCOLYSIS_27.png')

ax = sc.pl.dotplot(adata_spatial, genes_TARGET_PR, gene_symbols = "gene_name", groupby='leiden', use_raw = False, save = '_leiden2_TARGET_25.png')
ax = sc.pl.dotplot(adata_spatial, genes_TARGET2_PR, gene_symbols = "gene_name", groupby='leiden', use_raw = False, save = '_leiden2_TARGET2_25.png')
ax = sc.pl.dotplot(adata_spatial, genes_PROLINE_PR, gene_symbols = "gene_name", groupby='leiden', use_raw = False, save = '_leiden2_PROLINE_26.png')
ax = sc.pl.dotplot(adata_spatial, genes_GLYCOLYSIS_PR, gene_symbols = "gene_name", groupby='leiden', use_raw = False, save = '_leiden2_GLYCOLYSIS_27.png')

ax = sc.pl.dotplot(adata_spatial, genes_TARGET_PR, gene_symbols = "gene_name", groupby='celltype', use_raw = False, save = '_AUTHORS1_TARGET_25.png')
ax = sc.pl.dotplot(adata_spatial, genes_TARGET2_PR,  gene_symbols = "gene_name", groupby='celltype', use_raw = False, save = '_AUTHORS12_TARGET2_25.png')
ax = sc.pl.dotplot(adata_spatial, genes_PROLINE_PR, gene_symbols = "gene_name",  groupby='celltype', use_raw = False, save = '_AUTHORS1_PROLINE_26.png')
ax = sc.pl.dotplot(adata_spatial, genes_GLYCOLYSIS_PR,  gene_symbols = "gene_name", groupby='celltype', use_raw = False, save = '_AUTHORS1_GLYCOLYSIS_27.png')

ax = sc.pl.dotplot(adata_spatial, genes_TARGET_PR, gene_symbols = "gene_name", groupby='cellsubtype', use_raw = False, save = '_AUTHORS2_TARGET_25.png')
ax = sc.pl.dotplot(adata_spatial, genes_TARGET2_PR,  gene_symbols = "gene_name", groupby='cellsubtype', use_raw = False, save = '_AUTHORS2_TARGET2_25.png')
ax = sc.pl.dotplot(adata_spatial, genes_PROLINE_PR, gene_symbols = "gene_name",  groupby='cellsubtype', use_raw = False, save = '_AUTHORS2_PROLINE_26.png')
ax = sc.pl.dotplot(adata_spatial, genes_GLYCOLYSIS_PR,  gene_symbols = "gene_name", groupby='cellsubtype', use_raw = False, save = '_AUTHORS2_GLYCOLYSIS_27.png')

ax = sc.pl.dotplot(adata_spatial, genes_TARGET_PR, gene_symbols = "gene_name", groupby='CELLO', use_raw = False, save = '_CELLO_TARGET_25.png')
ax = sc.pl.dotplot(adata_spatial, genes_TARGET2_PR,  gene_symbols = "gene_name", groupby='CELLO', use_raw = False, save = '_CELLO_TARGET2_25.png')
ax = sc.pl.dotplot(adata_spatial, genes_PROLINE_PR, gene_symbols = "gene_name",  groupby='CELLO', use_raw = False, save = '_CELLO_PROLINE_26.png')
ax = sc.pl.dotplot(adata_spatial, genes_GLYCOLYSIS_PR,  gene_symbols = "gene_name", groupby='CELLO', use_raw = False, save = '_CELLO_GLYCOLYSIS_27.png')























