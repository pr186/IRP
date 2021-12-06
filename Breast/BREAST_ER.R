library(Seurat)
library(cowplot)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tibble)
library(purrr)
library(stringr)
library(gplots)
library(ggpubr)
library(Matrix)

#devtools::install_github('ZJUFanLab/scCATCH')
library(scCATCH)

INPUT <- "/home/pr186/Documents/IRP/BREAST/GSE161529_RAW"
OUTPUT <- "/home/pr186/Documents/IRP/FINAL/BREAST/N_ER+"


setwd(INPUT)

#######################################CONTROL1################################################################


CONTROL_01_matrix <- ReadMtx(mtx = "GSM4909254_N-PM0019-Total-matrix.mtx.gz", 
                             features = "features.tsv", 
                             cells = "GSM4909254_N-PM0019-Total-barcodes.tsv.gz")
CONTROL_01_matrix

CONTROL_01_object <- CreateSeuratObject(counts = CONTROL_01_matrix, project = "CONTROL_01", min.cells = 3, min.features = 200)
CONTROL_01_object$Type <- "CONTROL"

setwd(OUTPUT)
png(file = "GENE_COUNT_BF_FILTER.png")
VlnPlot(CONTROL_01_object, features = "nCount_RNA")
dev.off()

CONTROL_01_object <- subset(CONTROL_01_object, subset = nCount_RNA < 6500)
CONTROL_01_object <- NormalizeData(CONTROL_01_object)
CONTROL_01_object <- FindVariableFeatures(CONTROL_01_object, selection.method = "vst", nfeatures = 2000)
CONTROL_01_object

top10_CTRL_1 <- head(VariableFeatures(CONTROL_01_object), 10)

plot1 <- VariableFeaturePlot(CONTROL_01_object)
plot2 <- LabelPoints(plot = plot1, points = top10_CTRL_1, repel = TRUE)

png(file = "CTRL_1_VARIABLE.png", width = 1200, height = 700)
plot1 + plot2
dev.off()

#######################################CONTROL2################################################################

setwd(INPUT)
CONTROL_02_matrix <- ReadMtx(mtx = "GSM4909257_N-PM0095-Total-matrix.mtx.gz", 
                             features = "features.tsv", 
                             cells = "GSM4909257_N-PM0095-Total-barcodes.tsv.gz")
CONTROL_02_matrix

CONTROL_02_object <- CreateSeuratObject(counts = CONTROL_02_matrix, project = "CONTROL_02", min.cells = 3, min.features = 200)
CONTROL_02_object$Type <- "CONTROL"

setwd(OUTPUT)
png(file = "GENE_COUNT_2_BF_FILTER.png")
VlnPlot(CONTROL_02_object, features = "nCount_RNA")
dev.off()

CONTROL_02_object <- subset(CONTROL_02_object, subset = nCount_RNA < 6500)

CONTROL_02_object <- NormalizeData(CONTROL_02_object)
CONTROL_02_object <- FindVariableFeatures(CONTROL_02_object, selection.method = "vst", nfeatures = 2000)
CONTROL_02_object


top10_CTRL_2 <- head(VariableFeatures(CONTROL_02_object), 10)

plot3 <- VariableFeaturePlot(CONTROL_02_object)
plot4 <- LabelPoints(plot = plot3, points = top10_CTRL_2, repel = TRUE)

png(file = "CTRL_2_VARIABLE.png", width = 1200, height = 700)
plot3 + plot4
dev.off()


#######################################TUMOR################################################################
# ER


setwd(INPUT)
CANCER_01_matrix <- ReadMtx(mtx = "GSM4909315_ER-MH0167-T-matrix.mtx.gz", 
                            features = "features.tsv", 
                            cells = "GSM4909315_ER-MH0167-T-barcodes.tsv.gz")
CANCER_01_matrix

CANCER_01_object <- CreateSeuratObject(counts = CANCER_01_matrix, project = "ER+_TUMOR", min.cells = 3, min.features = 200)
CANCER_01_object$Type <- "ER+_TUMOR"

setwd(OUTPUT)
png(file = "GENE_COUNT_3_BF_FILTER.png")
VlnPlot(CANCER_01_object, features = "nCount_RNA")
dev.off()

CANCER_01_object <- subset(CANCER_01_object, subset = nCount_RNA < 6500)

CANCER_01_object <- NormalizeData(CANCER_01_object)
CANCER_01_object <- FindVariableFeatures(CANCER_01_object, selection.method = "vst", nfeatures = 2000)
CANCER_01_object


top10_TN_1 <- head(VariableFeatures(CANCER_01_object), 10)

plot5 <- VariableFeaturePlot(CANCER_01_object)
plot6 <- LabelPoints(plot = plot5, points = top10_TN_1, repel = TRUE)

png(file = "TUMOR_1_VARIABLE.png", width = 1200, height = 700)
plot5 + plot6
dev.off()



#######################################INTEGRATION&CLUSTERING################################################################


combined.anchors <- FindIntegrationAnchors(object.list = list(CONTROL_01_object,
                                                              CONTROL_02_object,
                                                              CANCER_01_object), dims = 1:20)

combined <- IntegrateData(anchorset = combined.anchors, dims = 1:20)

DefaultAssay(combined) <- "integrated"

#all.genes = rownames(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)

png(file = "PCA_ELBOW.png")
ElbowPlot(combined)
dev.off()

ElbowPlot(combined)
#PCAPlot(object = combined)
#PCHeatmap(combined)


combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- RunTSNE(combined, reduction = "pca", dims = 1:20)


library(dittoSeq)

png(file = "proportions.png", width = 700, height = 400)
dittoBarPlot(combined_annotated, "Type", group.by = "seurat_clusters", color.panel = c("#F8766D", "#00BFC4"))
dev.off()


#######################################VISUALISATION################################################################

#TSNE
p1 <- DimPlot(combined, reduction = "tsne", group.by = "Type", pt.size = 0.5)
p2 <- DimPlot(combined, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 6)
plot_grid(p1, p2)

png(file = "TSNE_SEP.png", width = 900, height = 400)
DimPlot(combined, reduction = "tsne",group.by = "Type", split.by = "Type", pt.size = 0.5)
dev.off()

png(file = "TSNE.png", width = 900, height = 400)
plot_grid(p1, p2)
dev.off()

png(file = "TSNE_2.png", width = 1200, height = 700)
plot_grid(p1, p2)
dev.off()

#UMAP
p3 <- DimPlot(combined, reduction = "umap", group.by = "Type", pt.size = 0.5)
p4 <- DimPlot(combined, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6)
plot_grid(p3, p4)

png(file = "UMAP_SEP.png", width = 800, height = 400)
DimPlot(combined, reduction = "umap",group.by = "Type", split.by = "Type", pt.size = 0.5)
dev.off()

png(file = "UMAP.png", width = 900, height = 400)
plot_grid(p3, p4)
dev.off()

png(file = "UMAP_2.png", width = 1200, height = 700)
plot_grid(p3, p4)
dev.off()


DefaultAssay(combined) <- "RNA"

TARGET <- c("SIX1", "SIX2", "PAX6", "DACH1", "DACH2", "EYA1", "EYA2", "EYA3")

png(file = "VLN_TARGET.png", width = 1200, height = 700)
VlnPlot(combined, features = c(TARGET), split.by = "Type", group.by = "Type")
dev.off()

png(file = "DOT_GROUP_TARGET.png", width = 600, height = 300)
DotPlot(combined, features = rev(TARGET), cols = c("blue", "red"), dot.scale = 10, group.by = "Type")
dev.off()

png(file = "VLN_TARGET_SMALL.png", width = 700, height = 500)
VlnPlot(combined, features = c(TARGET), split.by = "Type", group.by = "Type")
dev.off()

png(file = "DOT_TARGET.png", width = 1200, height = 700)
DotPlot(combined, features = rev(TARGET), cols = c("blue", "red"), dot.scale = 10, split.by = "Type")
dev.off()

png(file = "DOT_TARGET_SMALL.png", width = 700, height = 500)
DotPlot(combined, features = rev(TARGET), cols = c("blue", "red"), dot.scale = 10, split.by = "Type")
dev.off()

png(file = "FEATURE_TARGET_SMALL.png", width = 500, height = 1200)
FeaturePlot(combined, TARGET, split.by = "Type", pt.size = 1)
dev.off()

##FIND MARKERS WITH SEURAT
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
combined.markers %>%
  group_by(cluster) %>%
  slice_max(n=5, order_by = avg_log2FC) %>%
  print(n=30) -> TOP5_MARKERS

write.csv(combined.markers, "MARKERS.csv", quote = F)
write.csv(TOP5_MARKERS, "MARKERS_TOP5.csv", quote = F)


ann_1 <- scCATCH(combined.markers,
                 species = "Human",
                 tissue = c("Breast", "Mammary epithelium"))

ann_C <- scCATCH(combined.markers,
                 species = "Human",
                 cancer = "Breast Cancer",
                 tissue = "Breast")

##FIND MARKERS WITH scCHATCH
combined.markers_scCATCH <- findmarkergenes(combined,
                                            species = "Human",
                                            cluster = 'All',
                                            match_CellMatch = TRUE,
                                            cancer = "Breast Cancer",
                                            tissue = "Breast",
                                            cell_min_pct = 0.25,
                                            logfc = 0.25,
                                            pvalue = 0.05)


ann_2 <- scCATCH(combined.markers_scCATCH,
                 species = "Human",
                 tissue = c("Breast", "Mammary epithelium"))

ann_C_2 <- scCATCH(combined.markers_scCATCH,
                   species = "Human",
                   cancer = "Breast Cancer",
                   tissue = "Breast")


write.csv(ann_1, "ann_1_MARKERS.csv", quote = F)
write.csv(ann_C, "ann_C_MARKERS.csv", quote = F)
write.csv(ann_2, "ann_2_MARKERS.csv", quote = F)
write.csv(ann_C_2, "ann_C_2_MARKERS.csv", quote = F)

combined_annotated <- RenameIdents(combined, 
                                   '0' = "Helper T Cells",
                                   '1' = "Luminal Progenitor Cells",
                                   '2' = "Progenitor Cells",
                                   '3' = "Helper T Cells",
                                   '4' = "Epithelial Cells",
                                   '5' = "Luminal Epithelial Cells",
                                   '6' = "Progenitor Cells",
                                   '7' = "Basal Epithelial Cells",
                                   '8' = "Progenitor Cells",
                                   '9' = "Progenitor Cells",
                                   '10' = "Basal Epithelial Cells",
                                   '11' = "Progenitor Cells",
                                   '12' = "Luminal Progenitor Cells",
                                   '13' = "Cancer Stem Cells",
                                   '14' = "Epithelial Cells",
                                   '15' = "B Cells")



png(file = "DOT_ANN_TARGET_SMALL.png", width = 700, height = 500)
DotPlot(combined_annotated, features = rev(TARGET), cols = c("blue", "red"), dot.scale = 10, split.by = "Type")
dev.off()

##TSNE
p6 <- DimPlot(combined_annotated, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 6) + NoLegend()
plot_grid(p6)
p7 <- DimPlot(combined_annotated, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 6)
p11 <- DimPlot(combined_annotated, reduction = "tsne", label = FALSE, pt.size = 0.5, label.size = 6)

png(file = "TSNE_ANN.png", width = 1000, height = 500)
plot_grid(p2, p6)
dev.off()

png(file = "TSNE_ANN_2.png", width = 700, height = 500)
plot_grid(p7)
dev.off()

png(file = "TSNE_ANN_3.png", width = 700, height = 500)
plot_grid(p11)
dev.off()

#UMAP
p8 <- DimPlot(combined_annotated, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6) + NoLegend()
plot_grid(p8)
p9 <- DimPlot(combined_annotated, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6)
p10 <- DimPlot(combined_annotated, reduction = "umap", label = FALSE, pt.size = 0.5, label.size = 6)


png(file = "UMAP_ANN_3.png", width = 1000, height = 500)
plot_grid(p4, p8)
dev.off()

png(file = "UMAP_ANN_4.png", width = 700, height = 500)
plot_grid(p9)
dev.off()

png(file = "UMAP_ANN_5.png", width = 700, height = 500)
plot_grid(p10)
dev.off()

png(file = "UMAP_ANN_6.png", width = 1500, height = 500)
plot_grid(p4, p10)
dev.off()


TARGET2 <-  c("SIX1", "SIX2", "SIX3", "SIX4", "SIX5", "SIX6", "PAX4", "PAX6", 
              "DACH1", "DACH2", "EYA1", "EYA2", "EYA3", "EYA4", "EYA5", "EYA6")

png(file = "VLN_TARGET2.png", width = 1200, height = 700)
VlnPlot(combined, features = c(TARGET2), split.by = "Type", group.by = "Type")
dev.off()

png(file = "DOT_GROUP_TARGET2.png", width = 600, height = 300)
DotPlot(combined_annotated, features = rev(TARGET2), cols = c("blue", "red"), dot.scale = 10, group.by = "Type")
dev.off()

png(file = "VLN_TARGET2_SMALL.png", width = 700, height = 500)
VlnPlot(combined, features = c(TARGET2), split.by = "Type", group.by = "Type")
dev.off()

png(file = "DOT_TARGET2.png", width = 1200, height = 700)
DotPlot(combined, features = rev(TARGET2), cols = c("blue", "red"), dot.scale = 10, split.by = "Type")
dev.off()

png(file = "DOT_ANN_TARGET2_SMALL.png", width = 1200, height = 500)
DotPlot(combined_annotated, features = rev(TARGET2), cols = c("blue", "red"), dot.scale = 10, split.by = "Type")
dev.off()




png(file = "TSNE_SIX1.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "tsne", features = "SIX1", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "TSNE_SIX2.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "tsne", features = "SIX2", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "TSNE_SIX3.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "tsne", features = "SIX3", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "TSNE_SIX4.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "tsne", features = "SIX4", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "TSNE_SIX5.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "tsne", features = "SIX5", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "TSNE_SIX6.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "tsne", features = "SIX6", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "TSNE_PAX4.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "tsne", features = "PAX4", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "TSNE_PAX6.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "tsne", features = "PAX6", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "TSNE_DACH1.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "tsne", features = "DACH1", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "TSNE_DACH2.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "tsne", features = "DACH2", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "TSNE_EYA1.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "tsne", features = "EYA1", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "TSNE_EYA2.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "tsne", features = "EYA2", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "TSNE_EYA3.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "tsne", features = "EYA3", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "TSNE_EYA4.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "tsne", features = "EYA4", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

###UMAP

png(file = "UMAP_SIX1.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "umap", features = "SIX1", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "UMAP_SIX2.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "umap", features = "SIX2", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "UMAP_SIX3.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "umap", features = "SIX3", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "UMAP_SIX4.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "umap", features = "SIX4", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "UMAP_SIX5.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "umap", features = "SIX5", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "UMAP_SIX6.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "umap", features = "SIX6", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "UMAP_PAX4.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "umap", features = "PAX4", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "UMAP_PAX6.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "umap", features = "PAX6", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "UMAP_DACH1.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "umap", features = "DACH1", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "UMAP_DACH2.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "umap", features = "DACH2", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "UMAP_EYA1.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "umap", features = "EYA1", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "UMAP_EYA2.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "umap", features = "EYA2", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "UMAP_EYA3.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "umap", features = "EYA3", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()

png(file = "UMAP_EYA4.png", width = 700, height = 350)
FeaturePlot(combined, reduction = "umap", features = "EYA4", max.cutoff = 3, cols = c("grey", "red"), pt.size = 2, split.by = "Type")
dev.off()







