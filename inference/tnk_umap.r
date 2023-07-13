load("datasets/GSE131907-exp-mat-by-cell-types/CRC_primarytumor_TNK_intersectfeats_Seuratobj.RData")
load("datasets/GSE131907-exp-mat-by-cell-types/TNK_sprase.RData")

library(Seurat)
library(harmony)

lc_tnk_input <- query_df_sparse
colnames(lc_tnk_input) <- query_df_colnames
rownames(lc_tnk_input) <- query_df_rownames
lc_gene <- query_df_rownames

#####crc genes
crc_tnk_input <- CRC_TNK_Seuratobj@assays$RNA@counts
crc_gene <- rownames(crc_tnk_input)

#####get the intersection of gene
intersect_gene <- intersect(lc_gene, crc_gene)

#####get the gene index
crc_gene_index <- match(intersect_gene, crc_gene)
lc_gene_index <- match(intersect_gene, lc_gene)

#####prepare the crc expression of intersected genes
crc_sc_intersect <- crc_tnk_input[crc_gene_index, ]
lc_sc_intersect <- lc_tnk_input[lc_gene_index, ]

#####combine the crc dataset and breast dataset with same intersected genes
#####the dimension of the merged_sc_intersect is the number of intersected 
#####gene*combined cells(crc cells + breast cells)
merged_sc_intersect <- cbind(crc_sc_intersect, lc_sc_intersect)


####the input of the seurat should be genes*cells
#####use the merged expression matrix to create the seurat object

#####add the meta data to show the cell sample origin
crc_intersect_ncells <- ncol(crc_sc_intersect)
lc_intersect_ncells <- ncol(lc_sc_intersect)
total_intersect_ncells <- lc_intersect_ncells + crc_intersect_ncells

samples <- matrix(
    c(rep("crc", crc_intersect_ncells), rep("lc", lc_intersect_ncells)),
    total_intersect_ncells,
    1)
samples <- data.frame(samples)
rownames(samples) <- colnames(merged_sc_intersect)

####extract the count matrix of GEM for CRC data from the trained nHDP model directly
load("model/R_nHDPresults_binary_tnk_three_layer_20230620.RData")
load("inference/lung_cancer_TNK.RData")

# Extracting value matrix
crc_value_matrix <- nHDP_trained_mb$count_matrix
lc_value_matrix <- inferred_data$count_matrix

#####you need to set the gem value for breast cancer cells to be 0 (crc tnk gems)
n_dim <- dim(merged_sc_intersect)[2]

add_crc_tnk_gem_matrix <- matrix(0, n_dim, 85)
add_crc_tnk_gem_matrix[
    1: dim(crc_value_matrix)[2],
    ] <- t(crc_value_matrix)
colnames(add_crc_tnk_gem_matrix) <- paste("crc_tnk", seq(1, 85, 1), sep = "")

#####you need to set the gem value for crc cancer cells to be 0 (breast tnk gems)
add_lc_tnk_gem_matrix <- matrix(0, n_dim, 85)
add_lc_tnk_gem_matrix[
    (dim(crc_value_matrix)[2] + 1): n_dim,
    ] <- t(lc_value_matrix)
colnames(add_lc_tnk_gem_matrix) <- paste("lc_tnk", seq(1, 85, 1), sep = "")

#####this is where you incorporate the gem expression value into the seurat object
metadata <- cbind(samples, add_crc_tnk_gem_matrix, add_lc_tnk_gem_matrix)

# Creating a seurat object with merged_crc_lc
merged_crc_lc <- CreateSeuratObject(
    counts = merged_sc_intersect,
    project = "crc_lc_sc_dataset",
    min.cells = 3,
    min.features = 200,
    meta.data = metadata)

# Normalizing the data
merged_crc_lc <- NormalizeData(
    merged_crc_lc,
    normalization.method = "LogNormalize",
    scale.factor = 10000)

# Identifying highly variable features
merged_crc_lc <- FindVariableFeatures(
    merged_crc_lc,
    selection.method = "vst",
    nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(merged_crc_lc), 10)

# Scaling the data
all.genes <- rownames(merged_crc_lc)
merged_crc_lc <- ScaleData(merged_crc_lc, features = all.genes)

# Run Harmony
merged_crc_lc <- RunHarmony(merged_crc_lc, "samples")

VizDimLoadings(merged_crc_lc, dims = 1:2, reduction = "pca")
DimPlot(merged_crc_lc, reduction = "pca")
DimHeatmap(merged_crc_lc, dims = 1, cells = 500, balanced = TRUE)

# Heatmap of the most important features for each PC
DimHeatmap(merged_crc_lc, dims = 1:15, cells = 500, balanced = TRUE)

# Cluster the cells
merged_crc_lc <- FindNeighbors(merged_crc_lc, dims = 1:10)
merged_crc_lc <- FindClusters(merged_crc_lc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
# head(Idents(merged_crc_lc), 5)

# Run non-linear dimension reduction
merged_crc_lc <- RunUMAP(merged_crc_lc, dims = 1:30, reduction = "harmony")

###### Combined UPMAP is ready! ######


#####Plot UMAP using for loop (all the crc gems)
for (i in 1:dim(crc_value_matrix)[2]){
  jpeg(
    paste("inference/crc_tnk_umap/", "crc_tnk_cell_feature_gem_", i, ".jpeg", sep = ""),
    res = 300,
    width = 2000,
    height = 2000)
  print(FeaturePlot(
    object = merged_crc_lc, features = paste("crc_tnk", i, sep = ""),
    cols = c("grey85", "red"), raster = FALSE))
  dev.off()
}

#####example of all GEMs plot using for loop (all the lc gems)
for (i in 1:dim(lc_value_matrix)[2]){
  jpeg(
    paste("inference/lc_tnk_umap/", "lc_tnk_cell_feature_gem_", i, ".jpeg", sep = ""),
    res = 300,
    width = 2000,
    height = 2000)
  print(FeaturePlot(
    object = merged_crc_lc,
    features = paste("lc_tnk", i, sep = ""),
    cols = c("grey85", "red"), raster = FALSE))
  dev.off()
}