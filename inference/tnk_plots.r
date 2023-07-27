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

######samples is creating a single column with #crc cells + #lc cells
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

#####you need to set the gem value for lung cancer cancer cells to be 0 (crc tnk gems)
n_dim <- dim(merged_sc_intersect)[2]

# metadata requires rows as cells. So we are creating a matrix with n_dim rows
# dim(add_crc_tnk_gem_matrix) = 337787, 85
add_crc_tnk_gem_matrix <- matrix(0, n_dim, 85)
# filling first 246560 rows, (246560, 85) with crc value matrix (85, 246560)
# thus need to transpose
# NOTE: lc cell rows remain 0
add_crc_tnk_gem_matrix[
    1: dim(crc_value_matrix)[2],
    ] <- t(crc_value_matrix)
colnames(add_crc_tnk_gem_matrix) <- paste("crc_tnk", seq(1, 85, 1), sep = "")


######to do: set the maximum value to 200 and update the add_crc_tnk_gem_matrix
###### so that lower expression values in the umap are also visible clearly
add_crc_tnk_gem_matrix[add_crc_tnk_gem_matrix > 200] <- 200


######to do: only keep the GEMs expressed in more than 5% of cells
print(sprintf("add_crc_tnk_gem_matrix dimensions before filtering GEMs: %d", dim(add_crc_tnk_gem_matrix)))

# add_crc_tnk_gem_matrix has crc + lc cells. if we calculate nrows of this,
# we might calculate 5% of the total, and lose more gems than expected
crc_ncells <- ncol(crc_value_matrix)
threshold <- crc_ncells * 0.05
col_to_rm <- c()
for (col in colnames(add_crc_tnk_gem_matrix)){
    # count of genes with non-zero values less than equal to t will be invalid
    if (sum(add_crc_tnk_gem_matrix[, col] > 0) <= threshold) {
        col_to_rm <- c(col_to_rm, col)
    }
}
add_crc_tnk_gem_matrix <- add_crc_tnk_gem_matrix[, !colnames(add_crc_tnk_gem_matrix) %in% col_to_rm]
print(sprintf("add_crc_tnk_gem_matrix dimensions after filtering GEMs: %d", dim(add_crc_tnk_gem_matrix)))


#####you need to set the gem value for crc cancer cells to be 0 (breast tnk gems)
add_lc_tnk_gem_matrix <- matrix(0, n_dim, 85)
add_lc_tnk_gem_matrix[
    (dim(crc_value_matrix)[2] + 1): n_dim,
    ] <- t(lc_value_matrix)
colnames(add_lc_tnk_gem_matrix) <- paste("lc_tnk", seq(1, 85, 1), sep = "")


######to do: set the maximum value to 200 and update the add_breast_tnk_gem_matrix
add_lc_tnk_gem_matrix[add_lc_tnk_gem_matrix > 200] <- 200


######to do: only keep the GEMs expressed in more than 5% of cells 
######(save this for later. It seems that the new code leads to several 
######GEMs with relatively small values)
print(sprintf("add_lc_tnk_gem_matrix dimensions before filtering GEMs: %d", dim(add_lc_tnk_gem_matrix)))
lc_ncells <- ncol(lc_value_matrix)
threshold <- lc_ncells * 0.05
col_to_rm <- c()
for (col in colnames(add_lc_tnk_gem_matrix)){
    # count of genes with non-zero values less than equal to t will be invalid
    if (sum(add_lc_tnk_gem_matrix[, col] > 0) <= threshold) {
        col_to_rm <- c(col_to_rm, col)
    }
}
add_lc_tnk_gem_matrix <- add_lc_tnk_gem_matrix[, !colnames(add_lc_tnk_gem_matrix) %in% col_to_rm]
print(sprintf("add_lc_tnk_gem_matrix dimensions after filtering GEMs: %d", dim(add_lc_tnk_gem_matrix)))


#####this is where you incorporate the gem expression value into the seurat object
metadata <- cbind(samples, add_crc_tnk_gem_matrix, add_lc_tnk_gem_matrix)

# Saving col names
crc_tnk_gem_names <- colnames(add_crc_tnk_gem_matrix)
lc_tnk_gem_names <- colnames(add_lc_tnk_gem_matrix)

rm(add_crc_tnk_gem_matrix, crc_tnk_input, add_lc_tnk_gem_matrix, crc_sc_intersect, CRC_TNK_Seuratobj, crc_value_matrix, lc_sc_intersect, lc_tnk_input, lc_value_matrix, query_df_sparse, X)
gc()

# Creating a seurat object with merged_crc_lc
merged_crc_lc <- CreateSeuratObject(
    # counts has genes * cells
    counts = merged_sc_intersect,
    project = "crc_lc_sc_dataset",
    min.cells = 3,
    min.features = 200,
    # cells * features
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


# Run non-linear dimension reduction
merged_crc_lc <- RunUMAP(merged_crc_lc, dims = 1:30, reduction = "harmony")


# Plot 1: Merged DimPlot
jpeg("inference/results/merged_dimplot.jpeg",
  res = 300,
  width = 2000,
  height = 2000)
print(DimPlot(merged_crc_lc, reduction = "harmony", raster = FALSE))
dev.off()


# Plot 2: FP samples (mergerd crc & lc UMAP)
jpeg("inference/results/FeaturePlots/samples.jpeg",
  res = 300,
  width = 2000,
  height = 2000)
print(DimPlot(merged_crc_lc, group.by = "samples", raster = FALSE))
dev.off()


# Plot 3: Feature Plots for particular T-cell marker genes
features <- list("CD8A", "CD4", "MKI67", "CXCL13", "PDCD1", "CTLA4", "CCR7", "CXCR3", "CXCR5", "ICOS")
for (feat in features){
  jpeg(
    paste("inference/results/FeaturePlots/", feat, ".jpeg", sep = ""),
    res = 300,
    width = 2000,
    height = 2000)
  print(FeaturePlot(
  object = merged_crc_lc,
  features = feat,
  cols = c("grey85", "red"),
  raster = FALSE))
dev.off()
}


# Plot 4: UMAP for selected GEMs
#####Plot UMAP using for loop (all the crc gems)
for (gem_name in colnames(add_crc_tnk_gem_matrix)){
  jpeg(
    paste("inference/results/crc_tnk_umap/", gem_name, ".jpeg", sep = ""),
    res = 300,
    width = 2000,
    height = 2000)
  print(FeaturePlot(
    object = merged_crc_lc, features = gem_name,
    cols = c("grey85", "red"), raster = FALSE))
  dev.off()
}

#####example of all GEMs plot using for loop (all the lc gems)
for (gem_name in colnames(add_lc_tnk_gem_matrix)){
  jpeg(
    paste("inference/results/lc_tnk_umap/", gem_name, ".jpeg", sep = ""),
    res = 300,
    width = 2000,
    height = 2000)
  print(FeaturePlot(
    object = merged_crc_lc,
    features = gem_name,
    cols = c("grey85", "red"), raster = FALSE))
  dev.off()
}


# Plot 5: Wordcloud of top genes associated with a GEM

# Plot 6: Enriched Signaliing Pathway With a GEM

# Plot 7: Pick some GEMs to compare the difference between crc and breast