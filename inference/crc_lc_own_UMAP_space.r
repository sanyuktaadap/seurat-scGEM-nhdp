load("datasets/GSE131907-exp-mat-by-cell-types/CRC_primarytumor_TNK_intersectfeats_Seuratobj.RData")
load("datasets/GSE131907-exp-mat-by-cell-types/TNK_sprase.RData")

library(Seurat)
library(harmony)

#####crc genes
crc_tnk_input <- CRC_TNK_Seuratobj@assays$RNA@counts
crc_gene <- rownames(crc_tnk_input)

#####lc genes
lc_tnk_input <- query_df_sparse
colnames(lc_tnk_input) <- query_df_colnames
rownames(lc_tnk_input) <- query_df_rownames
lc_gene <- query_df_rownames

#####get the intersection of gene
intersect_gene <- intersect(lc_gene, crc_gene)

#####get the gene index
crc_gene_index <- match(intersect_gene, crc_gene)
lc_gene_index <- match(intersect_gene, lc_gene)

#####prepare the expression of intersected genes
crc_sc_intersect <- crc_tnk_input[crc_gene_index, ]
lc_sc_intersect <- lc_tnk_input[lc_gene_index, ]


####the input of the seurat should be genes*cells
#####use the merged expression matrix to create the seurat object
#####add the meta data to show the cell sample origin
crc_intersect_ncells <- ncol(crc_sc_intersect)
lc_intersect_ncells <- ncol(lc_sc_intersect)
# total_intersect_ncells <- lc_intersect_ncells + crc_intersect_ncells

######samples is creating a single column with #crc cells
samples_crc <- matrix(
    c(rep("crc", crc_intersect_ncells)),
    crc_intersect_ncells,
    1)
samples_crc <- data.frame(samples_crc)
rownames(samples_crc) <- colnames(crc_sc_intersect)

######samples is creating a single column with #lc cells
samples_lc <- matrix(
    c(rep("lc", lc_intersect_ncells)),
    lc_intersect_ncells,
    1)
samples_lc <- data.frame(samples_lc)
rownames(samples_lc) <- colnames(lc_sc_intersect)


####extract the count matrix of GEM for CRC data from the trained nHDP model directly
load("model/R_nHDPresults_binary_tnk_three_layer_20230620.RData")
load("inference/lung_cancer_TNK.RData")

# Extracting value matrix
crc_value_matrix <- nHDP_trained_mb$count_matrix
lc_value_matrix <- inferred_data$count_matrix

#####you need to set the gem value for lung cancer cancer cells to be 0 (crc tnk gems)
n_dim_crc <- dim(crc_sc_intersect)[2]
n_dim_lc <- dim(lc_sc_intersect)[2]

# metadata requires rows as cells. So we are creating a matrix with n_dim rows
# dim(add_crc_tnk_gem_matrix) = 337787, 85
add_crc_tnk_gem_matrix <- matrix(0, n_dim_crc, 85)
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
add_lc_tnk_gem_matrix <- matrix(0, n_dim_lc, 85)
add_lc_tnk_gem_matrix[
    1: (dim(lc_value_matrix)[2]),
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
metadata_crc <- cbind(samples_crc, add_crc_tnk_gem_matrix)
metadata_lc <- cbind(samples_lc, add_lc_tnk_gem_matrix)

# Saving col names
crc_tnk_gem_names <- colnames(add_crc_tnk_gem_matrix)
lc_tnk_gem_names <- colnames(add_lc_tnk_gem_matrix)

rm(add_crc_tnk_gem_matrix, crc_tnk_input, add_lc_tnk_gem_matrix, CRC_TNK_Seuratobj, crc_value_matrix, lc_tnk_input, lc_value_matrix, query_df_sparse, X)
gc()

# Creating a seurat object
crc_seurat_obj <- CreateSeuratObject(
    # counts has genes * cells
    counts = crc_sc_intersect,
    project = "crc_sc_dataset",
    min.cells = 3,
    min.features = 200,
    # cells * features
    meta.data = metadata_crc)

lc_seurat_obj <- CreateSeuratObject(
    # counts has genes * cells
    counts = lc_sc_intersect,
    project = "lc_sc_dataset",
    min.cells = 3,
    min.features = 200,
    # cells * features
    meta.data = metadata_lc)


# Normalizing the data
crc_seurat_obj <- NormalizeData(
    crc_seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000)
lc_seurat_obj <- NormalizeData(
    lc_seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000)


# Identifying highly variable features
crc_seurat_obj <- FindVariableFeatures(
    crc_seurat_obj,
    selection.method = "vst",
    nfeatures = 2000)

lc_seurat_obj <- FindVariableFeatures(
    lc_seurat_obj,
    selection.method = "vst",
    nfeatures = 2000)


# Identify the 10 most highly variable genes
top10_crc <- head(VariableFeatures(crc_seurat_obj), 10)
top10_lc <- head(VariableFeatures(lc_seurat_obj), 10)

# Scaling the data
all.genes <- rownames(crc_seurat_obj)
crc_seurat_obj <- ScaleData(crc_seurat_obj, features = all.genes)
all.genes <- rownames(lc_seurat_obj)
lc_seurat_obj <- ScaleData(lc_seurat_obj, features = all.genes)

# Run PCA
crc_seurat_obj <- RunPCA(
    crc_seurat_obj,
    features = VariableFeatures(object = crc_seurat_obj))
lc_seurat_obj <- RunPCA(
    lc_seurat_obj,
    features = VariableFeatures(object = lc_seurat_obj))


# Cluster the cells
# crc_seurat_obj <- FindNeighbors(crc_seurat_obj, dims = 1:10)
# crc_seurat_obj <- FindClusters(crc_seurat_obj, resolution = 0.5)
# lc_seurat_obj <- FindNeighbors(lc_seurat_obj, dims = 1:10)
# lc_seurat_obj <- FindClusters(lc_seurat_obj, resolution = 0.5)


# Run non-linear dimension reduction
crc_seurat_obj <- RunUMAP(crc_seurat_obj, dims = 1:30, reduction = "pca")
lc_seurat_obj <- RunUMAP(lc_seurat_obj, dims = 1:30, reduction = "pca")


# Feature Plots for particular T-cell marker genes
features <- list("MKI67", "CXCL13", "PDCD1", "CTLA4")
for (feat in features){
  jpeg(
    paste("inference/results/", feat, ".jpeg", sep = ""),
    res = 300,
    width = 2000,
    height = 2000)
  print(FeaturePlot(
  object = crc_seurat_obj,
  features = feat,
  cols = c("grey85", "red"), 
  raster = FALSE))
dev.off()
}

for (feat in features){
  jpeg(
    paste("inference/results/FeaturePlots/", feat, ".jpeg", sep = ""),
    res = 300,
    width = 2000,
    height = 2000)
  print(FeaturePlot(
  object = lc_seurat_obj,
  features = feat,
  cols = c("grey85", "red"), 
  raster = FALSE))
dev.off()
}