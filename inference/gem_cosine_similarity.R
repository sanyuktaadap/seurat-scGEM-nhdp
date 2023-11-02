library(lsa)

load("inference/R_nHDPresults_binary_tnk_three_layer_20230620.RData")
load("inference/lung_cancer_TNK.RData")

# Preparing CRC gene by gem matrix
crc_gene_by_gem_mat <- nHDP_trained_mb$centroids
rownames(crc_gene_by_gem_mat) <- nHDP_trained_mb$gene
# Assigning colnames
colnames(crc_gene_by_gem_mat) <- paste("crc_tnk_",
                                   seq(1, dim(crc_gene_by_gem_mat)[2], 1),
                                   sep = "")

# Preparing LC gene by gem matrix
lc_gene_by_gem_mat <- inferred_data$centroids
rownames(lc_gene_by_gem_mat) <- inferred_data$gene
# Assigning colnames
colnames(lc_gene_by_gem_mat) <- paste("lc_tnk_",
                                   seq(1, dim(lc_gene_by_gem_mat)[2], 1),
                                   sep = "")

# Create an empty data frame
cosine_values <- list()

for (i in seq(1,3)) {
    crc_idx <- order(crc_gene_by_gem_mat[, i], decreasing = TRUE)
    lc_idx <- order(lc_gene_by_gem_mat[, i], decreasing = TRUE)

    crc_top_50_idx <- crc_idx
    lc_top_50_idx <- lc_idx

    crc_sorted_values <- crc_gene_by_gem_mat[crc_top_50_idx, i]
    lc_sorted_values <- lc_gene_by_gem_mat[lc_top_50_idx, i]

    crc_sorted_values <- t(as.data.frame(crc_sorted_values))
    lc_sorted_values <- t(as.data.frame(lc_sorted_values))

    crc_sorted_gene_names <- colnames(crc_sorted_values)
    lc_sorted_gene_names <- colnames(lc_sorted_values)

    comb_gene_names <- union(crc_sorted_gene_names, lc_sorted_gene_names)

    crc_df <- data.frame(matrix(data = 0, ncol = length(comb_gene_names), nrow = 1))
    colnames(crc_df) <- comb_gene_names
    lc_df <- data.frame(matrix(data = 0, ncol = length(comb_gene_names), nrow = 1))
    colnames(lc_df) <- comb_gene_names

    for (gene_name in crc_sorted_gene_names) {
        crc_df[gene_name] <- crc_sorted_values[, gene_name]
    }
    for (gene_name in lc_sorted_gene_names) {
        lc_df[gene_name] <- lc_sorted_values[, gene_name]
    }

    crc_vec <- unlist(crc_df[1,], use.names=FALSE)
    lc_vec <- unlist(lc_df[1,], use.names=FALSE)

    result <- cosine(crc_vec, lc_vec)
    cosine_values <- append(cosine_values, result)
}