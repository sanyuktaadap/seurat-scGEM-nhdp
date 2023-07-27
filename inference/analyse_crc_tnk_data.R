load("inference/R_nHDPresults_binary_tnk_three_layer_20230620.RData")

library(msigdbr)

# GEMs v/s Cells
value_matrix <- nHDP_trained_mb$count_matrix

colnames(value_matrix) <- colnames(X)
rownames(value_matrix) <- paste("tnk_", seq(1, dim(nHDP_trained_mb$centroids)[2], 1), sep="")

# Create and Save top50 genes for each GEM
gene_name <- nHDP_trained_mb$gene
num_gem <- dim(nHDP_trained_mb$centroids)[2]
top_50_gene_name_gem <- matrix(0, num_gem, 50)
rownames(top_50_gene_name_gem) <- paste("tnk_", seq(1, num_gem, 1), sep="")

for (i in 1:num_gem) {
  # in each column of gene vs gem matrix
  temp_value <- nHDP_trained_mb$centroids[,i]
  # order the values in decreasing order
  temp_index <- order(temp_value, decreasing=TRUE)
  # extract the gene names based on the indices [1:50]
  temp_gene_name <- gene_name[temp_index[1:50]]
  # assign values to the i-th row on top_50_gene_name_gem matrix
  top_50_gene_name_gem[i,] <- temp_gene_name
}

write.csv(top_50_gene_name_gem,"inference/nhdp_3_layer_tnk_gem_top_50_genes.csv")