load("inference/lung_cancer_TNK.RData")

library(msigdbr)

value_matrix <- inferred_data$count_matrix
gene_by_gem <- inferred_data$centroids

# GEMs v/s Cells
colnames(value_matrix) <- inferred_data$cell
rownames(value_matrix) <- paste("tnk_", seq(1, dim(gene_by_gem)[2], 1), sep="")

# Create and Save top50 genes for each GEM
gene_name <- inferred_data$gene
num_gem <- dim(inferred_data$centroids)[2]
top_50_gene_name_gem <- matrix(0, num_gem, 50)
rownames(top_50_gene_name_gem) <- paste("tnk_", seq(1, num_gem, 1), sep="")

for (i in 1:num_gem) {
  # in each column of gene vs gem matrix
  temp_value <- inferred_data$centroids[,i]
  # order the values in decreasing order
  temp_index <- order(temp_value, decreasing=TRUE)
  # extract the gene names based on the indices [1:50]
  temp_gene_name <- gene_name[temp_index[1:50]]
  # assign values to the i-th row on top_50_gene_name_gem matrix
  top_50_gene_name_gem[i, ] <- temp_gene_name
}


# Saving files
write.csv(top_50_gene_name_gem, "inference/lc_tnk_gem_top_50_genes.csv")

save(list = ("value_matrix"),
     file = "inference/lung_cancer_TNK_value_matrix.csv")

save(list = ("gene_by_gem"),
     file = "inference/lung_cancer_TNK_gene_by_gem.csv")