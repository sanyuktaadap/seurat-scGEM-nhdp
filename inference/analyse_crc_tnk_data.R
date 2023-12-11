# Load the inferred model
load("inference/R_nHDPresults_binary_tnk_three_layer_20230620.RData")

# GEMs v/s Cells
value_matrix <- nHDP_trained_mb$count_matrix

colnames(value_matrix) <- colnames(X)
rownames(value_matrix) <- paste(
  "tnk_",
  seq(1, dim(nHDP_trained_mb$centroids)[2], 1),
  sep = "")

### Create and Save top50 genes for each GEM

# Extract gene names
gene_name <- nHDP_trained_mb$gene
# Extract total number of GEMs
num_gem <- dim(nHDP_trained_mb$centroids)[2]
# Create an empty matrix with #cols = #gems and #rows = 50
top_50_gene_name_gem <- matrix(0, num_gem, 50)
# Each row is a different TNK gem, so assign row names as tnk_1, tnk_2, tnk_3, etc.
rownames(top_50_gene_name_gem) <- paste("tnk_", seq(1, num_gem, 1), sep = "")

# Iterating over each gem
for (i in 1:num_gem) {
  # in each column of gene vs gem matrix
  temp_value <- nHDP_trained_mb$centroids[,i]
  # order the values in decreasing order
  temp_index <- order(temp_value, decreasing = TRUE)
  # extract the gene names based on the indices [1:50]
  temp_gene_name <- gene_name[temp_index[1:50]]
  # assign values to the i-th row on top_50_gene_name_gem matrix
  top_50_gene_name_gem[i, ] <- temp_gene_name
}

# Saving file into csv
write.csv(
  top_50_gene_name_gem,
  "inference/nhdp_3_layer_tnk_gem_top_50_genes.csv")

### Peform Signalling Pathway Analysis

<<<<<<< HEAD
# Import gene expression data
gene_sets <- read.csv('inference/nhdp_3_layer_tnk_gem_top_50_genes.csv', sep = ',')
gene_sets <- t(gene_sets)

# For converting gene symbols to Entrez gene IDs
library(org.Hs.eg.db)
# library(biomaRt)
# To perform pathway enrichment analysis
library(clusterProfiler)

convert_to_gene_ids <- function(gene_set) {
  gene_ids <- select(
    org.Hs.eg.db,
    keys = gene_set,
    keytype = "SYMBOL",
    column = "ENTREZID")
  return(gene_ids$ENTREZID)
}

perform_signaling_pathway_analysis <- function(gene_ids) {
  # Perform KEGG pathway analysis using enrichKEGG function
  kegg_result <- enrichKEGG(gene = gene_ids,
                            organism = "hsa",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)
  return(kegg_result)
  # return(kegg_result$Description)
}

n_row <- 2: nrow(gene_sets)
n_col <- 1: ncol(gene_sets)
spa_result <- list()
for (i in n_col){
  gene_names <- gene_sets[n_row, i]
  gene_ids <- convert_to_gene_ids(gene_names)
  pathway_result <- perform_signaling_pathway_analysis(gene_ids)
  spa_result[[i]] <- pathway_result
 }

max_length <- max(sapply(spa_result, length))

# Pad the shorter pathway results with NA values
padded_spa_result <- lapply(spa_result, function(x) {
  if (length(x) < max_length) {
    c(x, rep(NA, max_length - length(x)))
  } else {
    x
  }
})

# Convert the list to a data frame
spa_result_df <- as.data.frame(padded_spa_result)

# Write the data frame to a CSV file
write.csv(spa_result_df, file = "inference/results/spa_results.csv", row.names = FALSE)
