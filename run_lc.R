library(Seurat)
library(magrittr)
library(dplyr)

source("model_windows.R")

path <- ("datasets/GSE131907-exp-mat-by-cell-types/TNK_sprase.RData")

# Using escape character to split on period instead of regex.
path_split <- strsplit(path, split = "\\.")
ext <- tail(path_split[[1]], n=1)

if (ext == "RData") {
  srt <- miceadds::load.Rdata2(filename=path)
} else if (ext == "h5Seurat") {
  srt <- SeuratDisk::LoadH5Seurat(path)
} else {
  # throw error
  stop("Unexpected file extension")
}

set.seed(111)
# X = srt@assays$RNA@counts

# Making X binary.
# X is a sparse matrix with @x containing all non-zero values.
# Thus, can directly convert all of them to 1 to make it binary.
X_binary <- srt
X_binary@x <- rep(1, length(X_binary@x))

# umap = srt@reductions$umap@cell.embeddings

model_paras = list(b0 = 0.01, g1 = 5, g2 = 1, g3 = 1/3, g4 = 2/3)
# if you have in-house housekeeping gene list, set to FALSE
# require 'Housekeeping_GenesHuman.csv' and change path in < model.R >
nHDP_init = initTree(X_binary, 
                     num_topics = c(5, 4, 3), 
                     blacklist_genes = "^MT|^RP|^HSP|B2M|MALAT1", 
                     housekeeping_genes = TRUE)
# one-batch
nHDP_trained = minibatchInfer(X_binary, 
                              nHDP_init, 
                              model_paras, 
                              max_est = 50, 
                              subtree_size = c(1, 20),
                              batch_size = ncol(X), 
                              n_epoch = 50, 
                              learning_rate = 0.01)
# mini-batch
nHDP_trained_mb = minibatchInfer(X_binary, 
                                 nHDP_init, 
                                 model_paras, 
                                 max_est = 50, 
                                 subtree_size = c(1, 20),
                                 batch_size = 1500, 
                                 n_epoch = 50, 
                                 learning_rate = 0.01)


save(list = c("nHDP_init", "nHDP_trained", "nHDP_trained_mb", "X_binary"), 
     file = "LC_primarytumor_TNK_train.RData")