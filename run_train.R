library(Seurat)
library(magrittr)

source("model_windows.R")

path <- ("CRC_primarytumor_TNK_intersectfeats_Seuratobj.RData")

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
X = srt@assays$RNA@counts

if (length(X@x[X@x == 0]) > 0) { 
  X@x = log2(X@x + 1)
} else {
  X@x = log2(X@x)
}
#change the matrix to binary matrix: change X to be binary
X_binary <- X
X_binary@x[X_binary@x > 0] <- 1

umap = srt@reductions$umap@cell.embeddings


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
     file = "CRC_primarytumor_TNK_train.RData")
