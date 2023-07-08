source("model_windows.R")
source("utils/utils.R")

DF_PATH = "datasets/GSE131907-exp-mat-by-cell-types/TNK.rds"
SPARSE_PATH = "datasets/GSE131907-exp-mat-by-cell-types/TNK_sprase.RData"

if (file.exists(SPARSE_PATH)) {
  load(SPARSE_PATH)
} else {
  query_df <- readRDS(DF_PATH)
  
  # Storing rownames and colnames. Storing data as sparse matrix.
  query_df_colnames <- colnames(query_df)
  query_df_rownames <- rownames(query_df)
  
  # Converting the matrix to a sparse matrix
  query_df_sparse <- bigDataFrameSparsify(query_df)
  
  save(list = c("query_df_colnames", "query_df_rownames", "query_df_sparse"),
       file = SPARSE_PATH)
  
  rm(query_df)
}


load("inference/R_nHDPresults_binary_tnk_three_layer_20230620.RData")
rm(X, nHDP_init)

model_params = list(b0 = 0.01, g1 = 5, g2 = 1, g3 = 1/3, g4 = 2/3)
inferred_data = minibatchInfer(query_df_sparse, 
                               nHDP_trained_mb, 
                               model_params, 
                               max_est = 50, 
                               subtree_size = c(1, 20),
                               batch_size = 1500, 
                               n_epoch = 1,
                               learning_rate = 0.01)