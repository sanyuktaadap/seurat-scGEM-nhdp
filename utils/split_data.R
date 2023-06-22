source("utils/cell_annotation_utils.R")

PATH_PREFIX = "datasets/GSE131907/batches"
CELL_EXP_PATH <- "datasets/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds"
BATCH_SIZE = 1000

# Load exp data
cell_exp_data <- readRDS(CELL_EXP_PATH)

file_index = 1

for (i in seq(1, ncol(cell_exp_data), BATCH_SIZE)) {
  print(sprintf("Solving %d batch", file_index))
  end <- i + BATCH_SIZE - 1
  
  if (end > ncol(cell_exp_data)) {
    end <- ncol(cell_exp_data)
  }
  
  saveRDS(
    cell_exp_data[, i: end], 
    file = get_file_path(PATH_PREFIX, sprintf("%d", file_index)))
  
  file_index <- file_index + 1
}

