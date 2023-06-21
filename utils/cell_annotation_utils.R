library(dplyr)
library(tibble)

# Helper functions
get_file_path <- function(prefix, filename, ext = "rds") {
  dir.create(prefix)
  path <- sprintf("%s/%s.%s", prefix, filename, "rds")
  
  return(path)
}

save_exp_data_by_type <- function(cell_exp_data, cell_type_annot, filepath, filename) {
  # Transpose 
  cell_exp_data <- as.data.frame(t(cell_exp_data))
  
  # Join expression matrix and cell type annotations
  # by matching rownames in exp matrix to index in cell type annotations
  # This will add a column 'rowname' with cell names to exp matrix
  cell_exp_annot_data <- left_join(rownames_to_column(cell_exp_data), cell_type_annot, by = c("rowname" = "Index"))
  rownames(cell_exp_annot_data) <- cell_exp_annot_data[, 1]
  # Excluding the new column to only have columns:
  # 1. cell names, 2. GEMs, 3. cell type annotations
  cell_exp_annot_data <- subset(cell_exp_annot_data, select = -c(rowname))
  
  # Binarize all columns except cell type (last column).
  # Keep brackets or it will be evaluated as (1:n)-1.
  cell_exp_annot_data[, 1: (ncol(cell_exp_annot_data) - 1)][
    cell_exp_annot_data[, 1: (ncol(cell_exp_annot_data) - 1)] > 0] <- 1
  
  # splitting dataset based on cell types
  split_data <- split(
    cell_exp_annot_data,
    cell_exp_annot_data$Cell_type)
  
  
  # saving each cell type separately
  for (i in seq_along(split_data)) {
    path_i <- sprintf("%s/%s", filepath, names(split_data)[i])
    split_data_i <- split_data[i][[1]]
    split_data_i <- subset(split_data_i, select = -c(Cell_type))
    
    saveRDS(
      as.data.frame(t(split_data_i)), 
      get_file_path(path_i, filename))
  }  
}
