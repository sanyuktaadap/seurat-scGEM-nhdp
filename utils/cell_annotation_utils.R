library(dplyr)
library(tibble)

# Helper functions
get_file_path <- function(prefix, cell_type, ext = "rds") {
  dir.create(prefix)
  path <- sprintf("%s/%s.%s", prefix, cell_type, "rds")
  
  return(path)
}

# Global Variables
ANNOTATION_PATH <- "datasets/GSE131907_Lung_Cancer_cell_annotation.txt"
CELL_EXP_PATH <- "datasets/GSE131907_100.rds"
FILENAME_PREFIX <- "datasets/gse131907"

# Load and read cell annotations
annot <- read.table(ANNOTATION_PATH, header=TRUE, sep="\t")
cell_type_annot <- annot[ , c(1, 5)]

# load and read cell expression data
cell_exp_data <- readRDS(CELL_EXP_PATH)
cell_exp_data <- cell_exp_data[1:50,1:50]

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
data <- cell_exp_annot_data[, 1: (ncol(cell_exp_annot_data) - 1)]
data[data > 0] <- 1

# splitting dataset based on cell types
split_data <- split(
  cell_exp_annot_data,
  cell_exp_annot_data$Cell_type)

# saving each cell type separately
for (i in seq_along(split_data)) {
  saveRDS(
    split_data[i], 
    get_file_path(FILENAME_PREFIX, names(split_data)[i]))
}