source("utils/cell_annotation_utils.R")

# Global Variables
ANNOTATION_PATH <- "datasets/GSE131907_Lung_Cancer_cell_annotation.txt"
CELL_EXP_PATH <- "datasets/GSE131907/batches"
FILEPATH <- "datasets/GSE131907_res/"

# Load and read cell annotations
annot <- read.table(ANNOTATION_PATH, header=TRUE, sep="\t")
cell_type_annot <- annot[ , c(1, 5)]
cell_exp_filenames <- list.files(CELL_EXP_PATH)

i <- 1

# load and read cell expression data
for (filename in cell_exp_filenames) {
  print("\n")
  print(sprintf("Solving %d/%d files", i,length(cell_exp_filenames)))
  
  name <- strsplit(filename, "\\.")[[1]][1]
  cell_exp_batch <- readRDS(get_file_path(CELL_EXP_PATH, name))
  save_exp_data_by_type(cell_exp_batch, cell_type_annot, FILEPATH, name)
  
  i <- i + 1
}
