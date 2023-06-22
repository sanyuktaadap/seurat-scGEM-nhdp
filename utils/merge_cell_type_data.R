source("utils/cell_annotation_utils.R")

DATASET_PATH = "datasets/GSE131907_res"

cell_types <- list.files(DATASET_PATH)

for (i in seq(cell_types)) {
  print(sprintf("Solving %s", cell_types[i]))
  file_list <- list.files(sprintf("%s/%s", DATASET_PATH, cell_types[i]))
  
  z <- NULL

  for (j in seq(1, length(file_list) - 1)) {
    if (j == 1) {
      file_path_1 <- sprintf("%s/%s/%s", DATASET_PATH, cell_types[i], file_list[j])
      file_path_2 <- sprintf("%s/%s/%s", DATASET_PATH, cell_types[i], file_list[j+1])
      x1 <- readRDS(file_path_1)
      x2 <- readRDS(file_path_2)
      z <- cbind(x1, x2)
    }
    else {
      file_path <- sprintf("%s/%s/%s", DATASET_PATH, cell_types[i], file_list[j+1])
      x <- readRDS(file_path)
      z <- cbind(z, x)
    }
  }
  
  # Add saving logic.
  saveRDS(z, file = get_file_path(DATASET_PATH, cell_types[i]))
  # file_path_1 <- sprintf("%s/%s/%s", DATASET_PATH, cell_types[i], file_list[1])
  # z <- readRDS(file_path_1)
  # 
  # for (j in seq(2, file_list)) {
  #   file_path <- sprintf("%s/%s/%s", DATASET_PATH, cell_types[i], file_list[j])
  #   
  #   x <- readRDS(file_path_2)
  #   z <- cbind(z, x)
  # }
  
}
