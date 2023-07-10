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

original <- X[1:50, 1:50]
#convert sparse matrix (dots to zeros)
original_matrix <- as.matrix(original)
original_matrix[original_matrix>0] <- 1

#write.csv(original_matrix, "test.csv")


comparison <- X[1:50, 1:50]
comparison@x <- rep(1, length(comparison@x))
comparison_matrix <- as.matrix(comparison)


identical(original_matrix, comparison_matrix)

#write.csv(comparison_matrix, "test2.csv")