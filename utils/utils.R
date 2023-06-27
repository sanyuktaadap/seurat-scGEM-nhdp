
get2dIndexFrom1d <- function(index, nr) {
  i <- ((index - 1) %% nr) + 1
  j <- ((index - 1) %/% nr) + 1
  
  return(list("i"=i, "j"=j))
}

get1dIndexFrom2d <- function(i, j, nr) {
  return(((j - 1) * nr) + i) 
}

# Converts to binary sparse by taking list of indices where value is
# 1. Also, takes number of rows.
convertGEMsDfToSparse <- function(df) {
  # Gets 1d indices where df is 1.
  indices_1d <- which(df > 0)
  
  if (length(indices_1d) == 0) {
    return(NULL)
  }
  
  index_2d <- get2dIndexFrom1d(indices_1d, nrow(df))
  
  return(Matrix::sparseMatrix(
    index_2d$i,
    index_2d$j,
    x=rep(1, length(indices_1d)),
    dims=c(nrow(df), ncol(df))))
}

areMatricesEqual <- function(a, b) {
  
  if (max(abs(a - b)) == 0) {
    return(TRUE)
  }
  
  return(FALSE)
}

# Sparsifies dataframe one step at a time.
bigDataFrameSparsify <- function(df, split = 1000) {
  
  sparse_df <- convertGEMsDfToSparse(df[1:split])
  
  for(i in seq(split + 1, ncol(df), split)) {
    end <- i + split - 1
    
    if (end > ncol(df)) {
      end <- ncol(df)
    }
    
    print(sprintf("Sparsifying from %d to %d", i, end))
    
    sparse_df_i <- convertGEMsDfToSparse(df[i: end])
    
    if (is.null(sparse_df)) {
      sparse_df <- sparse_df_i
    } else if (!is.null(sparse_df_i)) {
      sparse_df <- cbind(sparse_df, sparse_df_i)
    }
  }
  
  return(sparse_df)
}