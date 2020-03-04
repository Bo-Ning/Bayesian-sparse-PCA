output <- function(data, pcs, selection, sig2) {
  
  pcs <- as.matrix(pcs)
  rank <- dim(pcs)[2]
  
  ### Loadings matrix 
  colnames(pcs) <- paste(1:rank, "pc")
  
  ### standard deviation
  eigenvalues <- colSums(pcs^2) + sig2
  
  
  ### Sparse component
  if (rank > 1) {
    active.rows <- lapply(1:rank, FUN = function(k) {which(selection[, k] == 1)})
  } else {
    active.rows <- which(selection == 1)
  }
  
  loadings <- svd(pcs)$u
  # score functions
  scores <- data %*% loadings
  
  return(list(eigenvalues = eigenvalues, active.rows = active.rows, 
              scores = scores))
  
}