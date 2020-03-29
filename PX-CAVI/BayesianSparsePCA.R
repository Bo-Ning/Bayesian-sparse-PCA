BayesianSparsePCA <- function(x, r, lambda = 0.01, max.iter = 1000, eps = 0.001,
                              theta.int = NULL, sig2 = "unknown", rho = 0.5,
                              alpha = 1e-9, sparsity = "regular") {
  
  if (sparsity == "regular") {
    
    source("BSPCA_regular.R")
    res <- BSPCA_regular(x, r, lambda, max.iter, eps, theta.int, sig2, rho, alpha)
    iter <- res$iter
    selection.vec <- res$selection.vec
    theta.vec <- res$theta.vec
    sig2.vec <- res$sig2.vec
    obj.fn.vec <- res$obj.fn.vec
    
  } else if (sparsity == "jointlyRowGroup") {
    
    source("BSPCA_jointlyRowSparse.R")
    res <- BSPCA_jointlyRowSparse(x, r, lambda, max.iter, eps, theta.int, sig2, rho, alpha)
    iter <- res$iter
    selection.vec <- res$selection.vec
    theta.vec <- res$theta.vec
    sig2.vec <- res$sig2.vec
    obj.fn.vec <- res$obj.fn.vec

  }
  
  
  return(list(iter = iter, selection.vec = selection.vec, theta.vec = theta.vec, 
              sig2.vec = sig2.vec, obj.fn.vec = obj.fn.vec))
  
}
