# By Bo Ning
# Feb, 2020

## sPCA.PXL.CAVI: Bayesian sparse pca using PXL-CAVI method
#' @param x: n-by-p data scaled and centered
#' @param r: user input rank
#' @param lambda: parameter of the variance of the slab part
#' @param max.iter: maximum iteration for the method
#' @param eps: threshold for deciding the convergence of this method
#' @param theta.int: initial value of theta
#' @param sig2: the variance of the method, default is unknown

BSPCA_jointlyRowSparse <- function(x, r, lambda = 0.01, max.iter = 1000, eps = 0.001,
                          theta.int = NULL, sig2 = "unknown", rho = 0.5, alpha = 1e-9) {
  
  # lambda <- 0.01
  # max.iter <- 1000
  # eps <- 0.001
  # theta.int <- NULL
  # sig2 <- "unknown"
  
  #------ Dimension of input matrices ------#
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  # choose an initial value for theta.hat
  if (is.null(theta.int) == T) {
    library(sparsepca)
    spca.res <- 
      spca(x, k = r, alpha = alpha, center = F, scale = F, verbose = F)
    # obtain the loadings matrix
    theta.mean <- t(t(spca.res$loadings)*sqrt(spca.res$eigenvalues))
  } else {
    theta.mean <- theta.int
  }
  U.hat <- svd(theta.mean)$u 
  eigen.hat <- svd(theta.mean)$d
  z <- as.numeric(rowSums(abs(theta.mean)) > 1e-5)
  z.old <- z
  
  # initialize parameter values
  if (sig2 == "unknown") {
    sig2 <- 1/10
  } else {
    sig2 <- sig2
  }
  kappa.para1 <- 1 # the first parameter of the beta prior
  kappa.para2 <- p + 1 # the second parameter of the beta prior
  beta.mean <- theta.mean # initial value of the mean of the eigenvectors
  beta.var <- 1e-3*diag(r) # initial value of the variance of the eigenvectors
  
  # collect results
  theta.res <- array(NA, dim = c(p, r, max.iter))
  theta.var.res <- rep(NA, max.iter)
  selection.res <- matrix(NA, p, max.iter)
  sig2.res <- rep(NA, max.iter)
  obj.fn.res <- matrix(NA, max.iter)
  theta.mean.old <- matrix(0, p, r)
  
  # start the algorithm
  for (iter in 1:max.iter) {
    
    ########################## update w #############################
    z.beta.var <- lapply(1:p, FUN = function(j) {z[j] * beta.var})
    sum.z.beta.var <- Reduce("+", z.beta.var)
    var.w <- sig2*solve(crossprod(theta.mean) + sig2*sum.z.beta.var + sig2*diag(r)) # var of w
    
    # mean of w
    mean.w <- x %*% (theta.mean %*% var.w) / sig2
    
    # MSE of w
    MSE.w.i <- lapply(1:n, FUN = function(i) {
      if (r == 1) {
        sig2*(mean.w[i]^2 + var.w)
      } else {
        sig2*(tcrossprod(mean.w[i, ]) + var.w)
      }
    })
    MSE.w <- Reduce("+", MSE.w.i) # sum of n second moments of w_i
    
    ########################## update beta #############################
    ## Mean of beta
    if (r == 1) {
      beta.term1 <- 1 / sig2 *  (sum(mean.w^2) + n*var.w) + lambda 
      beta.term2 <- 1 / sig2 * t(mean.w %*% x)
      beta.mean <- as.matrix(c(beta.term2) / c(beta.term1))
      beta.var.inv <- lambda + (sum(mean.w^2) + n*var.w)
      beta.var <- 1/beta.var.inv
    } else {
      beta.term1 <- 1 / sig2 *  (crossprod(mean.w) + n*var.w) + lambda * diag(r)
      beta.term2 <- 1 / sig2 * t(crossprod(mean.w, x))
      beta.mean <- beta.term2 %*% solve(beta.term1)
      beta.var.inv <- lambda * diag(r) + (crossprod(mean.w) + n*var.w)
      beta.var <- solve(beta.var.inv)
    }
    
    ################ Update z ####################
    z.para <- sapply(1:p, FUN = function(j) {
      if (r == 1) {
        
        term1 <- - 1 /(sig2) * sum(crossprod(mean.w, x[, j]) * beta.mean[j, ])
        term2 <- 1 / (2*sig2) * beta.mean[j, ]^2 * (sum(mean.w^2) + n*var.w)
        term3.1 <- sig2 * beta.var * (sum(mean.w^2) + n*var.w)
        term3 <- term3.1 / (2*sig2)
        term4 <- lambda / (2*sig2) * (beta.var + (beta.mean[j])^2)
        term5 <-  - r * log (lambda) / 2 - log(beta.var)/2 -
          1/2 - log(sig2)/2 + log(kappa.para2/kappa.para1)
        
      } else {
        
        term1 <- - 1/(sig2) * sum(crossprod(mean.w, x[, j]) * beta.mean[j, ])
        term2 <- c(1/(2*sig2) * beta.mean[j, ] %*% 
                     (crossprod(mean.w) + n*var.w) %*% beta.mean[j, ] )
        term3.1 <- sig2 * beta.var %*% (crossprod(mean.w) + n*var.w)
        term3 <- sum(diag(term3.1)) / (2*sig2)
        term4 <- lambda / (2*sig2) * (sum(diag(beta.var)) + crossprod(beta.mean[j, ]))
        term5 <- - r * log (lambda) / 2 - log(det(beta.var))/2 -
          1/2 - log(sig2)/2 + log(kappa.para2/kappa.para1)
      }
      
      summation <- -(term1 + term2 + term3 + term4 + term5)
    })
    
    z <- sigmoid(z.para) 
    z[z <= 0.5] <- 0
    z[z > 0.5] <- 1
    
    beta.mean[z < rho, ] <- 0 # if z is less than 0.5, then the variable would not be selected
    
    ################# Update sig2 ####################
    sig2.term1 <- sapply(1:p, FUN = function(j) {
      
      if (r == 1) {
        term1 <- sum(x[, j]^2)
        term2 <- - 2* sum(crossprod(mean.w, x[, j]) * beta.mean[j, ]) 
        term3 <- beta.mean[j, ]^2 * (sum(mean.w^2) + n*var.w)  
        term4 <- lambda * (beta.var + beta.mean[j]^2)
        summation <- term1 + term2 + term3 + term4
      } else {
        term1 <- sum(x[, j]^2)
        term2 <- - 2* sum(crossprod(mean.w, x[, j]) * beta.mean[j, ]) 
        term3 <- beta.mean[j, ] %*% (crossprod(mean.w) + n*var.w) %*% beta.mean[j, ] 
        term4 <- lambda * (sum(diag(beta.var)) + crossprod(beta.mean[j, ]))
        summation <- term1 + term2 + term3 + term4
      }
    })
    
    sig2 <- (sum(sig2.term1) + 4)/ (n*p - 4 + p*r)
    
    ##### rotate back ########
    chol.MSE.w <- t(chol(MSE.w))
    beta.mean <- beta.mean %*% (chol.MSE.w/sqrt(n))
    beta.var <- t(chol.MSE.w/sqrt(n)) %*% beta.var %*% chol.MSE.w/sqrt(n)
    
    ################# From beta to theta ###############
    beta.svd <- svd(beta.mean)
    U.mean <- beta.svd$u
    if (r == 1) {
      theta.mean <- beta.svd$u %*% beta.svd$d
    } else {
      theta.mean <- beta.svd$u %*% diag(beta.svd$d)
    }
    theta.var <- beta.var
    
    ################### Define objective function #########################
    obj.fn <- max(abs(tcrossprod(theta.mean) - tcrossprod(theta.mean.old)))
    
    ##################### Collect result ####################
    selection.res[, iter] <- z
    theta.res[, , iter] <- theta.mean
    sig2.res[iter] <- sig2
    obj.fn.res[iter] <- obj.fn
    
    ############### Stop the algorithm or go to the next iteration ########
    if (obj.fn < eps) {
      break
    } else {
      theta.mean.old <- theta.mean
      z.old <- z
    }
  }
  
  ############ Return value #################
  return(list(iter = iter, selection.vec = selection.res[, 1:iter], 
              theta.vec = theta.res[,,1:iter], sig2.vec = sig2.res[1:iter],
              obj.fn.vec = obj.fn.res[1:iter]))
}
