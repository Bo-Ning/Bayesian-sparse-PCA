BSPCA_regular <- function(x, r, lambda = 0.01, max.iter = 1000, eps = 0.001,
                          theta.int = NULL, sig2 = "unknown", rho = 0.5,
                          alpha = 1e-9) {

  
  # load packages
  library(sigmoid)
  
  #------ Dimension of input matrices ------#
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  # choose an initial value for theta.hat
  if (is.null(theta.int) == T) {
    library(sparsepca)
    spca.res <- 
      spca(x, k = r, alpha = alpha, center = F, scale = F, verbose = F)
    # obtain the loadings matrix
    theta.int <- t(t(spca.res$loadings)*sqrt(spca.res$eigenvalues))
  }
  theta.mean <- theta.int
  U.hat <- svd(theta.mean)$u # estimated value for the loadings matrix
  eigen.hat <- svd(theta.mean)$d # estimated values of eigenvalues
  z <- matrix(as.numeric(abs(theta.mean) > 1e-5), p, r) # obtain number of non-zero coordinates
  
  # initialize parameter values
  if (sig2 == "unknown") {
    sig2 <- 1e-5
  } else {
    sig2 <- sig2
  }
  kappa.para1 <- 1 # the first parameter of the beta prior
  kappa.para2 <- p + 1 # the second parameter of the beta prior
  beta.mean <- theta.mean # initial value of the mean of the eigenvectors
  beta.var <- 1e-3*diag(r) # initial value of the variance of the eigenvectors
  
  # create a matrix to collect results
  theta.res <- array(NA, dim = c(p, r, max.iter))
  theta.var.res <- rep(NA, max.iter)
  selection.res <- array(NA, dim = c(p, r, max.iter))
  sig2.res <- rep(NA, max.iter)
  obj.fn.res <- matrix(NA, max.iter)
  # theta.mean.old <- matrix(0, p, r)
  ELBO.res <- rep(NA, max.iter)
  exp.ELBO.old <- -1e10
  z.old <- z
  obj.fn <- -1e8
  
  # start the algorithm
  for (iter in 1:max.iter) {
    
    ########################## update w #############################
    # update the variance of w
    var.w.term1 <- lapply(1:p, FUN = function(j) {
      if (r == 1) {
        tcrossprod(theta.mean[j,]) + sig2 *z[j, ] * beta.var
      } else {
        tcrossprod(theta.mean[j,]) + sig2 * diag(z[j, ]) * beta.var
      }
    })
    var.w <- sig2 * solve(Reduce("+", var.w.term1) + sig2 * diag(r))
    
    # update the mean of w
    mean.w <- x %*% (theta.mean %*% var.w) / sig2
    
    # calculate the MSE of w, it is useful if parameter expansion is applied twice
    MSE.w.i <- lapply(1:n, FUN = function(i) {
      if (r == 1) {
        mean.w[i]^2 + var.w
      } else {
        tcrossprod(mean.w[i, ]) + var.w
      }
    })
    MSE.w <- Reduce("+", MSE.w.i) # sum of n second moments of w_i
    
    ########################## update beta #############################
    ## update the mean of beta and the variance of beta
    if (r == 1) {
      beta.mean.term1 <- sum(mean.w^2) + n*var.w + lambda 
      beta.mean.term2 <- crossprod(x, mean.w)
      beta.mean <- beta.mean.term2/c(beta.mean.term1)
      beta.var <- 1/beta.mean.term1
    } else {
      beta.mean.term1 <- crossprod(mean.w) + n*var.w + lambda * diag(r)
      beta.mean.term2 <- t(crossprod(mean.w, x))
      beta.mean <- beta.mean.term2 %*% solve(beta.mean.term1)
      beta.var <- diag(1/diag(beta.mean.term1))
    }
    
    ################ Update z ####################
    z.para <- sapply(1:p, FUN = function(j) {
      if (r == 1) {
        
        term1 <- - 1 /(sig2) * sum(crossprod(mean.w, x[, j]) * beta.mean[j, ])
        term2 <- 1 / (2*sig2) * beta.mean[j, ]^2 * (sum(mean.w^2) + n*var.w)
        term3 <- beta.var * (sum(mean.w^2) + n*var.w)/2
        term4 <- lambda / (2*sig2) * (sig2*beta.var + (beta.mean[j])^2)
        term5 <-  - r * log (lambda) / 2 - log(beta.var)/2 -
          1/2 - r*log(sig2)/2 + log((kappa.para2+kappa.para1)/kappa.para1)
        
      } else {
        
        term1 <- - 1/(sig2) * crossprod(mean.w, x[, j]) * beta.mean[j, ]
        term2.1 <- crossprod(mean.w) + n*var.w
        term2.2 <- sapply(1:r, FUN = function(k) {
          2*sum(term2.1[k, -k] * beta.mean[j, -k] ) *beta.mean[j, k]
        })
        term2 <- 1/(2*sig2) * (beta.mean[j, ]^2 * diag(term2.1) + term2.2)
        term3 <- diag(beta.var) * diag(term2.1) / 2
        term4 <- lambda / (2*sig2) * (sig2*diag(beta.var) + beta.mean[j, ]^2)
        term5 <- - log (lambda) / 2 - log(diag(beta.var))/2 - 1/2 - 
          log(sig2)/2 + log((kappa.para2+kappa.para1)/kappa.para1)
      }
      summation <- -(c(term1) + term2 + term3 + term4 + term5)
    })
    
    if (r == 1) {
      z <- matrix(sigmoid(z.para), p, 1)
    } else {
      z <- t(sigmoid(z.para))
    }
    
    for (k in 1:r) {
      beta.mean[z[, k] < rho, k] <- 0
    } # if z is less than 0.5, then this variable would not be selected
    
    ################# Update sig2 ####################
    sig2.term1 <- sapply(1:p, FUN = function(j) {
      
      if (r == 1) {
        term1 <- sum(x[, j]^2)
        term2 <- - 2* sum(crossprod(mean.w, x[, j]) * beta.mean[j]) 
        term3 <- beta.mean[j]^2 * (sum(mean.w^2) + n*var.w)  
        term4 <- lambda * beta.mean[j]^2
        summation <- term1 + term2 + term3 + term4
      } else {
        term1 <- sum(x[, j]^2)
        term2 <- - 2* sum(crossprod(mean.w, x[, j]) * beta.mean[j, ]) 
        term3 <- beta.mean[j, ] %*% (crossprod(mean.w) + n*var.w) %*% beta.mean[j, ] 
        term4 <- lambda * crossprod(beta.mean[j, ])
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
    
    ################# calculate ELBO ###############
    exp.ELBO.term1 <- sapply(1:p, FUN = function(j) {
      
      if (r == 1) {
        
        term1 <- n/2 * log(2*pi*sig2)
        term2 <- sum(x[, j]^2) / (2 * sig2)
        term3 <- sum(crossprod(mean.w, x[, j]) * beta.mean[j]) / sig2
        term4 <- beta.mean[j]^2 * (sum(mean.w^2) + n*var.w) / (2*sig2)
        term5 <- sum(mean.w^2) * beta.var / 2
        term6 <- var.w * beta.var / 2  
        term7 <- r * log(2*pi*sig2/lambda) / 2
        term8 <- lambda * (sig2*beta.var + beta.mean[j]^2) / (2*sig2)
        term9 <- 1/2 + r/2*log(2*pi*sig2) + z[j] * log(beta.var)
        if (sig2 == "unknown") {
          term10 <- -2*log(sig2) - 2/sig2
          summation <- - term1 - term2 - term3 -
            z[j]*(term4 + term5 + term6 + term7 + term8) + term9 + term10
        } else {
          summation <- - term1 - term2 - term3 - 
            z[j]*(term4 + term5 + term6 + term7 + term8) + term9
        }
        
      } else {
        
        term1 <- n/2 * log(2*pi*sig2)
        term2 <- sum(x[, j]^2) / (2 * sig2)
        term3 <- sum(crossprod(mean.w, x[, j]) * beta.mean[j, ]) / sig2
        term4 <- beta.mean[j, ] %*% (crossprod(mean.w) + n*var.w) %*% beta.mean[j, ] / (2*sig2)
        term5 <- sum(diag(mean.w %*% beta.var %*% t(mean.w))) / 2
        term6 <- sum(diag(var.w %*% beta.var)) / 2  
        term7 <- r * log(2*pi*sig2/lambda) / 2
        term8 <- lambda * (sig2*sum(diag(beta.var)) + crossprod(beta.mean[j, ])) / (2*sig2)
        term9 <- 1/2 + r/2*log(2*pi*sig2) + sum(z[j, ] * log(diag(beta.var)))
        if (sig2 == "unknown") {
          term10 <- -2*log(sig2) - 2/sig2
          summation <- - term1 - term2 - term3 -
            z[j]*(term4 + term5 + term6 + term7 + term8) + term9 + term10
        } else {
          summation <- - term1 - term2 - term3 - 
            z[j]*(term4 + term5 + term6 + term7 + term8) + term9
        }
      }
    }) 
    
    exp.ELBO.term2 <- - sum(rowSums(mean.w^2))/2 - n*sum(diag(var.w)) / 2
    exp.ELBO <- sum(exp.ELBO.term1) + exp.ELBO.term2
    
    ################### objective function #########################
    obj.fn <- abs(exp.ELBO - exp.ELBO.old) # the expected value of ELBO
    print(obj.fn)
    
    ##################### Collect result ####################
    selection.res[, , iter] <- z
    theta.res[, , iter] <- theta.mean
    sig2.res[iter] <- sig2
    obj.fn.res[iter] <- obj.fn
    ELBO.res[iter] <- exp.ELBO
    
    ############### Stop the algorithm or go to the next iteration ########
    if (abs(obj.fn) < eps) {
      break
    } else {
      theta.mean.old <- theta.mean
      z.old <- z
      exp.ELBO.old <- exp.ELBO
    }
  }
  
  ############ Return value #################
  return(list(iter = iter, selection.vec = selection.res[, ,1:iter], 
              theta.vec = theta.res[,,1:iter], sig2.vec = sig2.res[1:iter],
              obj.fn.vec = obj.fn.res[1:iter], ELBO.vec = ELBO.res[1:iter]))
}
