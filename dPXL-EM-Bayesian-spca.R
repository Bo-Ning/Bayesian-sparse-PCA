dPXL_EM <- function(x, r, lambda1 = 0.1, max.iter = 1000, eps = 1e-3,
                    theta.int = NULL, lambda0.int = NULL, alpha = 1e-3, 
                    density = "l2", double.rotation = TRUE, sig2 = "known") {
  
  # lambda1 = 0.1
  # max.iter = 1000
  # eps = 1e-3
  # alpha = 1e-2
  # density = "l1"
  # double.rotation = T
  power <- 0.1
  
  # load packages and functions
  require(MASS)
  require(mvtnorm)
  require(Matrix)
  require(boot)
  require(rmutil)
  source("auxiliary_fns_EMVS.R")
  
  #------ Dimension of input matrices ------#
  p <- dim(x)[1]
  n <- dim(x)[2]
  
  # default value for lambda0 is p*log(p)
  if (lambda0.int == "prior" || is.null(lambda0.int) == T){
    lambda0 <- p*log(p)
  } else {
    lambda0 <- lambda0.int
  }
  
  # choose an initial value for theta.hat
  if (is.null(theta.int) == T) {
    library(sparsepca)
    spca.res <- 
      spca(t(x), k = r, alpha = alpha, center = F, scale = F, verbose = F)
    # obtain the loadings matrix
    theta.hat <- t(t(spca.res$loadings)*sqrt(spca.res$eigenvalues))
  } else {
    theta.hat <- theta.int
  }
  U.hat <- svd(theta.hat)$u 
  eigen.hat <- svd(theta.hat)$d
  
  # initialize parameter values
  sig2.inv <- 10
  sig2 <- 1/sig2.inv
  
  # priors
  sig2.inv.a <- 2 # par1 for sigma^{-2}
  sig2.inv.b <- 2 # par2 for sigma^{-2}
  kappa.a <- 1 # par1 for prior kappa
  kappa.b <- p+1 # par2 for prior kappa
  kappa <- kappa.a / (kappa.a + kappa.b)
  
  # evaluate the values for gamma0 and gamma1
  if (r == 1){
    gamma1 <- 
      sapply(1:p, FUN = function(j){dlaplace(theta.hat[j], s = 1/lambda1, log = T)})
    gamma0 <- 
      sapply(1:p, FUN = function(j){dlaplace(theta.hat[j], s = 1/lambda0, log = T)})
  } else {
    if (density == "l2") {
      gamma1 <- 
        sapply(1:p, FUN = function(j){dl2norm(theta.hat[j,], r, lambda1, log = T)})
      gamma0 <- 
        sapply(1:p, FUN = function(j){dl2norm(theta.hat[j,], r, lambda0, log = T)})
    } else if (density == "l1") {
      gamma1 <- 
        colSums(
          sapply(1:p, FUN = function(j){dlaplace(theta.hat[j,], 
                                                 s = 1/lambda1, log = T)}))
      gamma0 <- 
        colSums(
          sapply(1:p, FUN = function(j){dlaplace(theta.hat[j,],
                                                 s = 1/lambda0, log = T)}))
    }
  }
  gamma.hat <- exp(gamma1)^power/(exp(gamma0)^power + exp(gamma1)^power)
  penalty <- gamma.hat * lambda1 + (1-gamma.hat) * lambda0
  
  # ---- create empty matrices to store results ---- #
  theta.res <- array(NA, dim = c(p, r, max.iter))
  sig2.res <- rep(NA, max.iter)
  selection.res <- matrix(NA, p, max.iter)
  kappa.res <- rep(NA, max.iter)
  obj.fn.res <- matrix(NA, 3, max.iter)
  log.post.res <- rep(NA, max.iter)
  log.post.old <- -1e-8
  #pb  <- txtProgressBar(1, max.iter, style=3) # progress bar
  #cat("\nStarting EM algorithm: \n")
  # ----- MFVB algorithm ----- #
  for (iter in 1:max.iter) {
    
    #setTxtProgressBar(pb, iter) # progress bar
    
    if (iter == 1) {
      theta <- theta.hat
      U.hat.old <- U.hat
      eigen.hat.old <- eigen.hat
      gamma.hat.old <- gamma.hat
    }
    
if (is.null(lambda0.int) == F) {
    if (lambda0.int == "prior") {
      if (iter == 1) {
        beta <- theta
      }
      
      a.term.part1 <- sapply(1:p, FUN = function(j) {
        if (r == 1) {
          abs(beta[j])
        } else {
          if (density == "l1") {
            sum(abs(beta[j, ]))
          } else {
            sqrt(beta[j, ]^2)
          }
        }
      })
      a.term <- sum((1-gamma.hat)*a.term.part1)
      if (a.term < 1e-6) {a.term <- 1e-6}
      b.term <- r*sum(1-gamma.hat)
      
      lambda0 <- (a.term + b.term) + sqrt( (a.term + b.term)^2 + 4*p^2 * a.term )/(2*a.term)
      # lambda0 <- p/sqrt(a.term)
      # lambda0 <- 1/kappa
      print(lambda0/p)
    }
    }   
    
    ## ------------- E step ---------------- ##
    var.w <- sig2*solve(crossprod(theta) + sig2*diag(r)) # var of w
    
    # mean of w
    mean.w <- sapply(1:n, FUN = function(i) { 
      x[,i] %*% (theta %*% var.w) / sig2
    })
    
    # MSE of w
    MSE.w.i <- lapply(1:n, FUN = function(i) {
      if (r == 1) {
        mean.w[i]^2 + var.w
      } else {
        tcrossprod(mean.w[, i]) + var.w
      }
    })
    MSE.w <- Reduce("+", MSE.w.i) # sum of n second moments of w_i
    
    ## ------------- M step ---------------- ##    
    # -------------- update B -------------- #
    # coordinate update for group lasso
    library(gglasso)
    chol.MSE.w <- t(chol(MSE.w))
    
    beta.transpose <- sapply(1:p, FUN = function(j) {
      interm.term.1 <- chol.MSE.w/sqrt(2*sig2)
      interm.term.2 <- mean.w %*% x[j, ]/(2*sig2)
      x.tilde <- interm.term.1
      y.tilde <- solve(interm.term.1, interm.term.2) 
      if (density == "l2") {
        obj <- gglasso(x=x.tilde, y=c(y.tilde), 
                       group=rep(1, r), loss="ls", nlambda = 1,
                       lambda = penalty[j], intercept = F)
        return(obj$beta)
      } else if (density == "l1") {
        obj <- gglasso(x=x.tilde, y=c(y.tilde), 
                       group=c(1:r), loss="ls", nlambda = 1,
                       lambda = penalty[j], intercept = F)
        return(obj$beta)
      }
    })
    if (r == 1) {
      beta <- as.matrix(beta.transpose)
    } else {
      beta <- t(beta.transpose)
    }
    
    # -------------- update gamma -------------- # 
    if (r == 1) {
      gamma1 <- sapply(1:p, 
                       FUN = function(j){dlaplace(beta[j], s = 1/lambda1, log = T)})
      gamma0 <- sapply(1:p, 
                       FUN = function(j){dlaplace(beta[j], s = 1/lambda0, log = T)})
    } else {
      if (density == "l2") {
        gamma1 <- 
          sapply(1:p, FUN = function(j){dl2norm(beta[j,], r, lambda1, log = T)})
        gamma0 <- 
          sapply(1:p, FUN = function(j){dl2norm(beta[j,], r, lambda0, log = T)})
      } else if (density == "l1") {
        gamma1 <- 
          colSums(
            sapply(1:p, FUN = function(j){dlaplace(beta[j,], s = 1/lambda1, log = T)}))
        gamma0 <- 
          colSums(
            sapply(1:p, FUN = function(j){dlaplace(beta[j,], s = 1/lambda0, log = T)}))
      }
    }
    gamma.hat <- exp(gamma1)^power / (exp(gamma1)^power + exp(gamma0)^power)
    penalty <- gamma.hat * lambda1 + (1-gamma.hat) * lambda0
    
    # -------------- update theta -------------- #
    kappa.par1 <- kappa.a + sum(gamma.hat) - 1
    kappa <- kappa.par1 / (kappa.a + kappa.b + p - 2)
    
    # -------------- update sig2.inv -------------- #
    if (sig2 == "Unknown") {
      sig2.inv.term1 <- sum(sapply(1:n, FUN = function(i) {
        if (r == 1) {
          -2*x[, i] %*% beta * mean.w[i]
        } else {
          -2*x[, i] %*% beta %*% mean.w[, i]
        }
      }))
      sig2.inv.term2 <- sum(sapply(1:p, FUN = function(j) {
        if (r == 1) {
          beta[j]^2 * (sum(mean.w^2) + n*var.w)
        } else {
          beta[j, ] %*% (tcrossprod(mean.w) + n*var.w) %*% beta[j, ]
        }
      }))
      sig2.inv.data.term <- sum(x^2)
      log.sig2.inv.lik <- sig2.inv.data.term + sig2.inv.term1 + 
        sig2.inv.term2
      sig2.inv <- (n*p/2 + sig2.inv.a + 1)/(log.sig2.inv.lik/2 + sig2.inv.b)
      sig2 <- 1/sig2.inv
      
      # collect sig2 values
      sig2.res[iter] <- sig2
    }
    
    # if use dPXL-EM algorithm
    if (iter <= 30) {
      if (double.rotation == TRUE) {
        beta <- beta %*% chol.MSE.w/sqrt(n)
      }
    }
    
    # apply SVO to transform rotate beta back to theta
    svd.beta <- svd(beta)
    theta <- beta %*% t(svd.beta$v)
    U.hat <- svd.beta$u
    eigen.hat <- svd.beta$d
    
    # -------------- calculate objective functions ------------- #
    frob.loss <- max((tcrossprod(U.hat) - tcrossprod(U.hat.old))^2)
    eigen.loss <- max(abs(eigen.hat - eigen.hat.old))
    mis.class <- sum( abs(gamma.hat - gamma.hat.old) )
    log.post.term1 <- sapply(1:p, FUN = function(j) {
      interm.term.1 <- chol.MSE.w/sqrt(sig2)
      interm.term.2 <- mean.w %*% x[j, ]/sig2
      x.tilde <- interm.term.1
      y.tilde <- solve(interm.term.1, interm.term.2)
      diff <- sum((x.tilde %*% beta[j, ] - y.tilde)^2)
      return(diff)
    })
    log.post.term2 <- penalty *
      sapply(1:p, FUN = function(j, density) {
        if (density == "l1") {
          sum(abs(beta[j, ]))
        } else {
          sqrt(sum(beta[j,]^2))
        }
      }, density = density)
    log.post.term3 <- (sum(gamma.hat) + kappa.a - 1) * log(kappa) + 
      (p - sum(gamma.hat) + kappa.b - 1) * log(1-kappa)
    log.post <- -sum(log.post.term1/2) - sum(log.post.term2) + log.post.term3
    
    
    # collect results
    selection.res[, iter] <- gamma.hat
    theta.res[, , iter] <- theta
    kappa.res[iter] <- kappa
    obj.fn.res[, iter] <- c(frob.loss, eigen.loss, mis.class)
    log.post.res[iter] <- log.post
    
    #if (max(abs(frob.loss), abs(eigen.loss), abs(mis.class)) < eps) {
    if (abs(log.post - log.post.old) < eps) {
      if (log.post >= log.post.old) {
        break
      }
    } else {
      log.post.old <- log.post 
      U.hat.old <- U.hat
      eigen.hat.old <- eigen.hat
      gamma.hat.old <- gamma.hat
    }
  }
  
  return(list(iter = iter, selection.vec = selection.res[, 1:iter], 
              theta.vec = theta.res[, , 1:iter], 
              obj.fn.vec = obj.fn.res[, 1:iter],
              log.post.vec = log.post.res[1:iter]))
}
