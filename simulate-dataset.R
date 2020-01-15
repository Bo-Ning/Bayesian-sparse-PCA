# load pacakges
library(MASS)
library(pracma)
####################
# output function
####################
output <- function(U, s, s.star, p, U.hat, selection.hat) {
  frob.loss <- norm(tcrossprod(U) - tcrossprod(U.hat), type = "F")
  spectra.loss <- norm(tcrossprod(U) - tcrossprod(U.hat), type = "2")
  # FDR and FNR rate
  num1 <- sum(selection.hat[(s+1):p] == 1)
  num2 <- sum(selection.hat[1:s] == 1)
  num3 <- sum(selection.hat[(s+1):p] == 0)
  num4 <- sum(selection.hat[1:s] == 0)
  FDR <- num1 / (num1+num2)
  FNR <- num4 / (num3+num4)
  if ((num1 + num4) == 0) {
    correct <- 1
  } else {
    correct <- 0
  }
  HAM <- sum(abs(selection.hat - s.star))/p
  final.result <- c(frob.loss, spectra.loss, FDR, FNR, correct, HAM)
  names(final.result) <- c("frob", "spectra", "FDR", "FNR", "correct", "HAM")
  return(final.result)
}
##################################
### Function simulation study ####
##################################
sim.study.bayesian.spca <- function(n = 100, p = 1000, s = s, r = r, lambda0 = NULL, 
                                    alpha, density = "l2") {
  
  ############ Simulate dataset ################
  # generate orthogonal matrix
  U.s <- randortho(s, type = c("orthonormal"))
  if (r == 1) {
    U <- rep(0, p)
    U[1:s] <- as.vector(U.s[, 1:r])
  } else {
    U <- matrix(0, p, r)
    U[1:s, ] <- U.s[, 1:r]
  }
  s.star <- rep(0, p)
  s.star[1:s] <- 1
  eigenvalue <- seq(20, 10, length.out = r)
  sig2 <- 0.1
  # generate Sigma
  if (r == 1) {
    theta.true <- U * sqrt(eigenvalue)
    Sigma <- crossprod(theta.true) + sig2
  } else {
    theta.true <- U %*% sqrt(diag(eigenvalue))
    Sigma <- tcrossprod(theta.true) + sig2 * diag(p)
  }
  
  # generate n*p dataset
  X <- t(mvrnorm(n, mu = rep(0, p), Sigma = Sigma))
  
  ##############################
  # run the dPXL-EM algorithm
  source("dPXL-EM-Bayesian-spca.R")
  res.dPXL <- dPXL_EM(x = X, r = r, lambda0.int = lambda0, alpha = alpha,
                      density = density, double.rotation = T)
  iter.dPXL <- res.dPXL$iter
  selection.mtx.dPXL <- res.dPXL$selection.vec
  theta.mtx.dPXL <- res.dPXL$theta.vec
  selection.dPXL <- selection.mtx.dPXL[,iter.dPXL]
  selection.dPXL[selection.dPXL <= 0.5] <- 0
  selection.dPXL[selection.dPXL > 0.5] <- 1
  theta.dPXL <- theta.mtx.dPXL[,,iter.dPXL]
  U.hat.dPXL <- svd(theta.dPXL)$u
  output.sim.dPXL <- output(U, s,s.star, p,  U.hat.dPXL, selection.dPXL)
  return(list(output.sim.dPXL = c(output.sim.dPXL, iter.dPXL)))
}

####################################
####### select.initial.value #######
####################################
select.initial.value <- function(n = 100, p = 1000, s = s, r = r, lambda0 = NULL, 
                                 density = "l2") {
  
  library(MASS)
  library(pracma)
  
  ############ Simulate dataset ################
  # generate orthogonal matrix
  U.s <- randortho(s, type = c("orthonormal"))
  if (r == 1) {
    U <- rep(0, p)
    U[1:s] <- as.vector(U.s[, 1:r])
  } else {
    U <- matrix(0, p, r)
    U[1:s, ] <- U.s[, 1:r]
  }
  s.star <- rep(0, p)
  s.star[1:s] <- 1
  eigenvalue <- seq(20, 10, length.out = r)
  sig2 <- 0.1
  # generate Sigma
  if (r == 1) {
    theta.true <- U * sqrt(eigenvalue)
    Sigma <- crossprod(theta.true) + sig2
  } else {
    theta.true <- U %*% sqrt(diag(eigenvalue))
    Sigma <- tcrossprod(theta.true) + sig2 * diag(p)
  }
  
  # generate n*p dataset
  X <- t(mvrnorm(n, mu = rep(0, p), Sigma = Sigma))
  
  # initial alpha value
  alpha.vec <- 10^(-seq(1,10,length.out = 100))
  log.post.res <- rep(0, 100)
  source("dPXL-EM-Bayesian-spca.R")
  for (i in 1:100) {
    res.dPXL <- dPXL_EM(x = X, r = r, lambda0.int = lambda0, alpha = alpha.vec[i],
                        density = density, double.rotation = T)
    iter <- res.dPXL$iter
    log.post.res[i] <- res.dPXL$log.post.vec[iter]
  }
  index <- which(log.post.res == max(log.post.res))
  alpha.select <- alpha.vec[index]
  return(alpha.select)
}

