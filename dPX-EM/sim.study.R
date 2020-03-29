rm(list = ls())

setwd("~/dPXL-EM/")
## simulate dataset

sim.data <- function(n, p, s, r) {
  
  library(pracma)
  library(MASS)
  
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
    Sigma <- tcrossprod(theta.true) + sig2*diag(p)
  } else {
    theta.true <- U %*% sqrt(diag(eigenvalue))
    Sigma <- tcrossprod(theta.true) + sig2 * diag(p)
  }
  
  # generate n*p dataset
  x <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  
  return(x)
}


set.seed(2019)
n <- 200
p <- 1000
s <- n * 0.1
r <- 2

source("dPXL-EM-Bayesian-spca.R")
x <- sim.data(n, p, s, r)
res <- BayesianSparsePCA(x, r)
iter <- res$iter
pcs <- res$theta.vec[,,iter]
sig2 <- res$sig2.vec[iter]
selection <- res$selection.vec[,,iter]
selection[selection > 0.5] <- 1
selection[selection <= 0.5] <- 0

source("output.R")
output <- output(x, pcs, selection, sig2)
eigenvalues <- output$eigenvalues
active.rows <- output$active.rows
scores <- output$scores
