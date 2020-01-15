setwd("~/Bayesian-sparse-pca")

rm(list = ls())

### Conduct simulation study 100 per each case ###
s <- 60 # number of nonzero rows
r <- 2 # rank
p <- 1000 # number of features
n <- 200 # number of observations
density <- "l2" # specify which norm wants to be used
source("sim.study.bayesian.spca")
res <- sim.study.bayesian.spca(n, p, s = s, r = r, 
                               density = density, lambda0 = "prior")