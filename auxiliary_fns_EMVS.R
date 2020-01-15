dl2norm <- function(b, r, lambda, log = T) {
  if (r > 1) {
    log.density <- r*log(pi)/2 + r*log(lambda) - lgamma(r+1) + lgamma(r/2 + 1) - 
      lambda*sqrt(sum(b^2))
  } else {
    log.density <- dlaplace(b, s = 1/lambda, log = log)
  }
  if (log == TRUE) {
    return(log.density)
  } else {
    return(exp(log.density))
  }
}

objective.function <- 
  function(x, B, r, mean.w, S.w, theta, sig2.inv, lambda0, lambda1,
           theta.a, theta.b, sig2.inv.a, sig2.inv.b, density = "l2") {
  
  p <- dim(x)[1]
  n <- dim(x)[2]
  
  E.log.lik.S.w.part <- sapply(1:p, FUN = function(j) {
    - B[j, ] %*% S.w %*% B[j, ] * sig2.inv / 2
  })
  E.log.lik.mean.w.part <- sapply(1:n, FUN = function(i) {
    if (r == 1) {
      sum(crossprod(x[, i] * sig2.inv, B) * mean.w[i]) 
    } else {
      sum(crossprod(x[, i] * sig2.inv, B) * mean.w[, i]) 
    }
  })
  E.log.lik <- sum(n*p*log(sig2.inv))/2 - sum(crossprod(x*sqrt(sig2.inv)))/2 +
    sum(E.log.lik.mean.w.part) - sum(E.log.lik.S.w.part)
  
  if (theta == 0) {
    theta <- 1e-8
  }
  
  if (density == "l2") {
    
    log.dl2norm.lambda0 <- 
      sapply(1:p, FUN = function(j) {dl2norm(B[j, ], r, lambda0, log=T)})
    log.dl2norm.lambda1 <- 
      sapply(1:p, FUN = function(j) {dl2norm(B[j, ], r, lambda1, log=T)})
    log.prior.B <- 
      sum(log(theta)+log.dl2norm.lambda1 + log(1-theta)+log.dl2norm.lambda0)
  } else if (density == "l1") {
    if (r == 1) {
      log.dl1norm.lambda0 <- 
        sapply(1:p, FUN = function(j) {dlaplace(B[j, ], s = 1/lambda0, log=T)})
      log.dl1norm.lambda1 <- 
        sapply(1:p, FUN = function(j) {dlaplace(B[j, ], s = 1/lambda1, log=T)})
    } else {
      log.dl1norm.lambda0 <- 
        colSums(sapply(1:p, FUN = function(j) {dlaplace(B[j, ], s = 1/lambda0, log=T)}))
      log.dl1norm.lambda1 <- 
        colSums(sapply(1:p, FUN = function(j) {dlaplace(B[j, ], s = 1/lambda1, log=T)}))
    }
      log.prior.B <- 
        sum(log(theta)+log.dl1norm.lambda1 + log(1-theta)+log.dl1norm.lambda0)
  }

  log.prior.theta <- theta.a * log(theta) + theta.b * log(1-theta)
  
  log.prior.sig2.inv <- 
    sum(dgamma(sig2.inv, shape = sig2.inv.a, rate = sig2.inv.b, log = T))
  
  obj.fn <- E.log.lik + log.prior.B + log.prior.theta + log.prior.sig2.inv
  
  return(obj.fn)
}
