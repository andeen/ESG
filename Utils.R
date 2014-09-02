
getCorrelatedVariates <- function(nSims, corrMat) { 
  n <- nrow(corrMat);
  cholDecomp <- chol(corrMat);
  noise <- matrix(rnorm(n * nSims, 0, 1), nrow = n, ncol = nSims);
  
  corrNoise <- t(cholDecomp) %*% noise;
  return(corrNoise)
}

sim.CIR.exact <- function(params, n, T, dt) {
  ## returns nx(T/dt+1) matrix containing exact simulation of n paths of the
  ## short rate process of length T sampled at points that are dt years apart
  ##
  kappa <- params[1]
  theta <- params[2]
  sigma <- params[3]
  x0 <- params[4]
  
  c <- sigma^2*(1-exp(-kappa*dt)) / (4*kappa)   # chisq multiplier
  npm <- exp(-kappa*dt)/c                       # Non-centrality multiplier
  degf <- 4*theta*kappa/sigma^2               # chisq degrees of freedom
  
  t <- seq(0, T, dt)
  m <- length(t)
  
  x <- matrix(nrow=n, ncol=m)
  x[,1] <- x0
  
  for(i in 1:(m-1)){
    ncp <- x[,i]*npm          # chisq non-centrality
    x[,i+1] <- c*rchisq(n, df=degf, ncp=ncp)
  }
  
  return(x)
}

simEquity <- function(r, sigmaEquity, Sprev, noise, dt) {
  return(Sprev * exp((r - 0.5 * sigmaEquity)*dt + sqrt(sigmaEquity*dt)*noise))
}

simLogEquity <- function(r, sigmaEquity, Xprev, noise, dt) {
  return(Xprev + (r - 0.5 * sigmaEquity)*dt + sqrt(sigmaEquity*dt)*noise)
}

SimCIR <- function(params, xprev, noise, dt) {
  
  kappa <-params[1];
  theta <-params[2];
  sigma <-params[3];
  x0 <- params[4];  
  x <- xprev + kappa * (theta - xprev) * dt + sigma * sqrt(xprev * dt) * noise;  
  return(x);
}
