corrMat <- matrix(c(1, 0, 0.3, 0, 1, -0.5, 0.3, -0.5, 1), nrow = 3);
nSims <- 100000;
horizon <- 1;
dt <- 1/252;
nTime <- horizon/dt;
K <- 100
S0 <- 100
alpha_v <- c(3, 0.4, 0.1, 0.2);
alpha_r <- c(3, 0.06, 0.1, 0.08);

S <- matrix(S0, nrow = nSims, ncol = nTime);
R <- matrix(alpha_r[4], nrow = nSims, ncol = nTime);
V <- matrix(alpha_v[4], nrow = nSims, ncol = nTime);
Beta <- matrix(1, nrow = nSims, ncol = nTime);
C_xr <- NULL;
C_xv <- NULL;
for(t in 2:nTime) {
  noise <- getCorrelatedVariates(nSims, corrMat); 
  R[, t] <- SimCIR(alpha_r, R[, t - 1], noise[1,], dt); 
  V[, t] <- SimCIR(alpha_v, V[, t - 1], noise[2,], dt); 
  S[,t] <- simEquity(R[, t], V[, t], S[, t - 1], noise[3,], dt); 
  
  Beta[,t] <- Beta[, t- 1] * exp(0.5 * (R[, t] + R[, t - 1]) * dt);
  
  C_xr[t-1] <- cor(R[,t], log(S[,t]));
  C_xv[t-1] <- cor(V[,t], log(S[,t]));
}

mean(1/Beta[, nTime] * pmax(S[,nTime] - K, 0));
#mean(pmax(S[,nTime] - K, 0));

plot(C_xr, type = "l", col = "blue", ylim = c(-0.6, 0.6))
lines(C_xv, type = "l", col = "red", ylim = c(-0.6, 0.6))
abline(h = -0.5)
abline(h = 0.3)

test <- sim.CIR.exact(alpha_r, 1000, 1, dt) 


simEquity <- function(r, sigmaEquity, Sprev, noise, dt) {
  return(Sprev * exp((r - 0.5 * sigmaEquity)*dt + sqrt(sigmaEquity*dt)*noise))
}

SimCIR <- function(params, xprev, noise, dt) {
  
  kappa <-params[1];
  theta <-params[2];
  sigma <-params[3];
  x0 <- params[4];  
  x <- xprev + kappa * (theta - xprev) * dt + sigma * sqrt(xprev * dt) * noise;  
  return(x);
}

getCorrelatedVariates <- function(nSims, corrMat) { 
  n <- nrow(corrMat);
  cholDecomp <- chol(corrMat);
  noise <- matrix(rnorm(n * nSims, 0, 1), nrow = n, ncol = nSims);
  
  corrNoise <- cholDecomp %*% noise;
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

