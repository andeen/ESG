source("Utils.R")
library(NMOF)

corrMat <- matrix(c(1,-0.5, -0.5, 1), nrow = 2);
nSims <- 10000;
horizon <- 1;
dt <- 1/252;
nTime <- horizon/dt;
K <- 100
S0 <- 100
alpha_v <- c(3, 0.4, 0.1, 0.2);
r <- 0.05

rho <- -0.5

S <- matrix(S0, nrow = nSims, ncol = nTime);
V <- matrix(alpha_v[4], nrow = nSims, ncol = nTime);

for(t in 2:nTime) {
  noise <- getCorrelatedVariates(nSims, corrMat);   
  noise1 <- rnorm(nSims, 0, 1)
  noise2 <- rho * noise1 + sqrt(1 - rho*rho)*rnorm(nSims, 0, 1)
  
  V[, t] <- SimCIR(alpha_v, V[, t - 1], noise[1,], dt); 
  S[,t] <- simEquity(r, V[, t-1], S[, t - 1], noise[2,], dt);   
  #V[, t] <- SimCIR(alpha_v, V[, t - 1], noise1, dt); 
  #S[,t] <- simEquity(r, V[, t - 1], S[, t - 1], noise2, dt);   
}

mean(exp(-r * nTime * dt) * pmax(S[,nTime] - K, 0));
#mean(pmax(S[,nTime] - K, 0));
plot(apply(V,2,mean),type = "l")

callHestoncf(S0, K, nTime * dt, r, 0, alpha_v[4], alpha_v[2], -0.5, alpha_v[1], alpha_v[3], implVol = FALSE)
