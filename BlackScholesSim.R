source("Utils.R")

nSims <- 100000;
horizon <- 1;
dt <- 1/252;
nTime <- horizon/dt;
K <- 100
S0 <- 100
r <- 0.05
sigma <- 0.2
S <- matrix(S0, nrow = nSims, ncol = nTime);
sigma2 <- sigma*sigma
for(t in 2:nTime) {
  noise <- rnorm(nSims, 0, 1);   
  S[,t] <- simEquity(r, sigma2, S[, t - 1], noise, dt);   
}

mean(exp(-r * nTime * dt) * pmax(S[,nTime] - K, 0));

GBSOption("c", S0, K, nTime*dt, r, r, sigma,title = NULL, description = NULL)
