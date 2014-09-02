source("Utils.R")

corrMat <- matrix(c(1, 0, 0.3, 0, 1, -0.5, 0.3, -0.5, 1), nrow = 3);
nSims <- 10000;
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

C_w1 <- NULL;
C_w2 <- NULL;

for(t in 2:nTime) {
  noise <- getCorrelatedVariates(nSims, corrMat); 
  R[, t] <- SimCIR(alpha_r, R[, t - 1], noise[1,], dt); 
  V[, t] <- SimCIR(alpha_v, V[, t - 1], noise[2,], dt); 
  S[,t] <- simEquity(R[, t-1], V[, t-1], S[, t - 1], noise[3,], dt); 
  
  Beta[,t] <- Beta[, t- 1] * exp(0.5 * (R[, t] + R[, t - 1]) * dt);
  
  C_w1[t-1] <- cor(noise[1,], noise[3,]);
  C_w2[t-1] <- cor(noise[2,], noise[3,]);

  C_xr[t-1] <- cor(R[,t], log(S[,t]));
  C_xv[t-1] <- cor(V[,t], log(S[,t]));
}

mean(1/t(Beta[, nTime]) * pmax(S[,nTime] - K, 0));
#mean(pmax(S[,nTime] - K, 0));

plot(C_xr, type = "l", col = "blue", ylim = c(-0.6, 0.6))
lines(C_xv, type = "l", col = "red", ylim = c(-0.6, 0.6))
abline(h = -0.5)
abline(h = 0.3)

plot(C_w1, type = "l", col = "blue", ylim = c(-0.6, 0.6))
lines(C_w2, type = "l", col = "red", ylim = c(-0.6, 0.6))
abline(h = -0.5)
abline(h = 0.3)

test <- sim.CIR.exact(alpha_r, 1000, 1, dt) 
plot(apply(R,2,mean), type = "l", col = "red")
lines(apply(test,2,mean), type = "l", col = "blue")

plot(apply(Beta,2,mean), type = "l", col = "red")
plot(apply(S,2,mean), type = "l", col = "red")
