
# Simulate data
X <- matrix(rnorm(10000), ncol = 10)
beta <- c(0, 1, 1, 0.1, -1, -2, -0.05, 2, -0.5, 0)
y <- X %*% beta + rnorm(nrow(X))

# Estimate
source("1_estimation.R")

out1 <- clm(cbind(y, X))
apply(out1$beta, 2, mean)

out2 <- clm(cbind(y, X), shrinkage = "NG")
apply(out2$beta, 2, mean)

out3 <- clm(cbind(y, X), shrinkage = "DL")
apply(out3$beta, 2, mean)

out4 <- clm(cbind(y, X), shrinkage = "HS")
apply(out4$beta, 2, mean)

out5 <- clm(cbind(y, X), shrinkage = "SSVS", tau0 = 1e-4, tau1 = 10)
apply(out5$beta, 2, mean)

# Plot
source("9_plot_beta.R")
source("9_plot_beta_shrink.R")
source("9_plot_hyper.R")

plot_beta(out1, type = "density")

plot_beta_shrink(list(out1, out2), c("none", "NG"), type = "violin",
  subs = c("beta_1", "beta_2"))

plot_beta_shrink(list(out1, out2, out3, out4),
  c("none", "NG", "DL", "HS"), type = "density", subs = c("beta_2"))

plot_beta_shrink(list(out1, out2, out5),
  c("none", "NG", "SSVS"), type = "density", subs = c("beta_1"))

plot_hyper(out2, type = "boxplot")
