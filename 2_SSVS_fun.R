
# SSVS Prior
# Author: Lukas Vashold

## inputs for function:
  # X: explanatory variables in a matrix, one variable in one column
  # Y: dependent variable
  # cons: amount of constant terms (including FEs), will not be sampled in
  #       SSVS step
  # n_burn, n_draw, n_thin: number of burned/drawn iterations and thining amount
  # tau0: prior for SSVS search; will be scaled using OLS SD estimates of explanatory variables
  # tau1: prior for SSVS search; will be scaled using OLS SD estimates of explanatory variables
  # s0, S0: Hyperparameters for Gamma-distribution for sigma2
  # cons_var: prior variance on constant terms (and fixed effects) if cons_slct = FALSE
  # delta.init: a vector indicating which of the explanatory variables to
  #             include initially; default is a vector of 1s, i.e. all
  #             explanatory variables will be included
  # p.prior: a vector containing the prior inclusion probabilities; default
  #          is a vector with p.prior=1/2 for all explanatory variables
  # verbose: print ouputs

SSVSPrior <- function(Y, X, cons = 0,
                     ndraw = 5000, nburn = 2500, nthin = 1,
                     tau0 = 1e-4, tau1 = 10,
                     s0 = 0.01, S0 = 0.01,
                     cons_var = 10 ^ 8,
                     delta.init = NULL, p.prior = NULL,
                     verbose = TRUE){

  # additional inputs
  nsave <- (ndraw - nburn) / nthin
  N <- nrow(Y)
  K <- ncol(X)
  K_wo_cons <- K - cons
  # OLS quantities
  XX <- crossprod(X)
  XY <- crossprod(X, Y)
  A.OLS <- solve(XX) %*% XY
  SSE <- crossprod(Y - X %*% A.OLS)
  SIG.OLS <- SSE / (N - K)
  sigma2.draw <- as.numeric(SIG.OLS)

  # SSVS prior now with scaling with the standard error of the OLS estimator
  var.beta <- sigma2.draw * solve(XX)
  tau0 <- as.numeric(tau0 * sqrt(diag(var.beta)[1:K_wo_cons]))
  tau1 <- as.numeric(tau1 * sqrt(diag(var.beta)[1:K_wo_cons]))

  # prior and initial values
  if(is.null(delta.init)){
    delta.draw <- matrix(1, K_wo_cons, 1)
  } else {delta.draw <- delta.init}
  V.prior <- diag(c(as.numeric(delta.draw * tau1 + (1 - delta.draw) * tau0),
                    rep(cons_var, cons)))
  if(is.null(p.prior)){
    p.prior <- matrix(0.5, K_wo_cons, 1)
  }

  # store stuff
  BETA.store <- matrix(NA, nsave, K)
  SIGMA.store <- matrix(NA, nsave, 1)
  delta.store <- matrix(NA, nsave, K_wo_cons)

  if(verbose) {pb <- txtProgressBar(min = 0, max = ndraw, style = 3)}
  for (irep in  (1 - nburn):(ndraw - nburn)){
    #Draw BETA given rest from multivariate normal
    V.post <- solve(XX * 1 / sigma2.draw + diag(1 / diag(V.prior)))
    A.post <- V.post %*% (XY * 1/sigma2.draw) # +solve(V.prior)%*%A.OLS
    A.draw <- A.post + t(chol(V.post)) %*% rnorm(K)

    #Draw indicators conditional on BETA
    for (jj in 1:K_wo_cons){
      p0 <- dnorm(A.draw[[jj]], 0, sqrt(tau0[[jj]])) # draw from zj in the notes
      p1 <- dnorm(A.draw[[jj]], 0, sqrt(tau1[[jj]])) # draw from qj in the notes
      if(is.null(p.prior)){
        p11 <- p1 / (p1 + p0)
      }else{
        p11 <- (p.prior[jj] * p1) / ((1 - p.prior[jj]) * p0 + p.prior[jj] * p1)
      }
      if (p11 > runif(1)) delta.draw[[jj]] <- 1 else delta.draw[[jj]] <- 0
    }
    #Construct prior VC matrix conditional on delta
    V.prior <- diag(c(as.numeric(delta.draw * tau1 + (1 - delta.draw) * tau0),
                      rep(cons_var, cons)))

    #Simulate sigma2 from inverse gamma
    S.post <- crossprod(Y - X %*% A.draw) / 2 + S0
    s.post <- S0 + N / 2
    sigma2.draw <- 1 / rgamma(1, s.post, S.post)

    if (irep > 0 && irep %% nthin == 0){
      BETA.store[irep / nthin, ] <- A.draw
      SIGMA.store[irep / nthin, ] <- sigma2.draw
      delta.store[irep / nthin, ] <- delta.draw
    }

    if(verbose) {setTxtProgressBar(pb, (irep + nburn))}

  }

  if(verbose) {close(pb)}

  out <- list(
    "beta" = BETA.store,
    "sigma" = SIGMA.store,
    "delta" = delta.store
  )
  return(out)
}
