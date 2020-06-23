
# Normal-Gamma prior
# Author: Lukas Vashold

# Input:
#   - coef: estimated coefficients from the regression as a vector
#   - tau2: values of tau2 from the last draw as a vector
#   - lambda: value of lambda from last draw
#   - theta: Hyperparameter value of scale & shape for Gamma prior on tau2
#   - c.ng: Hyperparameter value of scale for Gamma prior on lambda
#   - d.ng: Hyperparameter value of shape for Gamma prior on lambda
#   - cons: Number of constant terms included (fixed effects here as well?)
#   - prior: prior values, defaults to 0 i.e. shrinkage towards 0
#   - sample_theta: List. Contains stuff needed for sampling theta through
#         a RWMH step:
#         - scale: scale for RWMH step for drawing theta
#         - accept: counter for RWMH step for drawing theta
#         - irep: current iteration, needed for RWMH step of theta
#         - nburn: number of burns, needed for RWMH step of theta

# Output:
#   - V_prior_inv: vector with values used as inverse of the prior variance
#                  of the vcov-matrix of the coefficients
#   - tau.ng: posterior values of tau.ng, needed for next iteration
#   - lambda.ng: posterior values of lambda.ng, needed for next iteration
#   - if theta.ng is sampled through RWMH step:
#         - theta.ng: current value of theta.ng from RWMH step
#         - scale: current value of scale, needed for next iteration
#         - accept: current values of accepted draws of theta.ng, needed for
#                   next iteration

NGPrior <- function(coef, tau.ng, lambda.ng, theta.ng, c.ng, d.ng, cons,
                    prior = NULL, sample_theta = NULL) {
  require(GIGrvg, quietly = TRUE)
  K <- length(coef) - cons

  if(is.null(prior)) prior <- matrix(0, K, 1)

  # draw lambda
  cl <- c.ng + K * theta.ng
  dl <- d.ng + 0.5 * theta.ng * sum(tau.ng)
  lambda.ng <- rgamma(1, cl, dl)

  # draw tau2
  for(kk in 1:K) {
    tau.ng[kk] <- GIGrvg::rgig(1,
                               lambda = theta.ng - 0.5,
                               chi = (coef[kk] - prior[kk]) ^ 2,
                               psi = lambda.ng * theta.ng)
  }
  tau.ng[tau.ng < 1e-7] <- 1e-7
  V_prior_inv <- 1 / tau.ng

  # Sample theta.ng through a simple RWMH step
  if(!is.null(sample_theta)) {
    scale <- sample_theta[["scale"]]
    if(is.null(scale)) scale <- 0.5
    accept <- sample_theta[["accept"]]
    irep <- sample_theta[["irep"]]
    nburn <- sample_theta[["nburn"]]

    theta_prop <- exp(rnorm(1, 0, scale)) * theta.ng
    post_theta_prop <- theta_post(theta = theta_prop,
                                  tau2 = as.vector(tau.ng),
                                  lambda = lambda.ng)
    post_theta_old <- theta_post(theta = theta.ng,
                                  tau2 = as.vector(tau.ng),
                                  lambda = lambda.ng)
    post.diff <- post_theta_prop - post_theta_old +
                       log(theta_prop) - log(theta.ng)
    post.diff <- ifelse(is.nan(post.diff), -Inf, post.diff)
    if (post.diff > log(runif(1, 0, 1))){
      theta.ng <- theta_prop
      accept <- accept + 1
    }
    # Scale MH proposal during the first 75% of the burn-in stage
    if (irep < (0.75 * nburn)){
      if ((accept / irep) > 0.3)  scale <- 1.01 * scale
      if ((accept / irep) < 0.15) scale <- 0.99 * scale
    }
  }

  if(!is.null(sample_theta)) {
    return(list(V_prior_inv = V_prior_inv,
                tau.ng = tau.ng,
                lambda.ng = lambda.ng,
                theta.ng = theta.ng,
                scale = scale,
                accept = accept))
  } else {
    return(list(V_prior_inv = V_prior_inv,
                tau.ng = tau.ng,
                lambda.ng = lambda.ng))
  }
}


theta_post <- function(theta = theta, lambda = lambda, tau2 = tau2,
                       k = length(tau2), rat = 1){
  logpost <- sum(dgamma(tau2, theta, (theta * lambda / 2), log=TRUE)) +
    dexp(theta, rate=rat, log=TRUE)
  return(logpost)
}

