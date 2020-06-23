
# Dirichlet-Laplace prior
# Author: Lukas Vashold

# Input:
#   - coef: estimated coefficients from the regression as a vector
#   - a.dl: prior parameter on intensity parameter of Dirichlet distribution
#   - phi.dl: posterior value of phi of DL prior from last iteration
#   - psi.dl: posterior value of psi of DL prior from last iteration
#   - cons: Number of constant terms included (fixed effects here as well?)
#   - sample_a: List. Contains stuff needed for sampling a.dl through
#         a RWMH step:
#         - scale: scale for RWMH step for drawing a.dl
#         - accept: counter for RWMH step for drawing a.dl
#         - irep: current iteration, needed for RWMH step of a.dl
#         - nburn: number of burns, needed for RWMH step of a.dl
#
# Output:
#   - V_prior_inv: vector with values used as inverse of the prior variance
#                  of the vcov-matrix of the coefficients
#   - kappa.dl: posterior value of kappa of DL prior
#   - phi.dl: posterior value of phi of DL prior, needed for next iteration
#   - psi.dl: posterior value of psi of DL prior, needed for next iteration
#   - if a.dl is sampled through RWMH step:
#         - a.dl: current value of a.dl from RWMH step, needed for next iteration
#         - scale: current value of scale, needed for next iteration
#         - accept: current values of accepted draws of theta.ng, needed for
#                   next iteration
DLPrior <- function(coef, a.dl, phi.dl, psi.dl, cons, sample_a = NULL) {

  require(GIGrvg, quietly = TRUE)
  K <- length(coef) - cons

  kappa.dl <- rgig(n = 1, lambda = K * (a.dl - 1),
                   chi = 2 * sum(abs(coef[1:K]) / phi.dl) + 1e-15, psi = 1)
  if(kappa.dl < 1e-07) kappa.dl <- 1e-07
  for (jj in 1:K){
    mu.dl <- phi.dl[jj] * kappa.dl / abs(coef[jj])
    mu.dl[mu.dl < 1e-07] <- 1e-07
    psi.dl[jj] <- 1 / rgig(1, lambda = -0.5, chi = 1, psi = 1 / mu.dl^2)
  }
  psi.dl[psi.dl < 1e-07] <- 1e-07

  T_i <- sapply(coef[1:K], function(c) {
    rgig(n = 1, lambda = a.dl - 1,
         chi = 2 * abs(c) + 1e-15, psi = 1)
  })
  T_i[T_i < 1e-07] <- 1e-07
  phi.dl <- T_i / sum(T_i)

  if(!is.null(sample_a)) {
    # Sample a.dl through a simple RWMH step
    scale <- sample_a[["scale"]]
    if(is.null(scale)) scale <- 0.1
    accept <- sample_a[["accept"]]
    irep <- sample_a[["irep"]]
    nburn <- sample_a[["nburn"]]

    require(truncnorm)
    a.prop <- rtruncnorm(1, a = 1 / K, b = 1 / 2,
                         mean = a.dl, sd = scale)

    # Ratio q(a.prop | a)/q(a | a.prop)
    q.ratio <- (pnorm((1 - a.dl) / scale) - pnorm((1 / K - a.dl) / scale)) /
               (pnorm((1 - a.prop) / scale) - pnorm(1 / K - a.prop / scale))
    # Ratio pi(a.prop | rest)/pi(a | rest)
    pi.ratio <- 2 ^ (K * (a.dl - a.prop)) * (gamma(a.dl) / gamma(a.prop)) ^ K * prod(psi.dl ^ (a.prop - a.dl))
    # Acceptance probability
    alpha <- min(1, (q.ratio * pi.ratio))
    # Accept/reject algorithm
    if (runif(1) < alpha){
      a.dl <- a.prop
      accept <- accept + 1
    }
    # Scale MH proposal during the first 75% of the burn-in stage
    if (irep < (0.75 * nburn)){
      if ((accept / irep) > 0.3)  scale <- 1.01 * scale
      if ((accept / irep) < 0.15) scale <- 0.99 * scale
    }
  }

  V_prior_inv  <- 1 / (psi.dl * phi.dl^2 * kappa.dl^2)
  V_prior_inv[V_prior_inv > 1e+10] <- 1e+10

  if(!is.null(sample_a)) {
    return(list(V_prior_inv = V_prior_inv,
                kappa.dl = kappa.dl,
                psi.dl = psi.dl,
                phi.dl = phi.dl,
                a.dl = a.dl,
                scale = scale,
                accept = accept))
  } else {
    return(list(V_prior_inv = V_prior_inv,
                kappa.dl = kappa.dl,
                psi.dl = psi.dl,
                phi.dl = phi.dl))
  }

}

