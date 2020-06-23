
# Horseshoe prior
# Author: Lukas Vashold

# Input:
#   - coef: estimated coefficients from the regression as a vector
#   - tau.hs: posterior value of tau of HS prior from last iteration
#   - nu.hs: posterior value of nu of HS prior from last iteration
#   - zeta.hs: posterior value of zeta of HS prior from last iteration
#   - cons: Number of constant terms included (fixed effects here as well?)
#
# Output:
#   - V_prior_inv: vector with values used as inverse of the prior variance
#                  of the vcov-matrix of the coefficients
#   - lambda.hs: posterior value of lambda of HS prior, needed for next iteration
#   - tau.hs: posterior value of tau of HS prior, needed for next iteration
#   - nu.hs: posterior value of nu of HS prior, needed for next iteration
#   - zeta.hs: posterior value of zeta of HS prior, needed for next iteration
HSPrior <- function(coef, tau.hs, nu.hs, zeta.hs, cons) {

  K <- length(coef) - cons

  lambda.hs <- 1 / rgamma(K, shape = 1,
                          rate = 1 / nu.hs + coef[1:K] ^ 2 / (2 * tau.hs))
  tau.hs <- 1 / rgamma(1, shape = (K + 1) / 2,
                       rate = 1 / zeta.hs + sum(coef[1:K] ^ 2 / lambda.hs) / 2)
  nu.hs <- 1 / rgamma(K, shape = 1,
                      rate = 1 + 1 / lambda.hs)
  zeta.hs <- 1 / rgamma(1, shape = 1,
                        rate = 1 + 1 / tau.hs)

  V_prior_inv <- 1 / (lambda.hs * tau.hs)
  V_prior_inv[V_prior_inv > 1e+10] <- 1e+10

  return(list(V_prior_inv = V_prior_inv,
              tau.hs = tau.hs,
              lambda.hs = lambda.hs,
              nu.hs = nu.hs,
              zeta.hs = zeta.hs))
}

