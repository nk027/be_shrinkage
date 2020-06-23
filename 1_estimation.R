
# Linear model estimation
# Authors: Nikolas Kuschnig, Lukas Vashold

clm <- function(
  x, # Data (y, X)
  n_draw = 2000L, n_burn = 1000L, n_thin = 1L, # MCMC settings
  tfe = FALSE, ife = FALSE, # Fixed effect settings (ordered by year, individual)
  n_time, # Number of years
  shrinkage = c("none", "NG", "DL", "HS", "SSVS"), # Choose the prior
  sigma_a = 1, sigma_b = 1, # Prior parameters
  beta_mean = 0, beta_var = 10 ^ 8,
  tau_a = 0.1, tau_b = 0.01, # Normal-Gamma
  dl_a = 0.5, # Dirichlet-Laplace
  sample_hyper = FALSE, # Sample theta
  tau0 = 0.01, tau1 = 10, const_var = 10 ^ 8, # Stochastic Search
  const_shrink = FALSE,
  verbose = TRUE) {

  # Setup ----------

  cl <- match.call()
  start_time <- Sys.time()

  n_save <- (n_draw - n_burn) / n_thin

  y <- x[, 1] # Dependent in first column
  N <- length(y)
  X <- build_X(x, const = TRUE, tfe, ife, n_time)
  K <- ncol(X)

  M <- sum(colnames(X) == "beta")
  if(const_shrink) {M <- K}


  # Priors ----------

  beta_pr_mean <- matrix(beta_mean, K, 1)
  beta_pr_var <- diag(K) * beta_var
  beta_pr_prec <- diag(K) / beta_var

  shrinkage <- match.arg(shrinkage)
  if(shrinkage == "NG") {
    source("2_NG_fun.R")
    if(is.null(tau_b) || is.null(tau_a)) {
      c.ng <- d.ng <- 0.01
      theta.ng <- 0.1
    } else {
      c.ng <- d.ng <- tau_b
      theta.ng <- tau_a
    }

    tau.ng <- rep(10, M)
    lambda.ng <- 1

    beta_pr_var <- diag(K) * c(tau.ng, rep(beta_var, K - M))
    beta_pr_prec <- diag(K) / c(tau.ng, rep(beta_var, K - M))

    tau.ng_store <- matrix(NA, nrow = n_save, ncol = M)
    lambda.ng_store <- vector("numeric", n_save)

    if(sample_hyper) {
      sample_theta <- list(scale = 0.5, accept = 0, nburn = n_burn)
      theta.ng_store <- vector("numeric", n_save)
    } else {
      sample_theta <- NULL
    }

  } else if(shrinkage == "DL") {
    source("2_DL_fun.R")
    kappa.dl <- 1
    phi.dl <- psi.dl <- rep(1, M)
    if(missing(dl_a) | is.null(dl_a)) {
      dl_a <- 1 / M
    } else {
      dl_a <- dl_a
    }
    beta_pr_var <- diag(K) *
      c(psi.dl * phi.dl^2 * kappa.dl^2, rep(beta_var, K - M))
    beta_pr_prec <- diag(K) /
      c(psi.dl * phi.dl^2 * kappa.dl^2, rep(beta_var, K - M))

    kappa.dl_store <- vector("numeric", n_save)
    phi.dl_store <- psi.dl_store <- matrix(NA, nrow = n_save, ncol = M)

    if(sample_hyper) {
      sample_a <- list(scale = 0.1, accept = 0, nburn = n_burn)
      dl_a_store <- vector("numeric", n_save)
    } else {
      sample_a <- NULL
    }

  } else if(shrinkage == "HS") {
    source("2_HS_fun.R")
    tau.hs <- 1
    lambda.hs <- matrix(1, M, 1)
    nu.hs <- matrix(1, M, 1)
    zeta.hs <- 1
    beta_pr_var <- diag(K) *
      c(lambda.hs * tau.hs, rep(beta_var, K - M))
    beta_pr_prec <- diag(K) /
      c(lambda.hs * tau.hs, rep(beta_var, K - M))

    tau.hs_store <- vector("numeric", n_save)
    lambda.hs_store <- matrix(NA, nrow = n_save, ncol = M)
    nu.hs_store <- matrix(NA, nrow = n_save, ncol = M)
    zeta.hs_store <- vector("numeric", n_save)
  } else if(shrinkage == "SSVS") {
    source("2_SSVS_fun.R")

    ssvs <- SSVSPrior(Y = as.matrix(y), X = as.matrix(X),
      cons = K - M, ndraw = n_draw, nburn = n_burn, nthin = n_thin,
      tau0 = tau0, tau1 = tau1, s0 = sigma_a, S0 = sigma_b,
      cons_var = const_var, delta.init = NULL, p.prior = NULL,
      verbose = verbose)

    timer <- Sys.time() - start_time
    if(verbose) {
      cat("Finished after ", format(round(timer, 2)), ".\n", sep = "")
    }

    colnames(ssvs[["beta"]]) <- colnames(X)

    out <- list(
      "beta" = ssvs[["beta"]],
      "sigma" = ssvs[["sigma"]],
      "delta" = ssvs[["delta"]],
      "priors" = list(
        "sigma_a" = sigma_a, "sigma_b" = sigma_b,
        "tau0" = tau0, "tau1" = tau1,
        "cons_var" = cons_var
      ),
      "meta" = list(
        "timer" = timer, "y" = y, "X" = X, "N" = N, "K" = K,
        "n_draw" = n_draw, "n_burn" = n_burn, "n_save" = n_save, "n_thin" = n_thin
      )
    )
    class(out) <- "clm"

    return(out)
  }


  # Storage ----------

  beta_store <- matrix(NA, nrow = n_save, ncol = K)
  colnames(beta_store) <- colnames(X)
  sigma_store <- vector("numeric", n_save)
  ll <- vector("numeric", n_save)


  # Starting values
  beta_draw <- as.numeric(
    mvtnorm::rmvnorm(1L, mean = beta_pr_mean, sigma = beta_pr_var))
  sigma_draw <- 1 / rgamma(1, sigma_a / 2, sigma_b / 2)

  # Prep
  XX <- crossprod(X)
  Xy <- crossprod(X, y)

  if(verbose) {pb <- txtProgressBar(min = 0, max = n_draw, style = 3)}

  for(i in (1 - n_burn):(n_draw - n_burn)) {

    # Beta
    V <- solve(beta_pr_prec + 1 / sigma_draw * XX)
    b <- V %*% (beta_pr_prec %*% beta_pr_mean + 1 / sigma_draw * Xy)
    beta_draw <- as.numeric(mvtnorm::rmvnorm(1L, mean = b, sigma = V))

    # Sigma
    ess_draw <- crossprod(y - X %*% beta_draw)
    sigma_draw <- 1 / rgamma(1, sigma_a + N / 2,
      sigma_b + as.double(ess_draw) / 2)

    # Shrinkage
    if(shrinkage == "NG") {
      if(sample_hyper) {
        sample_theta[["irep"]] <- i + n_burn
        if(sample_theta[["irep"]] > 1) {
          sample_theta[["scale"]] <- shrink_obj[["scale"]]
          sample_theta[["accept"]] <- shrink_obj[["accept"]]
          theta.ng <- shrink_obj[["theta.ng"]]
        }
      }
      shrink_obj <- NGPrior(coef = beta_draw, tau.ng = tau.ng,
        lambda.ng = lambda.ng, theta.ng = theta.ng,
        c.ng = c.ng, d.ng = d.ng, cons = K - M,
        prior = NULL, sample_theta = sample_theta)
      beta_pr_prec <- diag(K) *
        c(shrink_obj[["V_prior_inv"]], rep(1 / beta_var, K - M))
      tau.ng <- shrink_obj[["tau.ng"]]
      lambda.ng <- shrink_obj[["lambda.ng"]]

    } else if(shrinkage == "DL") {
      if(sample_hyper) {
        sample_a[["irep"]] <- i + n_burn
        if(sample_a[["irep"]] > 1) {
          sample_a[["scale"]] <- shrink_obj[["scale"]]
          sample_a[["accept"]] <- shrink_obj[["accept"]]
          dl_a <- shrink_obj[["a.dl"]]
        }
      }
      shrink_obj <- DLPrior(coef = beta_draw, a.dl = dl_a, phi.dl = phi.dl,
        psi.dl = psi.dl, cons = K - M, sample_a = sample_a)
      beta_pr_prec <- diag(K) *
        c(shrink_obj[["V_prior_inv"]], rep(1 / beta_var, K - M))
      phi.dl <- shrink_obj[["phi.dl"]]
      psi.dl <- shrink_obj[["psi.dl"]]
    } else if(shrinkage == "HS") {
      shrink_obj <- HSPrior(coef = beta_draw, tau.hs = tau.hs,
        nu.hs = nu.hs, zeta.hs = zeta.hs, cons = K - M)
      beta_pr_prec <- diag(K) *
        c(shrink_obj[["V_prior_inv"]], rep(1 / beta_var, K - M))
      tau.hs <- shrink_obj[["tau.hs"]]
      nu.hs <- shrink_obj[["nu.hs"]]
      zeta.hs <- shrink_obj[["zeta.hs"]]
    }


    # Store
    if(i > 0 && i %% n_thin == 0) {

      beta_store[(i / n_thin), ] <- beta_draw
      sigma_store[(i / n_thin)] <- sigma_draw
      ll[(i / n_thin)] <- -ess_draw / (2 * sigma_draw)

      if(shrinkage == "NG") {
        tau.ng_store[(i / n_thin), ] <- shrink_obj[["tau.ng"]]
        lambda.ng_store[(i / n_thin)] <- shrink_obj[["lambda.ng"]]
        if(sample_hyper) {
          theta.ng_store[(i / n_thin)] <- shrink_obj[["theta.ng"]]
        }
      } else if(shrinkage == "DL") {
        kappa.dl_store[(i / n_thin)] <- shrink_obj[["kappa.dl"]]
        phi.dl_store[(i / n_thin), ] <- shrink_obj[["phi.dl"]]
        psi.dl_store[(i / n_thin), ] <- shrink_obj[["psi.dl"]]
        if(sample_hyper) {
          dl_a_store[(i / n_thin)] <- shrink_obj[["a.dl"]]
        }
      } else if(shrinkage == "HS") {
        tau.hs_store[(i / n_thin)] <- shrink_obj[["tau.hs"]]
        lambda.hs_store[(i / n_thin), ] <- shrink_obj[["lambda.hs"]]
        nu.hs_store[(i / n_thin), ] <- shrink_obj[["nu.hs"]]
        zeta.hs_store[(i / n_thin)] <- shrink_obj[["zeta.hs"]]
      }
    }

    if(verbose) {setTxtProgressBar(pb, (i + n_burn))}
  }

  timer <- Sys.time() - start_time

  if(verbose) {
    close(pb)
    cat("Finished after ", format(round(timer, 2)), ".\n", sep = "")
  }

  # Outputs ----------

  if(shrinkage == "NG") {
    if(sample_hyper) {
      shrink_out <- list(
        "tau.ng" = tau.ng_store,
        "lambda.ng" = lambda.ng_store,
        "theta.ng" = theta.ng_store,
        "MH.accept" = shrink_obj[["accept"]] / n_draw
      )
    } else {
      shrink_out <- list(
        "tau.ng" = tau.ng_store,
        "lambda.ng" = lambda.ng_store
      )
    }
  } else if(shrinkage == "DL") {
    if(sample_hyper) {
      shrink_out <- list(
        "kappa.dl" = kappa.dl_store,
        "phi.dl" = phi.dl_store,
        "psi.dl" = psi.dl_store,
        "dl_a" = dl_a_store,
        "MH.accept" = shrink_obj[["accept"]] / n_draw
      )
    } else {
      shrink_out <- list(
        "kappa.dl" = kappa.dl_store,
        "phi.dl" = phi.dl_store,
        "psi.dl" = psi.dl_store
      )
    }
  } else if(shrinkage == "HS") {
    shrink_out <- list(
      "tau.hs" = tau.hs_store,
      "lambda.hs" = lambda.hs_store,
      "nu.hs" = nu.hs_store,
      "zeta.hs" = zeta.hs_store
    )
  } else {
    shrink_out <- NULL
  }

  out <- list(
    "beta" = beta_store,
    "sigma" = sigma_store,
    "ll" = ll,
    "priors" = list(
      "sigma_a" = sigma_a, "sigma_b" = sigma_b,
      "beta_mean" = beta_mean, "beta_var" = beta_var
    ),
    "shrink" = shrink_out,
    "meta" = list(
      "timer" = timer, "y" = y, "X" = X, "N" = N, "K" = K,
      "n_draw" = n_draw, "n_burn" = n_burn, "n_save" = n_save, "n_thin" = n_thin
    )
  )
  class(out) <- "clm"

  return(out)
}


build_X <- function(
  x, const = TRUE,
  tfe = FALSE, ife = FALSE, n_time) {

  N <- nrow(x)
  X <- X_pre <- x[, -1]
  colnames(X) <- rep("beta", ncol(X))

  if(const) {X <- cbind(X, "alpha" = 1)}

  if(tfe | ife) {
    if(missing(n_time)) {stop("Please provide the number of time periods.")}
    n_ind <- N / n_time
    if(tfe) { # Time fixed effects
      TFE <- kronecker(diag(n_time), matrix(1, n_ind))[, -1]
      colnames(TFE) <- rep("tfe", ncol(TFE))
      X <- cbind(X, TFE)
    }
    if(ife) { # Individual fixed effects
      IFE <- kronecker(matrix(1, n_time), diag(n_ind))[, -1]
      colnames(IFE) <- rep("ife", ncol(IFE))
      X <- cbind(X, IFE)
    }
  }

  return(X)
}
