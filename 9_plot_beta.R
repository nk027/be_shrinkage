
# Plot posterior distributions of the beta coefficients function
# Author: Sebastian Luckeneder

# Input:
#   - out: output object of the clm function
#   - type: character, one of c("density", "joy", "boxplot", "violin")
#   - subs: character vector of subset of betas, e.g. c("beta_1", "beta_2"), default = FALSE
#   - median: T/F; draws vertical line for median, default = FALSE
#   - mean: T/F; vertical line for mean, default = FALSE
#   - cred_int: numeric (0, 1); vertical lines for hdi, default = NULL

library(tidyr)
library(dplyr)
library(ggplot2)
library(HDInterval)
library(ggjoy)

plot_beta <- function(out, type, subs = FALSE, median = FALSE, mean = FALSE, cred_int = NULL){

  df <- out$beta
  df <- df[, colnames(df) == "beta"]
  colnames(df) <- paste(colnames(df), c(1:ncol(df)), sep = "_")
  df <- df %>% as.data.frame()

  if (subs != FALSE){
    df <- df %>% dplyr::select(subs)
  }

  df <- df %>% tidyr::gather("beta", "value")

  if (type == "density"){

    p <- df %>% ggplot2::ggplot(aes(x = value, color=beta)) +
      geom_density() +
      ggplot2::theme_bw()
  }

  if (type == "joy") {
    p <- df %>% ggplot2::ggplot(aes(x = value, y = beta)) +
    ggjoy::geom_joy() +
      ggplot2::theme_bw()
  }

  if (type == "boxplot"){
    p <- df %>% ggplot2::ggplot(aes(x = beta, y = value)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw()
  }

  if (type == "violin"){
    p <- df %>% ggplot2::ggplot(aes(x = beta, y = value)) +
      ggplot2::geom_violin() +
      ggplot2::theme_bw()
    if (median == TRUE){
      med <- df %>% dplyr::group_by(beta) %>% dplyr::summarise(grp.median = median(value))
      p <- p + geom_point(data = med, aes(x = beta, y = grp.median), color = "red")
    }
  }

  if (median == TRUE & ! type %in% c("violin", "boxplot")){
    med <- df %>% dplyr::group_by(beta) %>% dplyr::summarise(grp.median = median(value))
    p <- p + geom_vline(data = med, aes(xintercept = grp.median, color = beta))
  }

  if (mean == TRUE & ! type %in% c("violin", "boxplot")){
    mea <- df %>% dplyr::group_by(beta) %>% dplyr::summarise(grp.mean = mean(value))
    p <- p + geom_vline(data = mea, aes(xintercept = grp.mean, color = beta))
  }

  if (! is.null(cred_int) & ! type %in% c("violin", "boxplot")) {
    cred <- df %>% dplyr::group_by(beta) %>% dplyr::summarise(grp.hdi.low = HDInterval::hdi(value, cred_int)[1],
                                                              grp.hdi.high = HDInterval::hdi(value, cred_int)[2])
    p <- p +
      geom_vline(data = cred, aes(xintercept = grp.hdi.low, color = beta), linetype="dashed") +
      geom_vline(data = cred, aes(xintercept = grp.hdi.high, color = beta), linetype="dashed")
  }

  return(p)

}

# # tests
# plot_beta(out1, type = "density", subs = c("beta_1", "beta_2"))
# plot_beta(out1, type = "density", subs = c("beta_2"), median = TRUE, cred_int = 0.95)
# plot_beta(out1, type = "joy")
# plot_beta(out1, type = "boxplot")
# plot_beta(out1, type = "violin", median = TRUE)


