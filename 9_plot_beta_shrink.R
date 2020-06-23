
# Function to compare beta posterior distributions across shrinkage priors
# Author: Sebastian Luckeneder

# Input:
#   - out: list of output objects of the clm function
#   - shrink: character vector of used shrinkage priors, e.g. c("NG", "DL", "HS")
#   - type: character, one of c("density", "boxplot", "violin")
#   - subs: character vector of subset of betas, e.g. c("beta_1", "beta_2"), default = FALSE
#   - median: T/F; draws points median in violin, default = FALSE

library(tidyr)
library(dplyr)
library(ggplot2)

plot_beta_shrink <- function(out, shrink, type, subs = FALSE, median = FALSE){

  dfs <- list()
  for (i in seq_along(out)){ # can someone to this w/o for loop?
    df <- out[[i]]$beta
    df <- df[, colnames(df) == "beta"]
    colnames(df) <- paste(colnames(df), c(1:ncol(df)), sep = "_")
    df <- df %>% as.data.frame() %>%
      dplyr::mutate(shrinkage = shrink[i])

    if (subs != FALSE){
      df <- df %>% dplyr::select(subs, shrinkage)
    }

    df <- df %>% tidyr::gather("beta", "value", -shrinkage)

    dfs[[i]] <- df
  }
  df <- do.call(rbind, dfs)


  if (type == "density"){
    p <- df %>% ggplot2::ggplot(aes(x = value, color = beta, linetype = shrinkage)) +
      geom_density() +
      ggplot2::theme_bw()
  }

  if (type == "boxplot"){
    p <- df %>% ggplot2::ggplot(aes(x = beta, y = value)) +
      ggplot2::geom_boxplot(aes(color = shrinkage)) +
      ggplot2::theme_bw()
  }

  if (type == "violin"){
    p <- df %>% ggplot2::ggplot(aes(x = beta, y = value)) +
      ggplot2::geom_violin(aes(color = shrinkage)) +
      ggplot2::theme_bw()
    if (median == TRUE){
      med <- df %>% dplyr::group_by(beta, shrinkage) %>% dplyr::summarise(grp.median = median(value))
      p <- p + geom_point(data = med, aes(x = beta, y = grp.median, color = shrinkage))
    }
  }

  return(p)

}

# # tests
# plot_beta_shrink(list(out1, out3), shrink = c("none", "NG"), type = "density", subs =  c("beta_1", "beta_2", "beta_6"))
# plot_beta_shrink(list(out1, out3, out4, out5), shrink = c("none", "NG", "DL", "HS"), type = "boxplot")
# plot_beta_shrink(list(out1, out3), shrink = c("none", "NG"), type = "boxplot", subs =  c("beta_1", "beta_2", "beta_6"))
# plot_beta_shrink(list(out1, out3), shrink = c("none", "NG"), type = "violin", subs =  c("beta_1", "beta_2", "beta_6"))
# plot_beta_shrink(list(out1, out3), shrink = c("none", "NG"), type = "violin", subs =  c("beta_1", "beta_2", "beta_6"), median = TRUE)


