
# Function to plot draws of hyperparameters for the priors
# Author: Sebastian Luckeneder

# Input:
#   - out: output object of the clm function
#   - type: character, one of c("boxplot", "violin")
#   - subs: numeric vector of subset where each coefficient has its own parameter, e.g. c(1, 3), default = FALSE
#   - subs2: character vector of subset for specific parameters, e.g. c("tau.ng"), default = FALSE
#   - median: T/F; draws points for median in violin, default = FALSE
#   - take_log: T/F; take log of the plotted values

library(tidyr)
library(dplyr)
library(ggplot2)

plot_hyper <- function(out, type, subs = FALSE, subs2 = FALSE, median = FALSE, take_log = FALSE){

  dfs <- list()
  for (i in seq_along(names(out$shrink))){ # can someone to this w/o for loop?
    df <- out$shrink[i]
    df <- df %>% as.data.frame() %>%
      dplyr::mutate(param = names(out$shrink)[i])

    if (length(names(df)) > 2 & subs != FALSE){
      df <- df %>% dplyr::select(subs, param)
    }

    df <- df %>% tidyr::gather("param", "value", -param) %>%
      dplyr::mutate(hp = gsub('\\.[0-9]+', '', param))

    dfs[[i]] <- df
  }
  df <- do.call(rbind, dfs)

  if (subs2 != FALSE){
    df <- df %>% dplyr::filter(hp %in% subs2)
  }

  if (take_log == TRUE){
    df <- df %>% dplyr::mutate(value = log(value))
  }

  if (type == "boxplot"){
    p <- df %>% ggplot2::ggplot(aes(x = param, y = value)) +
      ggplot2::geom_boxplot(aes(color = hp)) +
      ggplot2::theme_bw()
  }

  if (type == "violin"){
    p <- df %>% ggplot2::ggplot(aes(x = param, y = value)) +
      ggplot2::geom_violin(aes(color = hp)) +
      ggplot2::theme_bw()
    if (median == TRUE){
      med <- df %>% dplyr::group_by(param, hp) %>% dplyr::summarise(grp.median = median(value))
      p <- p + geom_point(data = med, aes(x = param, y = grp.median, color = hp))
    }
  }

  return(p)

}

# # tests
# plot_hyper(out3, type = "boxplot", subs2 = "tau.ng", take_log = TRUE)
# plot_hyper(out3, type = "violin", subs = c(1, 3), subs2 = "tau.ng", take_log = TRUE, median = TRUE)
# plot_hyper(out4, type = "boxplot")
# plot_hyper(out5, type = "boxplot", take_log = TRUE)
