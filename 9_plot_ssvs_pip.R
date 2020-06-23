
# plot posterior distributions of the beta coefficients function
# Author: Sebastian Luckeneder

# Input:
#   - out: output object of the ssvs function
#   - type: character, one of c("mean", violin"), default = "mean"
#   - subs: numeric vector of subset of betas, e.g. c(1:3), default = FALSE

library(dplyr)
library(ggplot2)

plot_ssvs_pip <- function(out, type = "mean", subs = FALSE){
  
  df <- out$delta
  colnames(df) <- rep("delta", ncol(out$delta))
  colnames(df) <- paste(colnames(df), c(1:ncol(df)), sep = "_")
  df <- df %>% as.data.frame()
  
  if (subs != FALSE){
    df <- df %>% dplyr::select(subs)
  }
  
  df <- df %>% tidyr::gather("beta", "value")
  
  df$beta <- factor(df$beta, levels = unique(df$beta))
  
  if (type == "mean"){
    p <- df %>% dplyr::group_by(beta) %>%
      dplyr::summarise(PIP_mean = mean(value)) %>%
      dplyr::mutate(group = ifelse(PIP_mean > 0.5, "high", "low")) %>%
      ggplot2::ggplot(aes(x = beta, y = PIP_mean, color=group)) + 
      geom_point(size = 2) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position="none")
  }
  
  if (type == "violin"){
    mea <- df %>% dplyr::group_by(beta) %>% 
      dplyr::summarise(grp.mean = mean(value)) %>%
      dplyr::mutate(group = ifelse(grp.mean > 0.5, "high", "low"))
    
    p <- df %>% dplyr::group_by(beta) %>%
      dplyr::mutate(group = ifelse(mean(value) > 0.5, "high", "low")) %>%
      ggplot2::ggplot(aes(x = beta, y = value, color=group)) + 
      geom_violin() +
      ggplot2::geom_point(data = mea, aes(y = grp.mean)) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position="none")
  }
  
  return(p)
  
}

# # tests
# plot_ssvs_pip(out_ssvs, type = "mean")
# plot_ssvs_pip(out_ssvs_fe, type = "mean") + ggplot2::theme(axis.text.x = element_text(angle = 90))
# plot_ssvs_pip(out_ssvs_fe, type = "mean", subs = c(1:10))
# plot_ssvs_pip(out_ssvs_fe, type = "violin") + ggplot2::theme(axis.text.x = element_text(angle = 90))


