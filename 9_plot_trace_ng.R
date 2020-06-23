
# plot posterior distributions of the beta coefficients function
# Author: Sebastian Luckeneder

# Input:
#   - out: output object of the clm function using shrinkage == "NG"
#   - subs: subset of  c("theta.ng", "lambda.ng", "tau.ng"), default = FALSE

library(dplyr)
library(ggplot2)

plot_trace_ng <- function(out, take_log = FALSE, subs = FALSE){
  
  dfs <- list()
  for (i in 1: (length((names(out$shrink)))-1)){ # can someone to this w/o for loop?
    df <- out$shrink[i]
    df <- df %>% as.data.frame() %>% 
      dplyr::mutate(param = names(out$shrink)[i])
    
    df <- df %>% tidyr::gather("param", "value", -param) %>%
      dplyr::mutate(hp = gsub('\\.[0-9]+', '', param)) %>%
      dplyr::group_by(param) %>%
      dplyr::mutate(id = row_number())
    
    dfs[[i]] <- df
  }
  df <- do.call(rbind, dfs)
  
  if(subs != FALSE){
    df <- df %>% dplyr::filter(hp %in% subs)
  }
  
  if (take_log == TRUE){
    df <- df %>% dplyr::mutate(value = log(value))
  }
  
  p <- df %>% ggplot2::ggplot(aes(x = id, y = value)) + 
    ggplot2::geom_line(aes(color = param)) +
    ggplot2::theme_bw()
  
  return(p)
  
}


# # tests
# plot_trace_ng(out_ng, subs = "theta.ng")
# plot_trace_ng(out_ng, subs = "tau.ng", take_log = TRUE)
