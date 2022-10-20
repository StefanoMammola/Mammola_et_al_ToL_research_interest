## ------------------------------------------------------------------------
## 'Species popularity and research interests across the Tree of Life' 
## ------------------------------------------------------------------------

##################
# Created by Stefano Mammola
# Updated 17/03/2022
##################

## ------------------------------------------------------------------------
# 'Custom functions'
## ------------------------------------------------------------------------

# Cuastom theme for ggplot2
custom_theme <- theme(#text = element_text(family = "Arial"),
  axis.text = element_text(size = 10), 
  axis.title = element_text(size = 12),
  axis.line.x = element_line(color="black"), 
  axis.line.y = element_line(color="black"),
  panel.border = element_blank(),
  panel.grid.major.x = element_blank(),                                          
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.y = element_blank(),  
  plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
  plot.title = element_text(size = 18, vjust = 1, hjust = 0),
  legend.text = element_text(size = 12),          
  legend.title = element_blank(),                              
  legend.position = c(0.95, 0.15), 
  legend.key = element_blank(),
  legend.background = element_rect(color = "black", 
                                   fill = "transparent", 
                                   size = 2, linetype = "blank"))

# Custom for GGally::ggpairs

LowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    #geom_point(colour = "blue") +
    geom_smooth(method = method, color = "grey10", ...)
  p
}

cor_func <- function(data, mapping, method, ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  corr <- cor(x, y, method=method, use='complete.obs')
  
  
  ggally_text(
    label = as.character(round(corr, 2)), 
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    color = 'black',
    ...
  )
}

#Custom function for R2
my.r2 <- function(model) {
  X <- model.matrix(model)
  n <- nrow(X)
  Beta <- fixef(model)$cond
  Sf <- var(X %*% Beta)
  Sigma.list <- VarCorr(model)
  Sl <- 
    sum(
      sapply(Sigma.list$cond,
             function(Sigma)
             {
               Z <-X[,rownames(Sigma)]
               sum(diag(Z %*% Sigma %*% t(Z)))/n
             }))
  Se <- attr(Sigma.list, "sc")^2
  Sd <- 0
  total.var <- Sf + Sl + Se + Sd
  R2.marginal <- Sf / total.var 
  R2.conditional <- (Sf + Sl) / total.var 

  return(list(R2.marginal = R2.marginal, R2.conditional = R2.conditional))
}
