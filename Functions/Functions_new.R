## ------------------------------------------------------------------------
## 'Species popularity and research interests across the Tree of Life' 
## ------------------------------------------------------------------------

##################
# Created by Stefano Mammola
# Updated 17/03/2022
##################

## ------------------------------------------------------------------------
# 'Custom themes ggplot2'
## ------------------------------------------------------------------------

# Colors
color_models <- c("violetred4","blue")
  
custom_color <- c("grey50",    # animalia
                  "purple",    # fungi
                  "darkgreen") # plantae

# Variable names
var.names <-  c("Intercept",
                "Colorful [yes]",
                "Common name [yes]",
                "Domain [freshwater]",
                "Domain [marine]",
                "Domain [terrestrial]",
                "Harmful to humans [yes]",
                "Human use [yes]",
                "IUCN [threatened]",
                "IUCN [non-threatened]",
                "Phylogenetic distance to humans",
                "Range size",
                "Organism size",
                "Taxonomic uniqueness (Family)")

var.order <- c("Intercept",
               "Organism size",
               "Colorful [yes]",
               "Range size",
               "Taxonomic uniqueness (Family)",
               "Domain [freshwater]",
               "Domain [marine]",
               "Domain [terrestrial]",
               "IUCN [threatened]",
               "IUCN [non-threatened]",
               "Common name [yes]",
               "Human use [yes]",
               "Harmful to humans [yes]",
               "Phylogenetic distance to humans")

var.names.sub <-  c("Intercept",
                    "Colorful [yes]",
                    "Common name [yes]",
                    "Domain [freshwater]",
                    "Domain [marine]",
                    "Domain [terrestrial]",
                    "Harmful to humans [yes]",
                    "Human use [yes]",
                    "IUCN [threatened]",
                    "IUCN [non-threatened]",
                    "Range size",
                    "Organism size",
                    "Taxonomic uniqueness (Genus)")

var.order.sub <- c("Intercept",
                   "Organism size",
                   "Colorful [yes]",
                   "Range size",
                   "Taxonomic uniqueness (Genus)",
                   "Domain [freshwater]",
                   "Domain [marine]",
                   "Domain [terrestrial]",
                   "IUCN [threatened]",
                   "IUCN [non-threatened]",
                   "Common name [yes]",
                   "Human use [yes]",
                   "Harmful to humans [yes]")

# Custom theme for ggplot2
custom_theme <- theme(#text = element_text(family = "Arial"),
  axis.text = element_text(size = 9), 
  axis.title = element_text(size = 10),
  axis.line.x = element_line(color="grey10"), 
  axis.line.y = element_line(color="grey10"),
  panel.border = element_blank(),
  panel.grid.major.x = element_blank(),                                          
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.y = element_blank(),  
  #plot.margin = unit(c(.2, .2, .2, .2), units = , "cm"),
  plot.title = element_text(size = 12, vjust = 1, hjust = 0),
  legend.text = element_text(size = 9),          
  legend.title = element_blank(),                              
  legend.key = element_blank(),
  legend.background = element_rect(color = "white", 
                                   fill = "transparent", 
                                   size = 2, linetype = "white"))


# axis.text=element_text(size=10, angle=0,hjust = 0.5,colour ="grey30"),
# axis.ticks.y = element_line(color = "grey20",size = 0.7),
# axis.line.y = element_line(color = "grey20",size = 0.7, linetype = "solid"),
# axis.ticks.x = element_line(color = "grey20",size = 0.7),
# axis.line.x = element_line(color = "grey20",size = 0.7, linetype = "solid"))

## ------------------------------------------------------------------------
# 'Custom functions'
## ------------------------------------------------------------------------

#SE
my.SE <- function(x) sd(x, na.rm = TRUE)/sqrt(length(x))

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
