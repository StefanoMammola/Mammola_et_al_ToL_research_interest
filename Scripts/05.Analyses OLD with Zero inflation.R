## ------------------------------------------------------------------------
## 'Species popularity and research interests across the Tree of Life' 
## ------------------------------------------------------------------------

##################
# Created by Stefano Mammola
# Updated 17/03/2022
##################

# R script to generate the analysis

####################################################################################
# Preparation ----------------------------------------------------------------------
####################################################################################

# Loading R packages ------------------------------------------------------

library("Amelia")
library("dplyr")
library("ggforce")
library("ggpubr")
library("glmmTMB")
library("modEvA")
library("parameters")
library("performance")
library("sjPlot")
library("tidylog")
library("tidyverse") 

# Custom functions & plot parameters --------------------------------------

source("Functions/Functions.r")

# Loading the database  ---------------------------------------------------

db <- read.csv2(file = "./Data/SampleTREE_191022.csv", sep = '\t', dec = ',', header = TRUE, as.is = FALSE)

# Converting factors to factors
db$model_organism   <- as.factor(db$model_organism)
db$harmful_to_human <- as.factor(db$harmful_to_human)
db$human_use        <- as.factor(db$human_use)
db$common_name      <- as.factor(db$common_name)
db$colorful         <- as.factor(db$colorful)
db$color_blu        <- as.factor(db$color_blu)
db$color_red        <- as.factor(db$color_red)

####################################################################################
# Data exploration --------------------------------------------------------
####################################################################################

# Variable distribution & collinearity ----------------------------------------------

##########################
## Continuous variables ##
##########################

## Variables Y ##

# Summary stats
range(db$Total_wos,na.rm = TRUE) ; median(db$Total_wos,na.rm = TRUE)

range(db$total_wiki_pgviews, na.rm = TRUE) ; median(db$total_wiki_pgviews, na.rm = TRUE)
range(db$wiki_langs, na.rm = TRUE) ; median(db$wiki_langs, na.rm = TRUE)
range(db$wiki_mean_month_pgviews, na.rm = TRUE) ; median(db$wiki_mean_month_pgviews, na.rm = TRUE)

# Distribution
ggplot(db, aes(x = Total_wos)) + geom_dotplot(binaxis='x', stackdir='center',dotsize=0.1) + theme_bw()

ggplot(db, aes(x = total_wiki_pgviews)) + geom_dotplot(binaxis='x', stackdir='center',dotsize=0.1) + theme_bw()
ggplot(db, aes(x = wiki_langs)) + geom_dotplot(binaxis='x', stackdir='center',dotsize=0.1) + theme_bw()
ggplot(db, aes(x = wiki_mean_month_pgviews)) + geom_dotplot(binaxis='x', stackdir='center',dotsize=0.1) + theme_bw()

# Correlations 
db %>% dplyr::select(Total_wos, total_wiki_pgviews, wiki_langs, wiki_mean_month_pgviews) %>% 
  GGally::ggpairs() + theme_classic()

## Variables X ##
# Distribution 
ggplot(db, aes(x = uniqueness_family)) + geom_dotplot(binaxis='x', stackdir='center',dotsize=0.1) + theme_bw()
ggplot(db, aes(x = range_size)) + geom_dotplot(binaxis='x', stackdir='center',dotsize=0.1) + theme_bw()
ggplot(db, aes(x = size_avg)) + geom_dotplot(binaxis='x', stackdir='center',dotsize=0.1) + theme_bw()
ggplot(db, aes(x = centroid_lat)) + geom_dotplot(binaxis='x', stackdir='center',dotsize=0.1) + theme_bw()
ggplot(db, aes(x = mean_divergence_time_Mya)) + geom_dotplot(binaxis='x', stackdir='center',dotsize=0.1) + theme_bw()

# Correlations 
db %>% dplyr::select(uniqueness_family, uniqueness_genus, range_size, size_avg, centroid_lat, mean_divergence_time_Mya) %>% 
  GGally::ggpairs() + theme_classic()

#########################
## Factorial variables ##
#########################

## Taxonomy ##
table(db$kingdom)
table(db$phylum)

nlevels(db$kingdom)
nlevels(db$phylum)
nlevels(db$class)
nlevels(db$order)
nlevels(db$family)

## Variables X ##
table(db$model_organism) #not really usable
table(db$harmful_to_human) #so so
table(db$human_use) #OK!
table(db$common_name) #OK!
table(db$colorful) #OK!
table(db$color_blu) #OK!
table(db$color_red) #OK!
table(db$IUCN) #so so
table(db$domain) # so so

# Checking relationships --------------------------------------------------

ggplot(db, aes(x = phylum, y = log(Total_wos+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = phylum, y = log(total_wiki_pgviews+1))) + geom_boxplot() + theme_bw() + coord_flip()

ggplot(db, aes(x = harmful_to_human, y = log(Total_wos+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = human_use, y = log(Total_wos+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = model_organism, y = log(Total_wos+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = common_name, y = log(Total_wos+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = colorful, y = log(Total_wos+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = color_blu, y = log(Total_wos+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = color_red, y = log(Total_wos+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = IUCN, y = log(Total_wos+1))) + geom_boxplot() + theme_bw() + coord_flip()

####################################################################################
# Data analysis ---------------------------------------------------------------------
####################################################################################

# Data transformation -----------------------------------------------------

# Balancing levels IUCN
db$IUCN_rec <- db$IUCN

levels(db$IUCN_rec) <- c("Endangered", "Unknown", "Endangered", "Endangered", "Least concern","Unknown","Least concern","Endangered")
table(db$IUCN_rec)

db$IUCN_rec <- relevel(db$IUCN_rec, "Unknown") #setting baseline

# Balancing levels Domain
db$domain_rec <- db$domain

levels(db$domain_rec) <- c("aquatic","aquatic","aquatic + terrestrial","aquatic + terrestrial",
                           "aquatic", "aquatic + terrestrial", "terrestrial", "aquatic", "terrestrial")

db$domain_rec <- relevel(db$domain_rec, "aquatic + terrestrial") #setting baseline

table(db$domain_rec)

# Homogenize distribution
db <- db %>% 
  dplyr::mutate(log_uniqueness_family = log(uniqueness_family+1),
                log_uniqueness_genus = log(uniqueness_genus+1),
                log_range_size = log(range_size+1),
                log_size_avg = log(size_avg+1),
                log_distance_hu = log(mean_divergence_time_Mya+1))

# # Invisibile??
# # scaling size and range size by group
# scaled_size <- c()
# for(i in 1:nlevels(db$phylum)) {
#   db2 <- db[db$phylum == levels(db$phylum)[i], ]
#   scaled_size <- append(scaled_size, scale(db2$log_size_avg)[,1]) }
# 
# scaled_range_size <- c()
# for(i in 1:nlevels(db$phylum)) {
#   db2 <- db[db$phylum == levels(db$phylum)[i], ]
#   scaled_range_size <- append(scaled_range_size, scale(db2$log_range_size)[,1]) }

# db <- data.frame(db, scaled_size, scaled_range_size) ; rm(db2, scaled_size, scaled_range_size, i)

db$scaled_uniqueness_family   <- scale(db$log_uniqueness_family, center = TRUE, scale = TRUE)
db$scaled_log_distance_human  <- scale(db$log_distance_hu, center = TRUE, scale = TRUE)
db$scaled_range_size          <- scale(db$log_range_size, center = TRUE, scale = TRUE)
db$scaled_size                <- scale(db$log_size_avg, center = TRUE, scale = TRUE)

# Assembling a final database ---------------------------------------------

#######################
## Research interets ##
#######################

dbWOS2 <- db %>% dplyr::select(WOS = Total_wos,
                               kingdom,
                               phylum,
                               class,
                               order,
                               family,
                               #y = centroid_lat,
                               #x = centroid_long,
                               harmful_to_human,
                               human_use,
                               common_name,
                               colorful,
                               color_blu,
                               color_red,
                               IUCN_rec,
                               domain_rec,
                               starts_with("scaled_")) 

# Missing data
Amelia::missmap(dbWOS2)
dbWOS <- na.omit(dbWOS2)
















# Setting formula ---------------------------------------------------------

# random structure
nlevels(dbWOS$phylum)
nlevels(dbWOS$class)
nlevels(dbWOS$order)
nlevels(dbWOS$family)

random <- "(1 | phylum) + (1 | class) + (1 | order) + (1 | family)"

#formula
model.formula.WOS <- as.formula(paste0("WOS ~",
                                   paste(colnames(dbWOS)[7:ncol(dbWOS)], collapse = " + "),
                                   "+",
                                   random))

# Fit the model -----------------------------------------------------------

# First model
M1 <- glmmTMB::glmmTMB(model.formula.WOS,
                       family = poisson, 
                       data = dbWOS)
diagnose(M1)
performance::check_overdispersion(M1) #Model is overdispersed
performance::check_zeroinflation(M1) #Underfitting zeros

# Negative binomial
M1_nbinom <- glmmTMB::glmmTMB(model.formula.WOS,
                              family = nbinom2, 
                              data = dbWOS)

diagnose(M1_nbinom) 

performance::check_zeroinflation(M1_nbinom)

# # Zero-inflated
# M1_zinb <- glmmTMB::glmmTMB(model.formula.WOS,
#                              ziformula = ~ 1,
#                              data = dbWOS,
#                              family = nbinom2)

#M1_hurdle <- update(M1_zinb, ziformula = ~ ., family = truncated_nbinom2) # Hurdle model

# Comparing models
AIC(M1, M1_nbinom)

# Model validation
performance::check_model(M1_hurdle)    
performance::check_collinearity(M1_hurdle)

# Summary table
(par.M1 <- parameters::parameters(M1_zinb))

table.M1 <- par.M1 %>% dplyr::select(Parameter,
                         Component,
                         Effects,
                         Beta = Coefficient,
                         SE,
                         CI_low,
                         CI_high,
                         p) %>% 
                         data.frame() %>% 
                         mutate_if(is.numeric, ~ round(.,3)) ; rm(par.M1)

table.M1 <- table.M1[table.M1$Effects == "fixed",] %>% 
            dplyr::select(-c(Effects)) %>% 
            na.omit()

table.M1$Parameter <- as.factor(as.character(table.M1$Parameter))
table.M1$Component <- as.factor(as.character(table.M1$Component))

names.M1 <-  c("Intercept",
               "Color blue",
               "Color red",
               "Colorful",
               "Common name [yes]",
               "Domain [aquatic]",
               "Domain [terrestrial]",
               "Harmful to humans [yes]",
               "Human use [yes]",
               "IUCN [endangered]",
               "IUCN [non-endangered]",
               "Phylogenetic distance to humans",
               "Range size",
               "Organism size",
               "Family uniqueness (N° species)")

levels(table.M1$Component) <- c("Conditional", "Zero-inflated")
levels(table.M1$Parameter) <- names.M1

order.M1 <- c("Intercept",
              "Organism size",
              "Color blue",
              "Color red",
              "Colorful",
              "Domain [aquatic]",
              "Domain [terrestrial]",
              "IUCN [endangered]",
              "IUCN [non-endangered]",
              "Family uniqueness (N° species)",
              "Range size",
              "Common name [yes]",
              "Harmful to humans [yes]",
              "Human use [yes]",
              "Phylogenetic distance to humans")

table.M1$Parameter <- factor(table.M1$Parameter, rev(order.M1)) #Sort

#Categorizing variables
var.type <- c("Intercept",
             rep("Cultural",3),
             rep("Morphological",3),
             rep("Ecological",4),
             "Morphological",
             rep("Ecological",2),
             "Cultural")
table.M1 <- cbind(Type = rep(var.type,2), table.M1)

# R^2
(M1.R2 <- my.r2(M1_hurdle))

# A general look
# sjPlot::plot_model(M1_hurdle, sort.est = FALSE, se = TRUE,
#                    vline.color ="grey70",
#                    show.values = TRUE, value.offset = .3) + theme_bw()

sign <- ifelse(table.M1$p > 0.05, "", ifelse(table.M1$p>0.01,"", " *")) #Significance
col_p <- ifelse(table.M1$p > 0.05, "grey5", ifelse(table.M1$Beta>0,"orange", "blue") )
col_type <- c("black",
              rep(color_models[2],3),
              rep(color_models[1],3),
              rep(color_models[3],4),
              color_models[1],
              rep(color_models[3],2),
              color_models[2])

#1 - inflated? Tutti insieme?
(M1.forest_plot <- 
    table.M1 %>% 
    ggplot2::ggplot(aes(x = Beta, y = Parameter)) + 
    facet_wrap(. ~ Component, nrow = 1, ncol = 2) +  
  geom_vline(lty = 3, size = 0.5, col = "grey50", xintercept = 0) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0, col = rep(col_type,2))+
  geom_point(size = 2, pch = 21, col = rep(col_type,2), fill = rep(col_type,2)) +
  geom_text(
    label = paste0(round(table.M1$Beta, 3), sign, sep = "  "),col = rep(col_type,2), vjust = - 1, size = 3) +
  labs(title = "Scientific interest [N° papers in the Web of Science]",
       #subtitle = paste0("[Sample size = ", nrow(dbWOS) ," observations]"),
       x = expression(paste("Estimated beta" %+-% "95% Confidence interval")),
       y = NULL) +
    geom_text(data = data.frame(x = 2.2, y = 2.4, Component = "Zero-inflated", 
                                label = paste0("R^2 ==",round(as.numeric(M1.R2[1]),2))), 
              aes(x = x, y = y, label = label), 
              size = 4, parse = TRUE)+
    geom_text(data = data.frame(x = 2.2, y = 3.4, Component = "Zero-inflated", 
                                label = paste0("N ==",nrow(dbWOS))), 
              aes(x = x, y = y, label = label), 
              size = 4, parse = TRUE) +
  theme_classic() + theme(axis.text.y = element_text(colour = rev(c("black", 
                                                                rep(color_models[1],4),
                                                                rep(color_models[3],6),
                                                                rep(color_models[2],4)))))
)

# Variance partitioning ---------------------------------------------------

#Grouping
morpho <- "colorful + color_blu + color_red + scaled_size"
antro  <- "harmful_to_human + human_use + common_name + scaled_log_distance_human"
eco    <- "IUCN_rec + domain_rec + scaled_range_size + scaled_uniqueness_family"

#Setting formulas
formula.morpho           <- as.formula(paste0("WOS ~ ",morpho,"+",random))
formula.antro            <- as.formula(paste0("WOS ~ ",antro,"+",random))
formula.eco              <- as.formula(paste0("WOS ~ ",eco,"+",random))
formula.morpho.antro     <- as.formula(paste0("WOS ~ ",morpho,"+",antro,"+",random))
formula.morpho.eco       <- as.formula(paste0("WOS ~ ",morpho,"+",eco,"+",random))
formula.antro.eco        <- as.formula(paste0("WOS ~ ",antro,"+",eco,"+",random))
formula.morpho.antro.eco <- as.formula(paste0("WOS ~ ",morpho,"+",antro,"+",eco,"+",random))

#Fitting models
M1_VPA1 <- glmmTMB::glmmTMB(formula.morpho,
                         ziformula = ~ ., family = truncated_nbinom2, data = dbWOS)
M1_VPA2 <- glmmTMB::glmmTMB(formula.antro,
                         ziformula = ~ ., family = truncated_nbinom2, data = dbWOS)
M1_VPA3 <- glmmTMB::glmmTMB(formula.eco,
                         ziformula = ~ ., family = truncated_nbinom2, data = dbWOS)
M1_VPA4 <- glmmTMB::glmmTMB(formula.morpho.antro,
                         ziformula = ~ ., family = truncated_nbinom2, data = dbWOS)
M1_VPA5 <- glmmTMB::glmmTMB(formula.morpho.eco,
                         ziformula = ~ ., family = truncated_nbinom2, data = dbWOS)
M1_VPA6 <- glmmTMB::glmmTMB(formula.antro.eco,
                         ziformula = ~ ., family = truncated_nbinom2, data = dbWOS)
M1_VPA7 <- glmmTMB::glmmTMB(formula.morpho.antro.eco,
                         ziformula = ~ ., family = truncated_nbinom2, data = dbWOS)

#VPA
M1.VPA <- modEvA::varPart(A   = as.numeric(my.r2(M1_VPA1)$R2.marginal),
                          B   = as.numeric(my.r2(M1_VPA2)$R2.marginal),
                          C   = as.numeric(my.r2(M1_VPA3)$R2.marginal),
                          AB  = as.numeric(my.r2(M1_VPA4)$R2.marginal),
                          AC  = as.numeric(my.r2(M1_VPA5)$R2.marginal),
                          BC  = as.numeric(my.r2(M1_VPA6)$R2.marginal),
                          ABC = as.numeric(my.r2(M1_VPA7)$R2.marginal),
                          A.name = "Morphological",
                          B.name = "Cultural",
                          C.name = "Ecological", 
                          plot = TRUE, 
                          plot.unexpl = TRUE)

M1.VPA$Proportion <- round(M1.VPA$Proportion,3)
M1.random <- round(as.numeric(my.r2(M1_VPA7)[2]) - as.numeric(my.r2(M1_VPA7)[1]),3)
M1.Unexplained <- M1.VPA$Proportion[8] - M1.random

df.venn <- data.frame(x = c(3.2, 1, 2),
                      y = c(1, 1, 2.8), 
                      labels = c(M1.VPA[1,1], M1.VPA[2,1], M1.VPA[3,1]),
                      col.c = color_models)

(M1.venn <- df.venn %>% ggplot2::ggplot() + 
      xlim(-2,6)+
      ylim(-2,6)+
      ggforce::geom_circle(aes(x0 = x, y0 = y, r = 1.5, fill = col.c, color = col.c), alpha = .2, size = 1, show.legend = FALSE) + 
      scale_colour_identity() + 
      scale_fill_identity()+
      annotate("text", x = df.venn$x , y = df.venn$y, label = df.venn$labels, size = 5)+ #ABC
      annotate("text", x = 2.1, y = 1, label = "0.0", size = 4)+ #AB
      annotate("text", x = 2.7, y = 2,label = "0.0" ,size = 4)+ #AC
      annotate("text", x = 1.35, y = 2,label = M1.VPA[5,1] ,size = 4)+ #BC
      annotate("text", x = 2.1, y = 1.6,label = "0.02", size = 3)+ #ABC
      annotate("text", x = 4.4, y = -0.8, label ="Morphological", color = df.venn$col.c[1], size = 4, fontface = "bold")+
      annotate("text", x = -0.2, y = -0.8, label ="Cultural", color = df.venn$col.c[2],size = 4, fontface = "bold")+
      annotate("text", x = 2, y = 4.7, label="Ecological", color = df.venn$col.c[3],size = 4, fontface = "bold") +
      annotate("text", x = 6, y = 3.8, label=paste("Unexplained = ", M1.Unexplained), color = "black",size = 4,hjust = 1) +
      annotate("text", x = 6, y = 3.5, label=paste("Random = ", M1.random), color = "black",size = 4,hjust = 1) +
      coord_fixed() + 
      theme_void())

######################
## Popular interets ##
######################

dbWIKI2 <- db %>% dplyr::select(wiki = total_wiki_pgviews,
                                kingdom,
                                phylum,
                                class,
                                order,
                                family,
                                #y = centroid_lat,
                                #x = centroid_long,
                                harmful_to_human,
                                human_use,
                                common_name,
                                colorful,
                                color_blu,
                                color_red,
                                IUCN_rec,
                                domain_rec,
                                starts_with("scaled_")) 

# Missing data
Amelia::missmap(dbWIKI2)
dbWIKI <- na.omit(dbWIKI2)


# Setting formula ---------------------------------------------------------

model.formula.WIKI <- as.formula(paste0("wiki ~",
                                   paste(colnames(dbWIKI)[7:ncol(dbWIKI)], collapse = " + "),
                                   "+ (1 | phylum) + (1 | class) + (1 | order) + (1 | family)")) #+ (1 | family)

# Fit the model -----------------------------------------------------------

# First model
M2 <- glmmTMB::glmmTMB(model.formula.WIKI,
                       family = poisson, 
                       data = dbWIKI)

performance::check_overdispersion(M2) #model is overdispersed!

# Negative binomial
M2_nbinom <- glmmTMB::glmmTMB(model.formula.WIKI,
                              family = nbinom2, 
                              data = dbWIKI)

performance::check_zeroinflation(M2_nbinom) # Yes!

# Zero-inflated
M2_zinb  <- glmmTMB::glmmTMB(model.formula.WIKI,
                              ziformula = ~ 1,
                              data = dbWIKI,
                              family = nbinom2)

# Comparing the models
AIC(M2, M2_nbinom, M2_zinb) #M2_zinb

# Model validation
performance::check_model(M2_zinb)    
performance::check_collinearity(M2_zinb)

# Summary table
(par.M2 <- parameters::parameters(M2_zinb))

table.M2 <- par.M2 %>% dplyr::select(Parameter,
                                     Component,
                                     Effects,
                                     Beta = Coefficient,
                                     SE,
                                     CI_low,
                                     CI_high,
                                     p) %>% 
  data.frame() %>% 
  mutate_if(is.numeric, ~ round(., 3)) ; rm(par.M2)

table.M2 <- table.M2[table.M2$Effects == "fixed",] %>% 
  dplyr::select(-c(Effects)) %>% 
  na.omit()

table.M2 <- table.M2[-16,]

table.M2$Parameter <- as.factor(as.character(table.M2$Parameter))

levels(table.M2$Parameter) <- names.M1
table.M2$Parameter <- factor(table.M2$Parameter, rev(order.M1)) #Sort

#Categorizing variables
var.type <- c("Intercept",
              rep("Cultural",3),
              rep("Morphological",3),
              rep("Ecological",4),
              "Morphological",
              rep("Ecological",2),
              "Cultural")
table.M2 <- cbind(Type = var.type, table.M2)

# R^2
(M2.R2 <- my.r2(M2_zinb))

# A general look
# sjPlot::plot_model(M1_hurdle, sort.est = FALSE, se = TRUE,
#                    vline.color ="grey70",
#                    show.values = TRUE, value.offset = .3) + theme_bw()
sign <- ifelse(table.M2$p > 0.05, "", ifelse(table.M2$p>0.01,"", " *")) #Significance
col_p <- ifelse(table.M2$p > 0.05, "grey5", ifelse(table.M2$Beta>0,"orange", "blue") )
col_type <- c("black",
              rep(color_models[2],3),
              rep(color_models[1],3),
              rep(color_models[3],4),
              color_models[1],
              rep(color_models[3],2),
              color_models[2])

(M2.forest_plot <- 
    table.M2 %>% 
    ggplot2::ggplot(aes(x = Beta, y = Parameter)) + 
    geom_vline(lty = 3, size = 0.5, col = "grey50", xintercept = 0) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0, col = col_type)+
    geom_point(size = 2, pch = 21, col = col_type, fill = col_type) +
    geom_text(
      label = paste0(round(table.M2$Beta, 3), sign, sep = "  "), col = col_type, vjust = - 1, size = 3) +
    labs(title = "General interest [N° views in Wikipedia]",
         subtitle = " ",
         x = expression(paste("Estimated beta" %+-% "95% Confidence interval")),
         y = NULL) +
    geom_text(data = data.frame(x = 5, y = 1.6, 
                                label = paste0("R^2 ==",round(as.numeric(M2.R2[1]),2))), 
              aes(x = x, y = y, label = label), 
              size = 4, parse = TRUE)+
    geom_text(data = data.frame(x = 5, y = 2.6,
                                label = paste0("N ==",nrow(dbWIKI))), 
              aes(x = x, y = y, label = label), 
              size = 4, parse = TRUE) +
    theme_classic() + theme(axis.text.y = element_text(colour = rev(c("black", 
                                                                      rep(color_models[1],4),
                                                                      rep(color_models[3],6),
                                                                      rep(color_models[2],4)))))
)

# Variance partitioning ---------------------------------------------------

#Grouping
morpho <- "colorful + color_blu + color_red + scaled_size"
antro  <- "harmful_to_human + human_use + common_name + scaled_log_distance_human"
eco    <- "IUCN_rec + domain_rec + scaled_range_size + scaled_uniqueness_family"

#Setting formulas
formula.morpho           <- as.formula(paste0("wiki ~ ",morpho,"+",random))
formula.antro            <- as.formula(paste0("wiki ~ ",antro,"+",random))
formula.eco              <- as.formula(paste0("wiki ~ ",eco,"+",random))
formula.morpho.antro     <- as.formula(paste0("wiki ~ ",morpho,"+",antro,"+",random))
formula.morpho.eco       <- as.formula(paste0("wiki ~ ",morpho,"+",eco,"+",random))
formula.antro.eco        <- as.formula(paste0("wiki ~ ",antro,"+",eco,"+",random))
formula.morpho.antro.eco <- as.formula(paste0("wiki ~ ",morpho,"+",antro,"+",eco,"+",random))

#Fitting models
M2_VPA1 <- glmmTMB::glmmTMB(formula.morpho,
                            ziformula = ~ 1, family = nbinom2, data = dbWIKI)
M2_VPA2 <- glmmTMB::glmmTMB(formula.antro,
                            ziformula = ~ 1, family = nbinom2, data = dbWIKI)
M2_VPA3 <- glmmTMB::glmmTMB(formula.eco,
                            ziformula = ~ 1, family = nbinom2, data = dbWIKI)
M2_VPA4 <- glmmTMB::glmmTMB(formula.morpho.antro,
                            ziformula = ~ 1, family = nbinom2, data = dbWIKI)
M2_VPA5 <- glmmTMB::glmmTMB(formula.morpho.eco,
                            ziformula = ~ 1, family = nbinom2, data = dbWIKI)
M2_VPA6 <- glmmTMB::glmmTMB(formula.antro.eco,
                            ziformula = ~ 1, family = nbinom2, data = dbWIKI)
M2_VPA7 <- glmmTMB::glmmTMB(formula.morpho.antro.eco,
                            ziformula = ~ 1, family = nbinom2, data = dbWIKI)

#VPA
M2.VPA <- modEvA::varPart(A   = as.numeric(my.r2(M2_VPA1)$R2.marginal),
                          B   = as.numeric(my.r2(M2_VPA2)$R2.marginal),
                          C   = as.numeric(my.r2(M2_VPA3)$R2.marginal),
                          AB  = as.numeric(my.r2(M2_VPA4)$R2.marginal),
                          AC  = as.numeric(my.r2(M2_VPA5)$R2.marginal),
                          BC  = as.numeric(my.r2(M2_VPA6)$R2.marginal),
                          ABC = as.numeric(my.r2(M2_VPA7)$R2.marginal),
                          A.name = "Morphological",
                          B.name = "Cultural",
                          C.name = "Ecological", 
                          plot = TRUE, 
                          plot.unexpl = TRUE)

M2.VPA$Proportion <- round(M2.VPA$Proportion,3)
M2.random <- round(as.numeric(my.r2(M2_VPA7)[2]) - as.numeric(my.r2(M2_VPA7)[1]),3)
M2.Unexplained <- M2.VPA$Proportion[8] - M1.random

df.venn <- data.frame(x = c(3.2, 1, 2),
                      y = c(1, 1, 2.8), 
                      labels = c("0.000",M2.VPA[2,1], M2.VPA[3,1]),
                      col.c = color_models)

(M2.venn <- df.venn %>% ggplot2::ggplot() + 
    xlim(-2,6)+
    ylim(-2,6)+
    geom_circle(aes(x0 = x, y0 = y, r = 1.5, fill = col.c, color = col.c), alpha = .2, size = 1, show.legend = FALSE) + 
    scale_colour_identity() + 
    scale_fill_identity()+
    annotate("text", x = df.venn$x , y = df.venn$y, label = df.venn$labels, size = 5)+ #A - B - C
    annotate("text", x = 2.1, y = 1, label = M2.VPA[4,1], size = 4)+ #AB
    annotate("text", x = 2.7, y = 2,label =  "<0.001",size = 4)+ #AC
    annotate("text", x = 1.35, y = 2,label = "0.000" ,size = 4)+ #BC
    annotate("text", x = 2.1, y = 1.6,label = "<0.01", size = 3)+ #ABC
    annotate("text", x = 4.4, y = -0.8, label ="Morphological", color = df.venn$col.c[1], size = 4, fontface = "bold")+
    annotate("text", x = -0.2, y = -0.8, label ="Cultural", color = df.venn$col.c[2],size = 4, fontface = "bold")+
    annotate("text", x = 2, y = 4.7, label="Ecological", color = df.venn$col.c[3],size = 4, fontface = "bold") +
    annotate("text", x = 6, y = 3.8, label=paste("Unexplained = ", M2.Unexplained), color = "black",size = 4,hjust = 1) +
    annotate("text", x = 6, y = 3.5, label=paste("Random = ", M2.random), color = "black",size = 4,hjust = 1) +
    coord_fixed() + 
    theme_void())

####################################################################################
# Figures --------------------------------------------------------------------------
####################################################################################

# Figure 1 ----------------------------------------------------------------

custom_color <- c("grey50","purple","green")

blankPlot <- db %>% ggplot(aes(x = log(Total_wos+1), 
                               y = log(total_wiki_pgviews+1), 
                               color = kingdom)) + 
  geom_point(alpha = 1, size = 5)+
  geom_point(color = "white", size = 6)+
  scale_color_manual("", values = custom_color)+
  geom_blank(aes(1,1))+
  theme(legend.position = c(0.5, 0.5),
        legend.background = element_rect(color = "white", 
                                         fill = "transparent", 
                                         size = 2, linetype = "blank"),
        plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )


# Summarize centroid

centroid = db %>% dplyr::group_by(phylum) %>% summarize(mean_WOS = mean( log(Total_wos+1), na.rm = TRUE),
                                              mean_WIKI = mean( log(total_wiki_pgviews+1), na.rm = TRUE))
cor_plot <- ggplot() + 
  xlab("N° papers in the Web of Science [logarithm]")+
  ylab("N° views in Wikipedia [logarithm]")+
  geom_point(data = db, aes(x = log(Total_wos+1), y = log(total_wiki_pgviews+1), color = kingdom),
             alpha = 0.4, size = 2)+
  geom_smooth(data = db, aes(x = log(Total_wos+1), y = log(total_wiki_pgviews+1)),
              method = "gam", se = T, fill = "grey30", col = "grey20", size = 0.7, alpha=0.4)+

  # geom_text(data = centroid,
  #           aes(x = mean_WOS, y = mean_WIKI,
  #               label = ifelse(mean_WOS > 0, as.character(phylum),'')),
  #           hjust = 1, vjust = -0.8, fontface="italic", size = 3)+
  # 
  # geom_point(data = centroid, 
  #            aes(x = mean_WOS, y = mean_WIKI),
  #            alpha = 1, size = 3, color = "black", shape = 2)+
  # 
  scale_x_continuous(  
    labels=c("0", "3", "6","9"),
    breaks=c(0,3,6,9))+
  
  scale_y_continuous(  
    labels=c("0", "5", "10","15","20"),
    breaks=c(0,5,10,15,20))+
  
  
  scale_color_manual("",values = custom_color)+
  theme_classic() + theme(legend.position = "none",
        axis.text.y=element_text(size=10, angle=0,hjust = 0.5,colour ="grey20"),
        axis.title.y=element_text(size=10, angle=90,colour ="black"),
        axis.ticks.y = element_line(color = "grey20",size = 0.7),
        axis.line.y = element_line(color = "grey20",size = 0.7, linetype = "solid"),
        
        axis.text.x = element_text(size = 10,angle=0,vjust=0.5,colour = "grey20"),
        axis.title.x=element_text(size=10,colour = "black"),
        axis.ticks.x = element_line(color = "grey20",size = 0.7),
        axis.line.x = element_line(color = "grey20",size = 0.7, linetype = "solid"))

density_WOS <- db %>% ggplot(aes(x = log(Total_wos+1), 
                             color = kingdom, fill = kingdom)) +
  geom_density(alpha = 0.1)+
  
  scale_x_continuous(  
    labels=c("0", "3", "6","9"),
    breaks=c(0,3,6,9))+
  
  scale_color_manual(values = custom_color)+
  scale_fill_manual(values = custom_color)+
  ylab("Density")+
  xlab(element_blank())+
  theme_classic()+
  theme(legend.position = "none", axis.text.x = element_text(size = 10, angle=0, vjust=0.5,colour = "grey30"),
        axis.ticks.y = element_line(color = "grey20",size = 0.7),
        axis.line.y = element_line(color = "grey20",size = 0.7, linetype = "solid"),
        axis.ticks.x = element_line(color = "grey20",size = 0.7),
        axis.line.x = element_line(color = "grey20",size = 0.7, linetype = "solid"))

density_WIKI <- db %>% ggplot(aes(x = log(total_wiki_pgviews+1), 
                                  color = kingdom, fill = kingdom)) +
  geom_density(alpha = 0.1)+
  scale_color_manual(values = custom_color)+
  scale_fill_manual(values = custom_color)+
  
  scale_x_continuous(  
    labels=c("0", "5", "10","15","20"),
    breaks=c(0,5,10,15,20))+

  ylab("Density")+
  xlab(element_blank())+
  theme_classic() + coord_flip() + 
  
  theme(legend.position = "none",
        axis.text.y = element_text(size = 10,angle=0,vjust=0.5),
        axis.ticks.y = element_line(color = "grey20",size = 0.7),
        axis.line.y = element_line(color = "grey20",size = 0.7, linetype = "solid"),
        axis.ticks.x = element_line(color = "grey20",size = 0.7),
        axis.line.x = element_line(color = "grey20",size = 0.7, linetype = "solid"))

pdf(file = "./Figures/Figure_1.pdf", width = 7, height = 5)

(plot_cor <- gridExtra::grid.arrange(density_WOS, blankPlot, cor_plot, density_WIKI, 
                          ncol=2, nrow=2, widths=c(3, 1.7), heights=c(1.7, 3)) )

dev.off()

# Figure 2 ----------------------------------------------------------------

#Merging
pdf(file = "Figures/Figure_2.pdf", width = 14, height = 12)

ggpubr::ggarrange(M1.forest_plot,M2.forest_plot,M1.venn,M2.venn,
                  common.legend = FALSE,
                  hjust = -5,
                  #align = "h",
                  labels = c("A", "B","C", "D"),
                  ncol=2, nrow=2) 

dev.off()

# Figure 3 ----------------------------------------------------------------
(plot3 <- db %>% 
   drop_na(kingdom, phylum, Total_wos) %>% 
   ggplot(aes(y = log(Total_wos+1), x = phylum,
              fill = kingdom, color = kingdom)) +
   geom_point(position = position_jitter(width = 0.35), size = 2, alpha = 0.3) +
   geom_boxplot(width = .8, outlier.shape = NA, alpha = 0.2, col = "grey20") +
   labs(y = "N° papers in the Web of Science [logarithm]", x = NULL) +
    scale_color_manual(values = custom_color)+
    scale_fill_manual(values = custom_color)+
    theme_classic() +
    custom_theme + theme(legend.position = "none", axis.text.y = element_text(size = 12)) + coord_flip())

(plot3 <- db %>% 
    drop_na(kingdom, phylum, Total_wos) %>% 
    ggplot(aes(y = log(total_wiki_pgviews+1), x = reorder(phylum,total_wiki_pgviews), fill = kingdom, fill = kingdom, color = kingdom)) +
    geom_point(position = position_jitter(width = 0.35), size = 2, alpha = 0.3) +
    geom_boxplot(width = .8, outlier.shape = NA, alpha = 0.2, col = "grey20") +
    labs(y = "N° views in Wikipedia [logarithm]", x = NULL) +
    scale_color_manual(values = custom_color)+
    scale_fill_manual(values = custom_color)+
    theme_classic() +
    custom_theme + theme(axis.text.y = element_text(size = 12)) + coord_flip())

