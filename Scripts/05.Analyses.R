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
library("ggthemes")
library("glmmTMB")
library("modEvA")
library("parameters")
library("performance")
library("png")
library("sjPlot")
library("tidylog")
library("tidyverse") 

# Custom functions & plot parameters --------------------------------------

source("Functions/Functions.r")

# Loading silhouettes ------------------------------------------------------

# Taken from PhyloPics: http://phylopic.org/
# (With open license)

animal_png <- png::readPNG("Phylopics/Animal.png")
fungi_png  <- png::readPNG("Phylopics/Fungi.png")
plant_png  <- png::readPNG("Phylopics/Plant.png")

# Loading the database  ---------------------------------------------------

db <- read.csv2(file = "./Data/SampleTREE_191022.csv", sep = '\t', dec = ',', header = TRUE, as.is = FALSE)

nrow(db)

length(na.omit(db$size_m))/nrow(db)
length(na.omit(db$size_f))/nrow(db)

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

levels(db$domain_rec) <- c("freshwater","multiple","multiple","multiple",
                           "marine", "multiple", "terrestrial", "multiple", "terrestrial")

db$domain_rec <- relevel(db$domain_rec, "multiple") #setting baseline

table(db$domain_rec)

# Homogenize distribution
db <- db %>% 
  dplyr::mutate(log_uniqueness_family = log(uniqueness_family+1),
                log_uniqueness_genus = log(uniqueness_genus+1),
                log_range_size = log(range_size+1),
                log_size_avg = log(size_avg+1),
                log_distance_hu = log(mean_divergence_time_Mya+1))

db <- db %>% 
  dplyr::mutate(scaled_uniqueness_family = scale(log_uniqueness_family, center = TRUE, scale = TRUE),
                scaled_log_distance_human = scale(log_distance_hu, center = TRUE, scale = TRUE),
                scaled_range_size = scale(log_range_size, center = TRUE, scale = TRUE),
                scaled_size = scale(log_size_avg, center = TRUE, scale = TRUE))

db$latitude <- scale(abs(db$centroid_lat))

############################################################################
############################################################################
# Modelling research interest ---------------------------------------------
############################################################################
############################################################################

dbWOS2 <- db %>% dplyr::select(WOS = Total_wos,
                               kingdom,
                               phylum,
                               class,
                               order,
                               family,
                               latitude,
                               scaled_size,
                               colorful,
                               color_blu,
                               color_red,
                               scaled_range_size,
                               domain_rec,
                               IUCN_rec,
                               scaled_uniqueness_family,
                               common_name,
                               human_use,
                               harmful_to_human,
                               scaled_log_distance_human) 

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

# Model validation
performance::check_model(M1_nbinom)    
performance::check_collinearity(M1_nbinom)

# R^2
(M1.R2 <- my.r2(M1_nbinom))
performance::r2(M1_nbinom)

############################################################################
############################################################################
# Modelling popular interest -----------------------------------------------
############################################################################
############################################################################

dbWIKI2 <- db %>% dplyr::select(wiki = total_wiki_pgviews,
                                kingdom,
                                phylum,
                                class,
                                order,
                                family,
                                latitude,
                                scaled_size,
                                colorful,
                                color_blu,
                                color_red,
                                scaled_range_size,
                                domain_rec,
                                IUCN_rec,
                                scaled_uniqueness_family,
                                common_name,
                                human_use,
                                harmful_to_human,
                                scaled_log_distance_human) 

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

diagnose(M2)
performance::check_overdispersion(M2) #model is overdispersed!
performance::check_zeroinflation(M2)

# Negative binomial
M2_nbinom <- glmmTMB::glmmTMB(model.formula.WIKI,
                              family = nbinom2, 
                              data = dbWIKI)

diagnose(M2_nbinom)
performance::check_zeroinflation(M2_nbinom)

# R^2
(M2.R2 <- my.r2(M2_nbinom))

############################################################################
############################################################################
# Saving the model output --------------------------------------------------
############################################################################
############################################################################

# Summary table
(par.M1 <- parameters::parameters(M1_nbinom))
(par.M2 <- parameters::parameters(M2_nbinom))

table.M1 <- par.M1 %>% dplyr::select(Parameter,
                         Effects,
                         Beta = Coefficient,
                         SE,
                         CI_low,
                         CI_high,
                         z,
                         p) %>% 
                         data.frame() %>% 
                         mutate_if(is.numeric, ~ round(.,3)) ; rm(par.M1)

table.M1 <- table.M1[table.M1$Effects == "fixed",] %>% 
            dplyr::select(-c(Effects)) %>% 
            na.omit()

table.M2 <- par.M2 %>% dplyr::select(Parameter,
                                     Effects,
                                     Beta = Coefficient,
                                     SE,
                                     CI_low,
                                     CI_high,
                                     z,
                                     p) %>% 
  data.frame() %>% 
  mutate_if(is.numeric, ~ round(.,3)) ; rm(par.M2)

table.M2 <- table.M2[table.M2$Effects == "fixed",] %>% 
  dplyr::select(-c(Effects)) %>% 
  na.omit()

table.M <- cbind(Model = c(rep("Web of Science",nrow(table.M1)),
                           rep("Wikipedia",nrow(table.M2))),
                 rbind(table.M1,table.M2)) ; rm(table.M1,table.M2)
                 
table.M$Parameter <- as.factor(as.character(table.M$Parameter))
table.M$Model <- as.factor(as.character(table.M$Model))

var.names <-  c("Intercept",
               "Color blue [yes]",
               "Color red [yes]",
               "Colorful [yes]",
               "Common name [yes]",
               "Domain [freshwater]",
               "Domain [marine]",
               "Domain [terrestrial]",
               "Harmful to humans [yes]",
               "Human use [yes]",
               "IUCN [endangered]",
               "IUCN [non-endangered]",
               "Latitude",
               "Phylogenetic distance to humans",
               "Range size",
               "Organism size",
               "Family uniqueness (N° species)")

levels(table.M$Parameter) <- var.names

var.order <- c("Intercept",
               "Organism size",
               "Colorful [yes]",
               "Color blue [yes]",
               "Color red [yes]",
               "Range size",
               "Latitude",
               "Family uniqueness (N° species)",
               "Domain [freshwater]",
               "Domain [marine]",
               "Domain [terrestrial]",
               "IUCN [endangered]",
               "IUCN [non-endangered]",
               "Common name [yes]",
               "Human use [yes]",
               "Harmful to humans [yes]",
               "Phylogenetic distance to humans")

table.M$Parameter <- factor(table.M$Parameter, rev(var.order)) #Sort

#Categorizing variables
var.type <- c("Intercept",
              "Ecological",
              rep("Morphological",4),
              rep("Ecological",7),
              rep("Cultural",4))
             
table.M <- cbind(Type = rep(var.type,2), table.M)

# Saving the table

write.csv(table.M,"Tables/TableS1.csv")

############################################################################
############################################################################
# Visualizing the effect ---------------------------------------------------
############################################################################
############################################################################

table.plot <- table.M[table.M$Parameter != "Intercept",] ; table.plot = droplevels(table.plot)

sign <- ifelse(table.plot$p > 0.05, "", ifelse(table.plot$p > 0.001,"", " *")) #Significance

color.axis <- c(rep(color_models[3],4),
                rep(color_models[2],8),
                rep(color_models[1],4))

(M1.2.forest_plot <- 
    table.plot %>%
    ggplot2::ggplot(aes(x = Beta, y = Parameter)) + 
    facet_wrap(. ~ Model, nrow = 1, ncol = 2) +  
    geom_vline(lty = 3, size = 0.5, col = "grey50", xintercept = 0) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high, col = Type), width = 0)+
    geom_point(aes(col = Type, fill = Type), size = 2, pch = 21) +
    geom_text(aes(col = Type),label = paste0(round(table.plot$Beta, 3), sign, sep = "  "), vjust = - 1, size = 2.5) +
    labs(x = expression(paste("Estimated beta" %+-% "95% Confidence interval")),
       y = NULL) +
    
    scale_color_manual(values = color_models)+
    scale_fill_manual(values = color_models)+
    
    #R^2
    geom_text(data = data.frame(x = 1.8, y = 14, Model = "Web of Science",
                                label = paste0("R^2 ==",round(as.numeric(M1.R2[1]),2))),
              aes(x = x, y = y, label = label),
              size = 4, parse = TRUE)+
    geom_text(data = data.frame(x = 1.8, y = 13, Model = "Web of Science",
                                label = paste0("N ==",nrow(dbWOS))),
              aes(x = x, y = y, label = label),
              size = 4, parse = TRUE) +

    geom_text(data = data.frame(x = 1.8, y = 14, Model = "Wikipedia",
                                label = paste0("R^2 ==",round(as.numeric(M2.R2[1]),2))),
              aes(x = x, y = y, label = label),
              size = 4, parse = TRUE)+
    geom_text(data = data.frame(x = 1.8, y = 13, Model = "Wikipedia",
                                label = paste0("N ==",nrow(dbWIKI))),
              aes(x = x, y = y, label = label),
              size = 4, parse = TRUE) +
    
  theme_classic() + theme(legend.position = "none",
                          axis.text = element_text(size = 12), 
                          axis.title = element_text(size = 14),
                          strip.text = element_text(size = 14),
                          axis.text.y = element_text(colour = rev(color.axis)))
)

# Variance partitioning ---------------------------------------------------

#Grouping
morpho <- "colorful + color_blu + color_red + scaled_size"
antro  <- "harmful_to_human + human_use + common_name + scaled_log_distance_human"
eco    <- "latitude + IUCN_rec + domain_rec + scaled_range_size + scaled_uniqueness_family"

#Setting formulas
formula.morpho           <- as.formula(paste0("WOS ~ ",morpho,"+",random))
formula.antro            <- as.formula(paste0("WOS ~ ",antro,"+",random))
formula.eco              <- as.formula(paste0("WOS ~ ",eco,"+",random))
formula.morpho.eco       <- as.formula(paste0("WOS ~ ",morpho,"+",eco,"+",random))
formula.morpho.antro     <- as.formula(paste0("WOS ~ ",morpho,"+",antro,"+",random))
formula.antro.eco        <- as.formula(paste0("WOS ~ ",antro,"+",eco,"+",random))
formula.morpho.antro.eco <- as.formula(paste0("WOS ~ ",morpho,"+",antro,"+",eco,"+",random))

#Fitting models
M1_VPA1 <- glmmTMB::glmmTMB(formula.morpho, family = nbinom2, data = dbWOS)#A
M1_VPA2 <- glmmTMB::glmmTMB(formula.antro, family = nbinom2, data = dbWOS)#B
M1_VPA3 <- glmmTMB::glmmTMB(formula.eco, family = nbinom2, data = dbWOS)#C
M1_VPA4 <- glmmTMB::glmmTMB(formula.morpho.eco, family = nbinom2, data = dbWOS)#AB
M1_VPA5 <- glmmTMB::glmmTMB(formula.morpho.antro, family = nbinom2, data = dbWOS)#AC
M1_VPA6 <- glmmTMB::glmmTMB(formula.antro.eco, family = nbinom2, data = dbWOS)#BC
M1_VPA7 <- glmmTMB::glmmTMB(formula.morpho.antro.eco, family = nbinom2, data = dbWOS)#ABC

#VPA
M1.VPA <- modEvA::varPart(A   = as.numeric(my.r2(M1_VPA1)$R2.marginal),
                          B   = as.numeric(my.r2(M1_VPA2)$R2.marginal),
                          C   = as.numeric(my.r2(M1_VPA3)$R2.marginal),
                          AB  = as.numeric(my.r2(M1_VPA4)$R2.marginal),
                          AC  = as.numeric(my.r2(M1_VPA5)$R2.marginal),
                          BC  = as.numeric(my.r2(M1_VPA6)$R2.marginal),
                          ABC = as.numeric(my.r2(M1_VPA7)$R2.marginal),
                          A.name = "Morphological",
                          B.name = "Ecological",
                          C.name = "Cultural", 
                          plot = TRUE, 
                          plot.unexpl = TRUE)

M1.VPA$Proportion <- round(M1.VPA$Proportion,3)
M1.VPA$Proportion  <- ifelse(M1.VPA$Proportion<0,0,M1.VPA$Proportion) #converting negative to zero
M1.random <- round(as.numeric(my.r2(M1_VPA7)[2]) - as.numeric(my.r2(M1_VPA7)[1]),3)
M1.Unexplained <- M1.VPA$Proportion[8] - M1.random

df.venn.M1 <- data.frame(x = c(3.2, 1, 2),
                      y = c(1, 1, 2.8), 
                      labels = c(M1.VPA[1,1], M1.VPA[2,1], M1.VPA[3,1]),
                      col.c = c("grey30","blue","violetred4"))

(M1.venn <- df.venn.M1 %>% ggplot2::ggplot() + 
      xlim(-3,6)+
      ylim(-1,5)+
      ggforce::geom_circle(aes(x0 = x, y0 = y, r = 1.5, fill = col.c, color = col.c), alpha = .2, size = 1, show.legend = FALSE) + 
      scale_colour_identity() + 
      scale_fill_identity()+
      annotate("text", x = df.venn.M1$x , y = df.venn.M1$y, label = df.venn.M1$labels, size = 5)+ #ABC
      annotate("text", x = 2.1, y = 1, label = M1.VPA[4,1], size = 4)+ #AB
      annotate("text", x = 1.35, y = 2,label = M1.VPA[5,1] ,size = 4)+ #BC
      annotate("text", x = 2.7, y = 2,label = M1.VPA[6,1] ,size = 4)+ #AC
      annotate("text", x = 2.1, y = 1.6,label = round(M1.VPA[7,1],2), size = 3)+ #ABC
      annotate("text", x = 4.4, y = -0.8, label ="Morphological", color = df.venn.M1$col.c[1], size = 6, fontface = "bold")+
      annotate("text", x = -0.2, y = -0.8, label ="Ecological", color = df.venn.M1$col.c[2],size = 6, fontface = "bold")+
      annotate("text", x = 2, y = 4.7, label="Cultural", color = df.venn.M1$col.c[3],size = 6, fontface = "bold") +
      annotate("text", x = 6, y = 3.8, label=paste("Unexplained = ", M1.Unexplained), color = "black",size = 4,hjust = 1) +
      annotate("text", x = 6, y = 3.5, label=paste("Random = ", M1.random), color = "black",size = 4,hjust = 1) +
      coord_fixed() + 
      theme_void())

# Variance partitioning ---------------------------------------------------

#Setting formulas
formula.morpho           <- as.formula(paste0("wiki ~ ",morpho,"+",random))
formula.antro            <- as.formula(paste0("wiki ~ ",antro,"+",random))
formula.eco              <- as.formula(paste0("wiki ~ ",eco,"+",random))
formula.morpho.eco       <- as.formula(paste0("wiki ~ ",morpho,"+",eco,"+",random))
formula.morpho.antro     <- as.formula(paste0("wiki ~ ",morpho,"+",antro,"+",random))
formula.antro.eco        <- as.formula(paste0("wiki ~ ",antro,"+",eco,"+",random))
formula.morpho.antro.eco <- as.formula(paste0("wiki ~ ",morpho,"+",antro,"+",eco,"+",random))

#Fitting models
M2_VPA1 <- glmmTMB::glmmTMB(formula.morpho, family = nbinom2, data = dbWIKI)#A
M2_VPA2 <- glmmTMB::glmmTMB(formula.antro, family = nbinom2, data = dbWIKI)#B
M2_VPA3 <- glmmTMB::glmmTMB(formula.eco, family = nbinom2, data = dbWIKI)#C
M2_VPA4 <- glmmTMB::glmmTMB(formula.morpho.eco, family = nbinom2, data = dbWIKI)#AB
M2_VPA5 <- glmmTMB::glmmTMB(formula.morpho.antro, family = nbinom2, data = dbWIKI)#AC
M2_VPA6 <- glmmTMB::glmmTMB(formula.antro.eco, family = nbinom2, data = dbWIKI)#BC
M2_VPA7 <- glmmTMB::glmmTMB(formula.morpho.antro.eco, family = nbinom2, data = dbWIKI)#ABC

#VPA
M2.VPA <- modEvA::varPart(A   = as.numeric(my.r2(M2_VPA1)$R2.marginal),
                          B   = as.numeric(my.r2(M2_VPA2)$R2.marginal),
                          C   = as.numeric(my.r2(M2_VPA3)$R2.marginal),
                          AB  = as.numeric(my.r2(M2_VPA4)$R2.marginal),
                          AC  = as.numeric(my.r2(M2_VPA5)$R2.marginal),
                          BC  = as.numeric(my.r2(M2_VPA6)$R2.marginal),
                          ABC = as.numeric(my.r2(M2_VPA7)$R2.marginal),
                          A.name = "Morphological",
                          B.name = "Ecological",
                          C.name = "Cultural", 
                          plot = TRUE, 
                          plot.unexpl = TRUE)

M2.VPA$Proportion <- round(M2.VPA$Proportion,3)
M2.VPA$Proportion  <- ifelse(M2.VPA$Proportion<0,0,M2.VPA$Proportion) #converting negative to zero
M2.random <- round(as.numeric(my.r2(M2_VPA7)[2]) - as.numeric(my.r2(M2_VPA7)[1]),3)
M2.Unexplained <- M2.VPA$Proportion[8] - M2.random

df.venn.M2 <- data.frame(x = c(3.2, 1, 2),
                      y = c(1, 1, 2.8), 
                      labels = c(M2.VPA[1,1], M2.VPA[2,1], M2.VPA[3,1]),
                      col.c = c("grey30","blue","violetred4"))

(M2.venn <- df.venn.M2 %>% ggplot2::ggplot() + 
    xlim(-3,6)+
    ylim(-1,5)+
    ggforce::geom_circle(aes(x0 = x, y0 = y, r = 1.5, fill = col.c, color = col.c), alpha = .2, size = 1, show.legend = FALSE) + 
    scale_colour_identity() + 
    scale_fill_identity()+
    annotate("text", x = df.venn.M2$x , y = df.venn.M2$y, label = df.venn.M2$labels, size = 5)+ #ABC
    annotate("text", x = 2.1, y = 1, label = M2.VPA[4,1], size = 4)+ #AB
    annotate("text", x = 1.35, y = 2,label = M2.VPA[5,1] ,size = 4)+ #BC
    annotate("text", x = 2.7, y = 2,label = M2.VPA[6,1] ,size = 4)+ #AC
    annotate("text", x = 2.1, y = 1.6,label = round(M2.VPA[7,1],2), size = 3)+ #ABC
    annotate("text", x = 4.4, y = -0.8, label ="Morphological", color = df.venn.M2$col.c[1], size = 6, fontface = "bold")+
    annotate("text", x = -0.2, y = -0.8, label = "Ecological", color = df.venn.M2$col.c[2],size = 6, fontface = "bold")+
    annotate("text", x = 2, y = 4.7, label="Cultural", color = df.venn.M2$col.c[3],size = 6, fontface = "bold") +
    annotate("text", x = 6, y = 3.8, label=paste("Unexplained = ", M2.Unexplained), color = "black",size = 4,hjust = 1) +
    annotate("text", x = 6, y = 3.5, label=paste("Random = ", M2.random), color = "black",size = 4,hjust = 1) +
    coord_fixed() + 
    theme_void())

############################################################################
############################################################################
# Modelling major phyla ----------------------------------------------------
############################################################################
############################################################################

#Modelling 
db.phyla <- db %>% dplyr::select(phylum,
                                 WOS = Total_wos,
                                 wiki = total_wiki_pgviews,
                                 class,
                                 order,
                                 family,
                                 latitude,
                                 scaled_size,
                                 colorful,
                                 color_blu,
                                 color_red,
                                 scaled_range_size,
                                 domain_rec,
                                 IUCN_rec,
                                 scaled_uniqueness_genus = log_uniqueness_genus,
                                 common_name,
                                 human_use,
                                 harmful_to_human)
                                # scaled_log_distance_human) 

db.phyla$scaled_uniqueness_genus <- scale(db.phyla$scaled_uniqueness_genus)

# Setting formula ---------------------------------------------------------

random.phyla <- "(1 | class) + (1 | order) + (1 | family)"

#Enough levels in Class?
nlevels(droplevels(db.phyla[db.phyla$phylum == "Chordata",]$class))
nlevels(droplevels(db.phyla[db.phyla$phylum == "Tracheophyta",]$class))
nlevels(droplevels(db.phyla[db.phyla$phylum == "Arthropoda",]$class))

# Setting formula
model.formula.phyla.WOS <- as.formula(paste0("WOS ~",
                                       paste(colnames(db.phyla)[7:ncol(db.phyla)], collapse = " + "),
                                       "+",
                                       random.phyla))

model.formula.phyla.wiki <- as.formula(paste0("wiki ~",
                                             paste(colnames(db.phyla)[7:ncol(db.phyla)], collapse = " + "),
                                             "+",
                                             random.phyla))

# WOS models --------------------------------------------------------------

M.WOS.chordata <- glmmTMB::glmmTMB(model.formula.phyla.WOS, family = nbinom2, 
                               data = db.phyla[db.phyla$phylum == "Chordata",],
                               control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

M.WOS.arthropoda <- glmmTMB::glmmTMB(model.formula.phyla.WOS, family = nbinom2, 
                               data = db.phyla[db.phyla$phylum == "Arthropoda",],
                               control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

M.WOS.tracheo <- glmmTMB::glmmTMB(model.formula.phyla.WOS, family = nbinom2, 
                                 data = db.phyla[db.phyla$phylum == "Tracheophyta",],
                                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

# Check
performance::check_model(M.WOS.chordata)    
performance::check_model(M.WOS.arthropoda)    
performance::check_model(M.WOS.tracheo)    

# R^2
(R2.WOS.chordata  <- my.r2(M.WOS.chordata))
(R2.WOS.arthropoda <- my.r2(M.WOS.arthropoda))
(R2.WOS.tracheo    <- my.r2(M.WOS.tracheo))

# WIKI models -------------------------------------------------------------

M.wiki.chordata <- glmmTMB::glmmTMB(model.formula.phyla.wiki, family = nbinom2, 
                                   data = db.phyla[db.phyla$phylum == "Chordata",],
                                   control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

M.wiki.arthropoda <- glmmTMB::glmmTMB(model.formula.phyla.wiki, family = nbinom2, 
                                     data = db.phyla[db.phyla$phylum == "Arthropoda",],
                                     control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

M.wiki.tracheo <- glmmTMB::glmmTMB(model.formula.phyla.wiki, family = nbinom2, 
                                  data = db.phyla[db.phyla$phylum == "Tracheophyta",],
                                  control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

# Check
performance::check_model(M.wiki.chordata)    
performance::check_model(M.wiki.arthropoda)    
performance::check_model(M.wiki.tracheo)    

# R^2
(R2.wiki.chordata   <- my.r2(M.wiki.chordata))
(R2.wiki.arthropoda <- my.r2(M.wiki.arthropoda))
(R2.wiki.tracheo    <- my.r2(M.wiki.tracheo))

# Setting tables WOS -------------------------------------------------------

# Summary table
(par.chordata.WOS   <- parameters::parameters(M.WOS.chordata))
(par.arthropoda.WOS <- parameters::parameters(M.WOS.arthropoda))
(par.tracheo.WOS    <- parameters::parameters(M.WOS.tracheo))

table.chordata.WOS <- par.chordata.WOS %>% 
  dplyr::select(Parameter,Effects,Beta = Coefficient, SE,CI_low,CI_high,z,p) %>% 
  data.frame() %>% mutate_if(is.numeric, ~ round(.,3)) ; rm(par.chordata.WOS)

table.arthropoda.WOS  <- par.arthropoda.WOS  %>% 
  dplyr::select(Parameter,Effects,Beta = Coefficient, SE,CI_low,CI_high,z,p) %>% 
  data.frame() %>% mutate_if(is.numeric, ~ round(.,3)) ; rm(par.arthropoda.WOS)

table.tracheo.WOS <- par.tracheo.WOS %>% 
  dplyr::select(Parameter,Effects,Beta = Coefficient, SE,CI_low,CI_high,z,p) %>% 
  data.frame() %>% mutate_if(is.numeric, ~ round(.,3)) ; rm(par.tracheo.WOS)

table.chordata.WOS   <- table.chordata.WOS[table.chordata.WOS$Effects == "fixed",] %>% dplyr::select(-c(Effects)) %>%  na.omit()
table.arthropoda.WOS <- table.arthropoda.WOS[table.arthropoda.WOS$Effects == "fixed",] %>% dplyr::select(-c(Effects)) %>%  na.omit()
table.tracheo.WOS    <- table.tracheo.WOS[table.tracheo.WOS$Effects == "fixed",] %>% dplyr::select(-c(Effects)) %>%  na.omit()

#adding missing factor

add <- table.chordata.WOS[8:9,]
add[1:2,2:7] <- NA

table.tracheo.WOS <- rbind(table.tracheo.WOS[1:7,],
                           add,
                           table.tracheo.WOS[8:14,]) ; rm(add)

table.sub.WOS <- cbind(Model = c(rep("Chordata",nrow(table.chordata.WOS)),
                                 rep("Arthropoda",nrow(table.arthropoda.WOS)),
                                 rep("Tracheophyta",nrow(table.tracheo.WOS))),
                       rbind(table.chordata.WOS,table.arthropoda.WOS,table.tracheo.WOS)) ; rm(table.chordata.WOS,table.arthropoda.WOS,table.tracheo.WOS)

table.sub.WOS$Parameter <- as.factor(as.character(table.sub.WOS$Parameter))
table.sub.WOS$Model     <- as.factor(as.character(table.sub.WOS$Model))

var.names.sub <-  c("Intercept",
                    "Color blue [yes]",
                    "Color red [yes]",
                    "Colorful [yes]",
                    "Common name [yes]",
                    "Domain [freshwater]",
                    "Domain [marine]",
                    "Domain [terrestrial]",
                    "Harmful to humans [yes]",
                    "Human use [yes]",
                    "IUCN [endangered]",
                    "IUCN [non-endangered]",
                    "Latitude",
                    "Range size",
                    "Organism size",
                    "Genus uniqueness (N° species)")

levels(table.sub.WOS$Parameter) <- var.names.sub

var.order.sub <- c("Intercept",
                   "Organism size",
                   "Colorful [yes]",
                   "Color blue [yes]",
                   "Color red [yes]",
                   "Range size",
                   "Latitude",
                   "Genus uniqueness (N° species)",
                   "Domain [freshwater]",
                   "Domain [marine]",
                   "Domain [terrestrial]",
                   "IUCN [endangered]",
                   "IUCN [non-endangered]",
                   "Common name [yes]",
                   "Human use [yes]",
                   "Harmful to humans [yes]")

table.sub.WOS$Parameter <- factor(table.sub.WOS$Parameter, rev(var.order.sub)) #Sort

#Categorizing variables
var.type.sub <- c("Intercept",
              "Ecological",
              rep("Morphological",4),
              rep("Ecological",7),
              rep("Cultural",3))

table.sub.WOS <- cbind(Type = rep(var.type.sub,3), table.sub.WOS)

# Saving the table
write.csv(table.M,"Tables/TableS2_subWOS.csv")

# Plotting WOS ------------------------------------------------------------------

color.axis.sub <- c(rep(color_models[3],4),
                    rep(color_models[2],8),
                    rep(color_models[1],3))

(M.WOS.sub.forest_plot <- 
   table.sub.WOS[table.sub.WOS$Parameter != "Intercept",] %>%
   ggplot2::ggplot(aes(x = Beta, y = Parameter)) + 
   facet_wrap(. ~ Model, nrow = 1, ncol = 3) +  
   geom_vline(lty = 3, size = 0.5, col = "grey50", xintercept = 0) +
   geom_errorbar(aes(xmin = CI_low, xmax = CI_high, col = Type), width = 0)+
   geom_point(aes(col = Type, fill = Type), size = 2, pch = 21) +
   #geom_text(aes(col = Type),label = paste0(round(table.sub.WOS$Beta, 3), sign, sep = "  "), vjust = - 1, size = 2.5) +
   labs(title = "N° papers in the Web of Science",
        x = expression(paste("Estimated beta" %+-% "95% Confidence interval")),
        y = NULL) +
   
   scale_color_manual(values = color_models)+
   scale_fill_manual(values = color_models)+
   
   # #R^2
   geom_text(data = data.frame(x = -4, y = 1, Model = "Arthropoda",
                               label = paste0("R^2 ==",round(as.numeric(R2.wiki.arthropoda[1]),2))),
             aes(x = x, y = y, label = label),
             size = 3, parse = TRUE)+
   
   geom_text(data = data.frame(x = -4, y = 1, Model = "Chordata",
                                label = paste0("R^2 ==",round(as.numeric(R2.wiki.chordata[1]),2))),
              aes(x = x, y = y, label = label),
              size = 3, parse = TRUE)+
    
    geom_text(data = data.frame(x = -4, y = 1, Model = "Tracheophyta",
                                label = paste0("R^2 ==",round(as.numeric(R2.wiki.tracheo[1]),2))),
              aes(x = x, y = y, label = label),
              size = 3, parse = TRUE)+
   
   
   theme_classic() + theme(legend.position = "none",
                           axis.text = element_text(size = 12), 
                           axis.title = element_text(size = 14),
                           strip.text = element_text(size = 14),
                           axis.text.y = element_text(colour = rev(color.axis.sub)))
)

# Setting tables WIKI -------------------------------------------------------

# Summary table
(par.chordata.wiki   <- parameters::parameters(M.wiki.chordata))
(par.arthropoda.wiki <- parameters::parameters(M.wiki.arthropoda))
(par.tracheo.wiki    <- parameters::parameters(M.wiki.tracheo))

table.chordata.wiki <- par.chordata.wiki %>% 
  dplyr::select(Parameter,Effects,Beta = Coefficient, SE,CI_low,CI_high,z,p) %>% 
  data.frame() %>% mutate_if(is.numeric, ~ round(.,3)) ; rm(par.chordata.wiki)

table.arthropoda.wiki  <- par.arthropoda.wiki  %>% 
  dplyr::select(Parameter,Effects,Beta = Coefficient, SE,CI_low,CI_high,z,p) %>% 
  data.frame() %>% mutate_if(is.numeric, ~ round(.,3)) ; rm(par.arthropoda.wiki)

table.tracheo.wiki <- par.tracheo.wiki %>% 
  dplyr::select(Parameter,Effects,Beta = Coefficient, SE,CI_low,CI_high,z,p) %>% 
  data.frame() %>% mutate_if(is.numeric, ~ round(.,3)) ; rm(par.tracheo.wiki)

table.chordata.wiki   <- table.chordata.wiki[table.chordata.wiki$Effects == "fixed",] %>% dplyr::select(-c(Effects)) %>%  na.omit()
table.arthropoda.wiki <- table.arthropoda.wiki[table.arthropoda.wiki$Effects == "fixed",] %>% dplyr::select(-c(Effects)) %>%  na.omit()
table.tracheo.wiki    <- table.tracheo.wiki[table.tracheo.wiki$Effects == "fixed",] %>% dplyr::select(-c(Effects)) %>%  na.omit()

#adding missing factor

add <- table.chordata.wiki[8:9,]
add[1:2,2:7] <- NA

table.tracheo.wiki <- rbind(table.tracheo.wiki[1:7,],
                           add,
                           table.tracheo.wiki[8:14,]) ; rm(add)

table.sub.wiki <- cbind(Model = c(rep("Chordata",nrow(table.chordata.wiki)),
                                 rep("Arthropoda",nrow(table.arthropoda.wiki)),
                                 rep("Tracheophyta",nrow(table.tracheo.wiki))),
                       rbind(table.chordata.wiki,table.arthropoda.wiki,table.tracheo.wiki)) ; rm(table.chordata.wiki,table.arthropoda.wiki,table.tracheo.wiki)

table.sub.wiki$Parameter <- as.factor(as.character(table.sub.wiki$Parameter))
table.sub.wiki$Model     <- as.factor(as.character(table.sub.wiki$Model))

levels(table.sub.wiki$Parameter) <- var.names.sub
table.sub.wiki$Parameter <- factor(table.sub.wiki$Parameter, rev(var.order.sub)) #Sort

#Categorizing variables

table.sub.wiki <- cbind(Type = rep(var.type.sub,3), table.sub.wiki)

# Saving the table
write.csv(table.M,"Tables/TableS2_sub_wiki.csv")

# Plotting wiki ------------------------------------------------------------------

color.axis.sub <- c(rep(color_models[3],4),
                    rep(color_models[2],8),
                    rep(color_models[1],3))

(M.wiki.sub.forest_plot <- 
    table.sub.wiki[table.sub.wiki$Parameter != "Intercept",] %>%
    ggplot2::ggplot(aes(x = Beta, y = Parameter)) + 
    facet_wrap(. ~ Model, nrow = 1, ncol = 3) +  
    geom_vline(lty = 3, size = 0.5, col = "grey50", xintercept = 0) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high, col = Type), width = 0)+
    geom_point(aes(col = Type, fill = Type), size = 2, pch = 21) +
    #geom_text(aes(col = Type),label = paste0(round(table.sub.wiki$Beta, 3), sign, sep = "  "), vjust = - 1, size = 2.5) +
    labs(title = "N° views in Wikipedia",
         x = expression(paste("Estimated beta" %+-% "95% Confidence interval")),
         y = NULL) +
    
    scale_color_manual(values = color_models)+
    scale_fill_manual(values = color_models)+
    
    # #R^2
    geom_text(data = data.frame(x = -1, y = 1, Model = "Arthropoda",
                                label = paste0("R^2 ==",round(as.numeric(R2.wiki.arthropoda[1]),2))),
              aes(x = x, y = y, label = label),
              size = 3, parse = TRUE)+
    
    geom_text(data = data.frame(x = -1, y = 1, Model = "Chordata",
                                label = paste0("R^2 ==",round(as.numeric(R2.wiki.chordata[1]),2))),
              aes(x = x, y = y, label = label),
              size = 3, parse = TRUE)+
    
    geom_text(data = data.frame(x = 3, y = 1, Model = "Tracheophyta",
                                label = paste0("R^2 ==",round(as.numeric(R2.wiki.tracheo[1]),2))),
              aes(x = x, y = y, label = label),
              size = 3, parse = TRUE)+
    
    
    theme_classic() + theme(legend.position = "none",
                            axis.text = element_text(size = 12), 
                            axis.title = element_text(size = 14),
                            strip.text = element_text(size = 14),
                            axis.text.y = element_text(colour = rev(color.axis.sub)))
)

####################################################################################
####################################################################################
####################################################################################
# Figures --------------------------------------------------------------------------
####################################################################################
####################################################################################
####################################################################################

# Figure 1 ----------------------------------------------------------------

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

centroid <- db %>% dplyr::group_by(phylum) %>% summarize(mean_WOS = mean( log(Total_wos+1), na.rm = TRUE),
                                              mean_WIKI = mean( log(total_wiki_pgviews+1), na.rm = TRUE))

#Can you include a dashed int=0,slope=1 line in the scatterplot between Wikipedia and WOS interest? (interesting result when looking at the colour of the points, meets the expectations).

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
pdf(file = "Figures/Figure_2.pdf", width = 11, height = 11)

ggpubr::ggarrange(M1.2.forest_plot,
                  ggpubr::ggarrange(M1.venn, M2.venn, 
                                    ncol = 2, hjust = -5, vjust = 4,
                                    labels = c("B", "C")),
                  common.legend = FALSE,
                  hjust = -5,
                  #align = "h",
                  labels = c("A", ""),
                  ncol=1, nrow=2) 

dev.off()

# Figure 3 ----------------------------------------------------------------

pdf(file = "Figures/Figure_3.pdf", width = 10, height = 11)
ggpubr::ggarrange(M.WOS.sub.forest_plot,
                  M.wiki.sub.forest_plot,
                  common.legend = FALSE,
                  hjust = -5,
                  align = "v",
                  labels = c("A", "B"),
                  ncol=1, nrow=2) 
dev.off()

# Figure 4 ----------------------------------------------------------------
(plot4a <- db %>% 
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

(plot4b <- db %>% 
    drop_na(kingdom, phylum, Total_wos) %>% 
    ggplot(aes(y = log(total_wiki_pgviews+1), x = phylum, fill = kingdom, fill = kingdom, color = kingdom)) +
    geom_point(position = position_jitter(width = 0.35), size = 2, alpha = 0.3) +
    geom_boxplot(width = .8, outlier.shape = NA, alpha = 0.2, col = "grey20") +
    labs(y = "N° views in Wikipedia [logarithm]", x = NULL) +
    scale_color_manual(values = custom_color)+
    scale_fill_manual(values = custom_color)+
    theme_classic() +
    custom_theme + theme(axis.text.y = element_text(size = 12)) + coord_flip())


pdf(file = "Figures/Figure_4.pdf", width = 10, height = 7)
ggpubr::ggarrange(plot4a,
                  plot4b,
                  common.legend = FALSE,
                  hjust = -5,
                  align = "v",
                  labels = c("A", "B"),
                  ncol=2, nrow=1) 
dev.off()

# Figure map --------------------------------------------------------------

# Load world map
world <- ggplot2::map_data("world")

map.total <- ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "grey10", size = 0.1) +
  geom_point(data = db,
             aes(centroid_long, centroid_lat, fill = kingdom), alpha = 0.7, shape =21, color = "black", size = 1.8) +
  scale_fill_manual("",values = custom_color) +
  ggthemes::theme_map() + theme(legend.position = "top")

map.animal <- ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "grey10", size = 0.1) +
  geom_point(data = db[db$kingdom == "Animalia",],
             aes(centroid_long, centroid_lat), alpha = 0.7, shape =21, color = "black", fill = custom_color[1], size = 1.8) +
  annotation_custom(grid::rasterGrob(animal_png), xmin = -140, xmax = -100, ymin = -10, ymax = -70)+
  ggthemes::theme_map() 

map.plantae <- ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "grey10", size = 0.1) +
  geom_point(data = db[db$kingdom == "Plantae",],
             aes(centroid_long, centroid_lat), alpha = 0.7, shape =21, color = "black", fill = custom_color[3], size = 1.8) +
  annotation_custom(grid::rasterGrob(plant_png), xmin = -140, xmax = -100, ymin = -10, ymax = -70)+
  ggthemes::theme_map()

map.fungi <- ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "grey10", size = 0.1) +
  geom_point(data = db[db$kingdom == "Fungi",],
             aes(centroid_long, centroid_lat), alpha = 0.7, shape =21, color = "black", fill = custom_color[2], size = 1.8) +
  annotation_custom(grid::rasterGrob(fungi_png), xmin = -140, xmax = -100, ymin = -10, ymax = -70)+
  ggthemes::theme_map()


pdf(file = "Figures/Figure_map.pdf", width = 10, height = 5)
ggpubr::ggarrange(map.total, map.animal, map.plantae, map.fungi,
                  common.legend = FALSE,
                  hjust = 0,
                  align = "hv",
                  labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2) 
dev.off()
