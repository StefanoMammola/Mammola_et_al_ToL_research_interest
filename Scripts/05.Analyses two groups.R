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
library("ggdist")
library("ggforce")
library("ggimage")
library("ggpubr")
library("ggthemes")
library("glmmTMB")
library("modEvA")
library("parameters")
library("performance")
library("png")
library("raster")
library("sf")
library("sjPlot")
library("tidylog")
library("tidyverse") 

# Custom functions & plot parameters --------------------------------------

source("Functions/Functions_new.r")

# Loading silhouettes ------------------------------------------------------

# Taken from PhyloPics: http://phylopic.org/
# (With open license)

animal_png <- png::readPNG("Phylopics/Animal.png")
fungi_png  <- png::readPNG("Phylopics/Fungi.png")
plant_png  <- png::readPNG("Phylopics/Plant.png")

# Loading the database  ---------------------------------------------------

db  <- read.csv(file = "./Data/full_data.csv", sep = ',', header = TRUE, as.is = FALSE)

str(db)
head(db,5)

length(na.omit(db$size_m))/nrow(db)
length(na.omit(db$size_f))/nrow(db)
length(na.omit(db$photo_google))/nrow(db)

#% missing centroids
length(na.omit(db[db$kingdom == "Animalia",]$centroid_lat))/nrow(db[db$kingdom == "Animalia",])
length(na.omit(db[db$kingdom == "Plantae",]$centroid_lat))/nrow(db[db$kingdom == "Plantae",])
length(na.omit(db[db$kingdom == "Fungi",]$centroid_lat))/nrow(db[db$kingdom == "Fungi",])
length(na.omit(db$centroid_lat))/nrow(db)

# Converting factors to factors
db$model_organism   <- as.factor(db$model_organism)
db$harmful_to_human <- as.factor(db$harmful_to_human)
db$human_use        <- as.factor(db$human_use)
db$common_name      <- as.factor(db$common_name)
db$colorful         <- as.factor(db$colorful)
db$color_blu        <- as.factor(db$color_blu)
db$color_red        <- as.factor(db$color_red)
db$biogeography     <- as.factor(db$biogeography)

####################################################################################
# Data exploration --------------------------------------------------------
####################################################################################

# Variable distribution & collinearity ----------------------------------------------

##########################
## Continuous variables ##
##########################

## Variables Y ##

# Summary stats
range(db$Total_wos,na.rm = TRUE)
median(db$Total_wos,na.rm = TRUE)
my.SE(db$Total_wos)

table(ifelse(db$Total_wos > 0, 1, 0))/nrow(db)

db[db$Total_wos == max(db$Total_wos,na.rm = T),]$genus #,most studied species: Ginko biloba
  
range(db$total_wiki_pgviews, na.rm = TRUE)
median(db$total_wiki_pgviews, na.rm = TRUE)
my.SE(db$total_wiki_pgviews)

# Distribution
ggplot(db, aes(x = Total_wos)) + geom_dotplot(binaxis='x', 
                                              stackdir='center', 
                                              dotsize = 0.5, binwidth = 50) + theme_bw()

ggplot(db, aes(x = total_wiki_pgviews)) + geom_dotplot(binaxis='x', 
                                                       stackdir='center', 
                                                       dotsize = 0.1) + theme_bw()

# Correlations 
db %>% dplyr::select(Total_wos, total_wiki_pgviews, wiki_langs, wiki_mean_month_pgviews) %>% 
  GGally::ggpairs() + theme_classic()

## Variables X ##
# Distribution 
ggplot(db, aes(x = uniqueness_family)) + geom_dotplot(binaxis='x', stackdir='center',dotsize=0.1,binwidth = 200) + theme_bw()
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
nlevels(db$biogeography)
table(db$biogeography)

## Variables X ##
table(db$model_organism) #not really usable
table(db$harmful_to_human) #so so
table(db$human_use) #OK!
table(db$common_name) #OK!
table(db$colorful) #OK!
table(db$color_blu) #so so, we use only colorful
table(db$color_red) #so so, we use only colorful
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

#Log of wos and wiki
db$log_wos  <- log(db$Total_wos+1)
db$log_wiki <- log(db$total_wiki_pgviews+1)

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
                               biogeography,
                               scaled_size,
                               colorful,
                               scaled_range_size,
                               domain_rec,
                               scaled_uniqueness_family,
                               IUCN_rec,
                               common_name,
                               human_use,
                               harmful_to_human,
                               scaled_log_distance_human) 

# Missing data
Amelia::missmap(dbWOS2)

dbWOS <- dbWOS2 %>% na.omit()

# Setting formula ---------------------------------------------------------

# random structure
random <- "(1 | phylum) + (1 | class) + (1 | order) + (1 | biogeography)"

#formula
model.formula.WOS <- as.formula(paste0("WOS ~",
                                   paste(colnames(dbWOS)[8:ncol(dbWOS)], collapse = " + "),
                                   "+",
                                   random))

# Fit the model -----------------------------------------------------------

# First model
M1 <- glmmTMB::glmmTMB(model.formula.WOS,
                       family = poisson, 
                       data = dbWOS)

# Model validation
diagnose(M1)
performance::check_overdispersion(M1) #Model is overdispersed

# Negative binomial
M1_nbinom <- glmmTMB::glmmTMB(model.formula.WOS,
                              family = nbinom2, 
                              data = dbWOS)

# Model validation
performance::check_collinearity(M1_nbinom)
performance::check_zeroinflation(M1_nbinom)
performance::check_model(M1_nbinom)

dev.off()
plot(fitted(M1_nbinom))
#identify(fitted(M1_nbinom)) #1836

#Refitting the model without the outlier
dbWOS[1836,] #elephant...
dbWOS <- dbWOS[-1836,]

M1_nbinom <- glmmTMB::glmmTMB(model.formula.WOS,
                              family = nbinom2, 
                              data = dbWOS)

# Model validation
performance::check_collinearity(M1_nbinom)
performance::check_zeroinflation(M1_nbinom)
performance::check_model(M1_nbinom)

# R^2
(M1.R2 <- my.r2(M1_nbinom))
summary(M1_nbinom)

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
                                biogeography,
                                scaled_size,
                                colorful,
                                scaled_range_size,
                                domain_rec,
                                scaled_uniqueness_family,
                                IUCN_rec,
                                common_name,
                                human_use,
                                harmful_to_human,
                                scaled_log_distance_human) 

# Missing data
Amelia::missmap(dbWIKI2)
dbWIKI <- na.omit(dbWIKI2)

# Setting formula ---------------------------------------------------------

model.formula.WIKI <- as.formula(paste0("wiki ~",
                                        paste(colnames(dbWIKI)[8:ncol(dbWIKI)], collapse = " + "),
                                        "+",
                                        random))

# Fit the model -----------------------------------------------------------

# First model
M2 <- glmmTMB::glmmTMB(model.formula.WIKI,
                       family = poisson, 
                       data = dbWIKI)

diagnose(M2)
performance::check_overdispersion(M2) #model is overdispersed!

# Negative binomial
M2_nbinom <- glmmTMB::glmmTMB(model.formula.WIKI,
                              family = nbinom2, 
                              data = dbWIKI)

# Model validation
performance::check_collinearity(M2_nbinom)
performance::check_zeroinflation(M2_nbinom)

# Zero inflated Negative binomial
M2_zinb <- glmmTMB::glmmTMB(model.formula.WIKI,
                              family = nbinom2(), 
                              zi = ~1,
                              data = dbWIKI)

AIC(M2_zinb, M2_nbinom) #improved

# Model validation
performance::check_collinearity(M2_zinb)
performance::check_model(M2_zinb)

dev.off()
plot(fitted(M2_zinb))
#identify(fitted(M2_zinb)) #1836

#Refitting the model without the outlier
dbWIKI[124,] #Crocodylia
dbWIKI <- dbWIKI[-124,]

# Refetting the Zero-inflated Negative binomial
M2_zinb <- glmmTMB::glmmTMB(model.formula.WIKI,
                            family = nbinom2(), 
                            zi = ~1,
                            data = dbWIKI)

# Model validation
performance::check_collinearity(M2_zinb)
performance::check_model(M2_zinb)

# R^2
(M2.R2 <- my.r2(M2_zinb))
summary(M2_zinb)

############################################################################
############################################################################
# Modelling residuals ------------------------------------------------------
############################################################################
############################################################################

dbRES <- db %>% dplyr::select(WOS = Total_wos,
                              WIKI = total_wiki_pgviews,
                              kingdom,
                              phylum,
                              class,
                              order,
                              family,
                              biogeography,
                              scaled_size,
                              colorful,
                              scaled_range_size,
                              domain_rec,
                              scaled_uniqueness_family,
                              IUCN_rec,
                              common_name,
                              human_use,
                              harmful_to_human,
                              scaled_log_distance_human) 

dbRES <- na.omit(dbRES)

M3.gam <- gam::gam(log(WIKI+1) ~ log(WOS+1), data = dbRES)
summary(M3.gam)
plot(M3.gam)

dbRES <- data.frame(res = residuals(M3.gam), dbRES) %>% dplyr::select(-c(WOS,WIKI))

# random structure
random <- "(1 | phylum) + (1 | class) + (1 | order) + (1 | biogeography)"

# #formula
model.formula.RES <- as.formula(paste0("res ~",
                                       paste(colnames(dbRES)[8:ncol(dbRES)], collapse = " + "),
                                       "+",
                                       random))

# Fit the model -----------------------------------------------------------

M0.Res <- glmmTMB::glmmTMB(model.formula.RES, family = gaussian,
                           data = dbRES)

#Check model
performance::check_model(M0.Res)

# R^2
(M0.R2 <- my.r2(M0.Res))
summary(M0.Res)

############################################################################
############################################################################
# Saving the model output --------------------------------------------------
############################################################################
############################################################################

# Summary table
(par.M1 <- parameters::parameters(M1_nbinom))
(par.M2 <- parameters::parameters(M2_zinb))
(par.M0 <- parameters::parameters(M0.Res))

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

table.M2 <- table.M2[1:14,]

table.M0 <- par.M0 %>% dplyr::select(Parameter,
                                     Effects,
                                     Beta = Coefficient,
                                     SE,
                                     CI_low,
                                     CI_high,
                                     z,
                                     p) %>% 
  data.frame() %>% 
  mutate_if(is.numeric, ~ round(.,3)) ; rm(par.M0)

table.M0 <- table.M0[table.M0$Effects == "fixed",] %>% 
  dplyr::select(-c(Effects)) %>% 
  na.omit()

table.M <- cbind(Model = c(rep("Web of Science",nrow(table.M1)),
                           rep("Wikipedia",nrow(table.M2)),
                           rep("Residuals",nrow(table.M0))),
                 rbind(table.M1,table.M2,table.M0)) ; rm(table.M1,table.M2,table.M0)
                 
table.M$Parameter <- as.factor(as.character(table.M$Parameter))
table.M$Model <- as.factor(as.character(table.M$Model))

levels(table.M$Parameter) <- var.names
table.M$Parameter <- factor(table.M$Parameter, rev(var.order)) #Sort

#Categorizing variables
var.type <- c("Intercept", rep("Species trait",7), 
              rep("Cultural trait",6))
             
table.M <- cbind(Type = rep(var.type,3), table.M)

# Saving the table
write.csv(table.M,"Tables/TableS1.csv")

############################################################################
############################################################################
# Visualizing the effect ---------------------------------------------------
############################################################################
############################################################################

table.plot <- table.M[table.M$Parameter != "Intercept",] ; table.plot = droplevels(table.plot)
table.plot.M1.2 <- table.plot[table.plot$Model != "Residuals",]

sign.M1.2 <- ifelse(table.plot.M1.2$p > 0.05, "", ifelse(table.plot.M1.2$p > 0.01,"", " *")) #Significance
color.axis <- c(rep(color_models[2],7),
                rep(color_models[1],6))

(M1.2.forest_plot <- 
    table.plot.M1.2 %>%
    ggplot2::ggplot(aes(x = Beta, y = Parameter)) + 
    facet_wrap(. ~ Model, nrow = 1, ncol = 2) +  
    geom_vline(lty = 3, size = 0.5, col = "grey50", xintercept = 0) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high, col = Type), width = 0)+
    geom_point(aes(col = Type, fill = Type), size = 2, pch = 21) +
    geom_text(aes(col = Type),label = paste0(round(table.plot.M1.2$Beta, 3), sign.M1.2, sep = "  "), 
              vjust = - 1, size = 3) +
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
    
  scale_y_discrete(drop=FALSE, expand=c(0.05,1))+
    
  theme_classic() + theme(legend.position = "none",
                          
                          strip.text.x = element_text(size = 12),
                          axis.text.y = element_text(colour = rev(color.axis),size = 12), 
                          axis.text.x = element_text(size = 11),
                          axis.title = element_text(size = 13))
    
)

# Visualizing the effect ---------------------------------------------------

table.plot.M0 <- table.plot[table.plot$Model == "Residuals",]

sign.M0 <- ifelse(table.plot.M0$p > 0.05, "", ifelse(table.plot.M0$p > 0.01,"", " *")) #Significance

(M0.forest_plot <- 
    table.plot.M0 %>%
    ggplot2::ggplot(aes(x = Beta, y = Parameter)) + 
    
    geom_vline(lty = 3, size = 0.5, col = "grey50", xintercept = 0) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high, col = Type), width = 0)+
    geom_point(aes(col = Type, fill = Type), size = 2, pch = 21) +
    geom_text(aes(col = Type),
              label = paste0(round(table.plot.M0$Beta, 3), sign.M0, sep = "  "), vjust = - 1, size = 2.5) +
    labs(x = expression(paste("Estimated beta" %+-% "95% Confidence interval")),
         y = NULL) +
    
    scale_color_manual(values = color_models)+
    scale_fill_manual(values = color_models)+
    
    #R^2
    geom_text(data = data.frame(x = 1, y = 2,
                                label = paste0("R^2 ==",round(as.numeric(M0.R2[1]),2))),
              aes(x = x, y = y, label = label),
              size = 4, parse = TRUE)+
    geom_text(data = data.frame(x = 1, y = 1,
                                label = paste0("N ==",nrow(dbRES))),
              aes(x = x, y = y, label = label),
              size = 4, parse = TRUE) +
    
    annotate("segment", x = 0.1, xend = 1, y = 14, yend = 14,
             color = "grey30",
             arrow = arrow(ends = "last", 
                           angle = 15, 
                           length = unit(.2,"cm")))+
    
    annotate("text", x = 0.2, y = 14.5, hjust = 0, vjust = 0.5,
             size = 3,
             color = "grey30",
             label = "Popular interest")+
    
    annotate("segment", x = -0.1, xend = -1, y = 14, yend = 14,
             color = "grey30",
             arrow = arrow(ends = "last", 
                           angle = 15, 
                           length = unit(.2,"cm")))+
    
    annotate("text", x = -0.2, y = 14.5, hjust =1, vjust = 0.5,
             size = 3,
             color = "grey30",
             label = "Scientific interest")+
    
    scale_y_discrete(drop=FALSE, expand=c(0.05,1.2))+
    
    theme_classic() + custom_theme + theme(legend.position = "none",
                            axis.text.y = element_text(colour = rev(color.axis)))
)

#ALL IN

table.plot <- table.M[table.M$Parameter != "Intercept",] ; table.plot = droplevels(table.plot)

sign <- ifelse(table.plot$p > 0.05, "", ifelse(table.plot$p > 0.01,"", " *")) #Significance

table.plot$Model <- factor(table.plot$Model, c("Web of Science",
                                              "Wikipedia",
                                              "Residuals")) #Sort
(M1.2.3.forest_plot <- 
    table.plot %>%
    ggplot2::ggplot(aes(x = Beta, y = Parameter)) + 
    facet_wrap(. ~ factor(Model, levels = c("Web of Science",
                                            "Wikipedia",
                                            "Residuals")), nrow = 1, ncol = 3) +  
    geom_vline(lty = 3, size = 0.5, col = "grey50", xintercept = 0) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high, col = Type), width = 0)+
    geom_point(aes(col = Type, fill = Type), size = 2, pch = 21) +
    #geom_text(aes(col = Type),label = paste0(round(table.plot$Beta, 3), sign, sep = "  "), vjust = - 1, size = 2.5) +
    labs(x = expression(paste("Estimated beta" %+-% "95% Confidence interval")),
         y = NULL) +
    
    scale_color_manual(values = color_models)+
    scale_fill_manual(values = color_models)+
    
    #R^2
    geom_text(data = data.frame(x = 1.8, y = 14, Model = "Web of Science",
                                label = paste0("R^2 ==",round(as.numeric(M1.R2[2]),2))),
              aes(x = x, y = y, label = label),
              size = 4, parse = TRUE)+
    geom_text(data = data.frame(x = 1.8, y = 13, Model = "Web of Science",
                                label = paste0("N ==",nrow(dbWOS))),
              aes(x = x, y = y, label = label),
              size = 4, parse = TRUE) +
    
    geom_text(data = data.frame(x = 1.8, y = 14, Model = "Wikipedia",
                                label = paste0("R^2 ==",round(as.numeric(M2.R2[2]),2))),
              aes(x = x, y = y, label = label),
              size = 4, parse = TRUE)+
    geom_text(data = data.frame(x = 1.8, y = 13, Model = "Wikipedia",
                                label = paste0("N ==",nrow(dbWIKI))),
              aes(x = x, y = y, label = label),
              size = 4, parse = TRUE) +
    
    geom_text(data = data.frame(x = 1.8, y = 14, Model = "Residuals",
                                label = paste0("R^2 ==",round(as.numeric(M0.R2[2]),2))),
              aes(x = x, y = y, label = label),
              size = 4, parse = TRUE) +
    
    geom_text(data = data.frame(x = 1.8, y = 13, Model = "Residuals",
                                label = paste0("N ==",nrow(dbRES))),
              aes(x = x, y = y, label = label),
              size = 4, parse = TRUE) +
    
    theme_classic() + custom_theme + theme(legend.position = "none",
                            axis.text.y = element_text(colour = rev(color.axis)))
)

# Variance partitioning ---------------------------------------------------

#Grouping
species  <- "colorful + scaled_size + domain_rec + scaled_range_size + scaled_uniqueness_family"
culture  <- "harmful_to_human + human_use + common_name + scaled_log_distance_human + IUCN_rec"

#Setting formulas
formula.species.WOS          <- as.formula(paste0("WOS ~ ",species,"+",random))
formula.culture.WOS          <- as.formula(paste0("WOS ~ ",culture,"+",random))
formula.species.culture.WOS  <- as.formula(paste0("WOS ~ ",species,"+",culture,"+",random))

#Fitting models
M1_VPA1 <- glmmTMB::glmmTMB(formula.species.WOS, family = nbinom2, data = dbWOS)#A
M1_VPA2 <- glmmTMB::glmmTMB(formula.culture.WOS, family = nbinom2, data = dbWOS)#B
M1_VPA3 <- glmmTMB::glmmTMB(formula.species.culture.WOS, family = nbinom2, data = dbWOS)#AB

# M1 VPA
M1.VPA <- modEvA::varPart(A   = as.numeric(my.r2(M1_VPA1)$R2.marginal),
                          B   = as.numeric(my.r2(M1_VPA2)$R2.marginal),
                          C   = NA,
                          AB  = as.numeric(my.r2(M1_VPA3)$R2.marginal),
                          AC  = NA,
                          BC  = NA,
                          ABC = NA,
                          A.name = "Species",
                          B.name = "Culture",
                          C.name = NA, 
                          plot = TRUE, 
                          plot.unexpl = TRUE)

M1.random      <- round(as.numeric(my.r2(M1_VPA3)[2]) - as.numeric(my.r2(M1_VPA3)[1]),3)
M1.Unexplained <- M1.VPA$Proportion[4] - M1.random

#Setting formulas
formula.species.WIKI          <- as.formula(paste0("wiki ~ ",species,"+",random))
formula.culture.WIKI          <- as.formula(paste0("wiki ~ ",culture,"+",random))
formula.species.culture.WIKI  <- as.formula(paste0("wiki ~ ",species,"+",culture,"+",random))

#Fitting models

M2_VPA1 <- glmmTMB::glmmTMB(formula.species.WIKI, family = nbinom2,zi = ~1, data = dbWIKI)#A
M2_VPA2 <- glmmTMB::glmmTMB(formula.culture.WIKI, family = nbinom2,zi = ~1, data = dbWIKI)#B
M2_VPA3 <- glmmTMB::glmmTMB(formula.species.culture.WIKI, family = nbinom2,,zi = ~1, data = dbWIKI)#AB

#VPA
M2.VPA <- modEvA::varPart(A   = as.numeric(my.r2(M2_VPA1)$R2.marginal),
                          B   = as.numeric(my.r2(M2_VPA2)$R2.marginal),
                          C   = NA,
                          AB  = as.numeric(my.r2(M2_VPA3)$R2.marginal),
                          AC  = NA,
                          BC  = NA,
                          ABC = NA,
                          A.name = "Species",
                          B.name = "Culture",
                          C.name = NA, 
                          plot = TRUE, 
                          plot.unexpl = TRUE)

M2.random      <- round(as.numeric(my.r2(M2_VPA3)[2]) - as.numeric(my.r2(M2_VPA3)[1]),3)
M2.Unexplained <- M2.VPA$Proportion[4] - M2.random

ord <- c("Unexplained", "Random","Species + Culture","Species","Culture")

db.VPA <- data.frame(Model = c(rep("Web of Science",5),rep("Wikipedia",5)),
                      Type = c(rep(c("Species","Culture","Species + Culture","Random", "Unexplained"),2)),
                    Values = c(M1.VPA$Proportion[c(1:3)],M1.random,M1.Unexplained,
                               M2.VPA$Proportion[c(1:3)],M2.random,M2.Unexplained)) %>% 
dplyr::mutate_if(is.character,as.factor) 

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
                                 biogeography,
                                 scaled_size,
                                 colorful,
                                 scaled_range_size,
                                 domain_rec,
                                 IUCN_rec,
                                 scaled_uniqueness_genus = log_uniqueness_genus,
                                 common_name,
                                 human_use,
                                 harmful_to_human)

db.phyla$scaled_uniqueness_genus <- scale(db.phyla$scaled_uniqueness_genus)

# Setting formula ---------------------------------------------------------

random.phyla <- "(1 | class) + (1 | order) + (1 | biogeography)"

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

M.WOS.chordata   <- glmmTMB::glmmTMB(model.formula.phyla.WOS, family = nbinom2, 
                               data = db.phyla[db.phyla$phylum == "Chordata",],
                               control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

M.WOS.arthropoda <- glmmTMB::glmmTMB(model.formula.phyla.WOS, family = nbinom2, 
                               data = db.phyla[db.phyla$phylum == "Arthropoda",],
                               control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

M.WOS.tracheo   <- glmmTMB::glmmTMB(model.formula.phyla.WOS, family = nbinom2, 
                                 data = db.phyla[db.phyla$phylum == "Tracheophyta",],
                                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

# Check
performance::check_model(M.WOS.chordata)    
performance::check_model(M.WOS.arthropoda)    
performance::check_model(M.WOS.tracheo)    

# R^2
(R2.WOS.chordata   <- my.r2(M.WOS.chordata))
(R2.WOS.arthropoda <- my.r2(M.WOS.arthropoda))
(R2.WOS.tracheo    <- my.r2(M.WOS.tracheo))

# WIKI models -------------------------------------------------------------

M.wiki.chordata <- glmmTMB::glmmTMB(model.formula.phyla.wiki, family = nbinom2, 
                                   data = db.phyla[db.phyla$phylum == "Chordata",],
                                   zi = ~1,
                                   control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

M.wiki.arthropoda <- glmmTMB::glmmTMB(model.formula.phyla.wiki, family = nbinom2, 
                                     data = db.phyla[db.phyla$phylum == "Arthropoda",],
                                     zi = ~1,
                                     control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

M.wiki.tracheo <- glmmTMB::glmmTMB(model.formula.phyla.wiki, family = nbinom2, 
                                  data = db.phyla[db.phyla$phylum == "Tracheophyta",],
                                  zi = ~1,
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

add <- table.chordata.WOS[5:6,]
add[1:2,2:7] <- NA

table.tracheo.WOS <- rbind(table.tracheo.WOS[1:4,],
                           add,
                           table.tracheo.WOS[5:nrow(table.tracheo.WOS),]) ; rm(add)

table.sub.WOS <- cbind(Model = c(rep("Chordata",nrow(table.chordata.WOS)),
                                 rep("Arthropoda",nrow(table.arthropoda.WOS)),
                                 rep("Tracheophyta",nrow(table.tracheo.WOS))),
                       rbind(table.chordata.WOS,table.arthropoda.WOS,table.tracheo.WOS)) ; rm(table.chordata.WOS,table.arthropoda.WOS,table.tracheo.WOS)

table.sub.WOS$Parameter <- as.factor(as.character(table.sub.WOS$Parameter))
table.sub.WOS$Model     <- as.factor(as.character(table.sub.WOS$Model))

levels(table.sub.WOS$Parameter) <- var.names.sub

table.sub.WOS$Parameter <- factor(table.sub.WOS$Parameter, 
                                  levels = rev(var.order.sub)) #Sort

#Categorizing variables
var.type.sub <- c("Intercept",
              rep("Species trait",6),
              rep("Cultural trait",2),
              "Species trait",rep("Cultural trait",3))

table.sub.WOS <- cbind(Type = rep(var.type.sub,3), table.sub.WOS)

# Saving the table
write.csv(table.sub.WOS,"Tables/TableS2_subWOS.csv")

# Plotting WOS ------------------------------------------------------------------
table.sub.WOS <- table.sub.WOS[table.sub.WOS$Parameter != "Intercept",] ; table.sub.WOS <- droplevels(table.sub.WOS)
sign.sub <- ifelse(table.sub.WOS$p > 0.05, "", ifelse(table.sub.WOS$p > 0.01,"", " *")) #Significance

color.axis.sub <- c(rep(color_models[2],7),
                    rep(color_models[1],6))

(M.WOS.sub.forest_plot <- 
   table.sub.WOS %>%
   ggplot2::ggplot(aes(x = Beta, y = Parameter)) + 
   facet_wrap(. ~ Model, nrow = 1, ncol = 3) +  
   geom_vline(lty = 3, size = 0.5, col = "grey50", xintercept = 0) +
   geom_errorbar(aes(xmin = CI_low, xmax = CI_high, col = Type), width = 0)+
   geom_point(aes(col = Type, fill = Type), size = 2, pch = 21) +
   geom_text(aes(col = Type),label = paste0(round(table.sub.WOS$Beta, 3), sign.sub, sep = "  "), vjust = - 1, size = 2.5) +
   labs(title = "N° papers in the Web of Science",
        x = expression(paste("Estimated beta" %+-% "95% Confidence interval")),
        y = NULL) +
   
   scale_color_manual(values = color_models)+
   scale_fill_manual(values = color_models)+
    
    xlim(-5,5)+
   
   # #R^2
   geom_text(data = data.frame(x = -4, y = 1, Model = "Arthropoda",
                               label = paste0("R^2 ==",round(as.numeric(R2.WOS.arthropoda[1]),2))),
             aes(x = x, y = y, label = label),
             size = 3, parse = TRUE)+
   
   geom_text(data = data.frame(x = -4, y = 1, Model = "Chordata",
                                label = paste0("R^2 ==",round(as.numeric(R2.WOS.chordata[1]),2))),
              aes(x = x, y = y, label = label),
              size = 3, parse = TRUE)+
    
   geom_text(data = data.frame(x = -4, y = 1, 
                               Model = "Tracheophyta",
                              label = paste0("R^2 ==",round(as.numeric(R2.WOS.tracheo[1]),2))),
              aes(x = x, y = y, label = label),
              size = 3, parse = TRUE)+
    
    ggimage::geom_phylopic(data = data.frame(x = -4, y = 10, Model = "Arthropoda",
                                             image = "593cd880-1440-4562-b589-264cc6f9e5f2"),  
                           aes(x = x, y = y, image = image),
                           size=.2, color = custom_color[1]) +
    
    ggimage::geom_phylopic(data = data.frame(x = -4, y = 10, Model = "Chordata",
                                             image = "6682f694-1c10-4386-a6e5-361581400f15"),  
                           aes(x = x, y = y, image = image),
                           size=.2, color = custom_color[1]) +
    
    ggimage::geom_phylopic(data = data.frame(x = -4, y = 10, Model = "Tracheophyta",
                                         image = "29762b5d-82b9-4fd5-908e-986b5340cadc"),  
                           aes(x = x, y = y, image = image),
                           size=.2, color = custom_color[3]) +
    scale_y_discrete(drop=FALSE, expand=c(0.05,.08))+
  
   theme_classic() + 
   custom_theme + 
   theme(legend.position = "none",
         axis.text.y = element_text(colour = rev(color.axis.sub)),
         axis.text.x = element_blank(),
         axis.title.x = element_blank())
)

# Setting tables WIKI -------------------------------------------------------

# Summary table
(par.chordata.wiki   <- parameters::parameters(M.wiki.chordata))
(par.arthropoda.wiki <- parameters::parameters(M.wiki.arthropoda))
(par.tracheo.wiki    <- parameters::parameters(M.wiki.tracheo))

table.chordata.wiki <- par.chordata.wiki %>% 
  dplyr::select(Parameter,Effects,Beta = Coefficient, SE,CI_low,CI_high,z,p) %>% 
  data.frame() %>% mutate_if(is.numeric, ~ round(.,3)) ; rm(par.chordata.wiki)

table.chordata.wiki <- table.chordata.wiki[1:13,]

table.arthropoda.wiki  <- par.arthropoda.wiki  %>% 
  dplyr::select(Parameter,Effects,Beta = Coefficient, SE,CI_low,CI_high,z,p) %>% 
  data.frame() %>% mutate_if(is.numeric, ~ round(.,3)) ; rm(par.arthropoda.wiki)

table.arthropoda.wiki <- table.arthropoda.wiki[1:13,]

table.tracheo.wiki <- par.tracheo.wiki %>% 
  dplyr::select(Parameter,Effects,Beta = Coefficient, SE,CI_low,CI_high,z,p) %>% 
  data.frame() %>% mutate_if(is.numeric, ~ round(.,3)) ; rm(par.tracheo.wiki)

table.tracheo.wiki <- table.tracheo.wiki[1:11,]

table.chordata.wiki   <- table.chordata.wiki[table.chordata.wiki$Effects == "fixed",] %>% dplyr::select(-c(Effects)) %>%  na.omit()
table.arthropoda.wiki <- table.arthropoda.wiki[table.arthropoda.wiki$Effects == "fixed",] %>% dplyr::select(-c(Effects)) %>%  na.omit()
table.tracheo.wiki    <- table.tracheo.wiki[table.tracheo.wiki$Effects == "fixed",] %>% dplyr::select(-c(Effects)) %>%  na.omit()

#adding missing factor
add <- table.chordata.wiki[7:8,]
add[1:2,2:7] <- NA

table.tracheo.wiki <- rbind(table.tracheo.wiki[1:4,],
                           add,
                           table.tracheo.wiki[5:nrow(table.tracheo.wiki),]) ; rm(add)

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
write.csv(table.sub.wiki,"Tables/TableS3_sub_wiki.csv")

# Plotting wiki ------------------------------------------------------------------
table.sub.wiki <- table.sub.wiki[table.sub.wiki$Parameter != "Intercept",] ; table.sub.wiki <- droplevels(table.sub.wiki)
sign.sub <- ifelse(table.sub.wiki$p > 0.05, "", ifelse(table.sub.wiki$p > 0.01,"", " *")) #Significance

color.axis.sub <- c(rep(color_models[2],7),
                    rep(color_models[1],5))

(M.wiki.sub.forest_plot <- 
    table.sub.wiki[table.sub.wiki$Parameter != "Intercept",] %>%
    ggplot2::ggplot(aes(x = Beta, y = Parameter)) + 
    facet_wrap(. ~ Model, nrow = 1, ncol = 3) +  
    geom_vline(lty = 3, size = 0.5, col = "grey50", xintercept = 0) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high, col = Type), width = 0)+
    geom_point(aes(col = Type, fill = Type), size = 2, pch = 21) +
    geom_text(aes(col = Type),label = paste0(round(table.sub.wiki$Beta, 3), sign.sub, sep = "  "), vjust = - 1, size = 2.5) +
    labs(title = "N° views in Wikipedia",
         x = expression(paste("Estimated beta" %+-% "95% Confidence interval")),
         y = NULL) +
    
    scale_color_manual(values = color_models)+
    scale_fill_manual(values = color_models)+
    
    xlim(-5,5)+
    
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
    
    scale_y_discrete(drop=FALSE, expand=c(0.05,.05))+
    theme_classic() + custom_theme + theme(legend.position = "none",
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

# f1.panel A

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

cor_plot <- ggplot() + 
  xlab("N° papers in the Web of Science [log-scaled axis]")+
  ylab("N° views in Wikipedia [log-scaled axis]")+
  geom_abline(intercept = 0, slope = 2, lty = 3, size = 0.5, col = "grey50")+
  geom_point(data = db, aes(x = Total_wos, y = total_wiki_pgviews, color = kingdom),
             alpha = 0.4, size = 2)+
  geom_smooth(data = db, aes(x = Total_wos, y = total_wiki_pgviews),
              method = "gam", se = TRUE, fill = "grey30", col = "grey20", size = 0.7, alpha=0.4)+
  scale_x_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 10, 100, 1000, 3000, 9000)) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 100, 1000, 100000, 50000000)) +
  scale_color_manual("",values = custom_color)+
  theme_classic() + theme(legend.position = "none") +
  custom_theme
                          
xplot <- db %>% 
  select(Total_wos, kingdom) %>%  na.omit() %>% 
  ggdensity("Total_wos", fill = "kingdom", color = "grey30",
            palette = custom_color) + 
  scale_x_continuous(trans=scales::pseudo_log_trans()) +
  theme(legend.position = "none") +
     clean_theme()

yplot <- db %>% 
  select(total_wiki_pgviews, kingdom) %>%  na.omit() %>% 
  ggdensity(x = "total_wiki_pgviews", fill = "kingdom", color = "grey30",
            palette = custom_color) + 
  scale_x_continuous(trans=scales::pseudo_log_trans()) +
  coord_flip() + 
  theme(legend.position = "none") + clean_theme()

f1.panelA <- ggarrange(xplot, blankPlot, cor_plot, yplot, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = FALSE)

# f1.panel B

dbRES <- db %>% dplyr::select(WOS = Total_wos,
                              WIKI = total_wiki_pgviews,
                              kingdom,
                              phylum) %>% na.omit()

M3.gam <- gam::gam(log(WIKI+1) ~ log(WOS+1), data = dbRES)
dbRES <- data.frame(res = residuals(M3.gam), dbRES) %>% dplyr::select(-c(WOS,WIKI))

f1.panelB <- dbRES %>%
  group_by(phylum) %>%
  mutate(median_res = median(res, na.rm = TRUE),
         n = n()) %>%
  ungroup() %>%
  arrange(desc(kingdom),median_res,phylum) %>% 
  mutate(phylum = factor(phylum, levels = unique(.$phylum))) %>% 
  ggplot(aes(x = res, y = phylum, fill = kingdom, color = kingdom)) +
  geom_point(data = ~.x %>%  filter(n == 1),
             aes(x = median_res, fill = kingdom), 
             colour = "gray30", shape = 21, alpha = .9)+
  stat_slab(
    data = ~ .x %>% filter(n > 1),
    aes(fill_ramp = stat(abs(x))),
    #, color_ramp = stat(-dnorm(x, 0, .5))),
    color = "gray15",
    size = .3,
    alpha = .9,
    expand = FALSE,
    trim = TRUE,
    height = 3
  ) +
  labs(x = "Residuals from fitted line in A", y = NULL) +
  #xlim(-7,10)+
  annotate("segment", x = 0.5, xend = 3.5, y = 30.3, yend = 30.3,
           color = "grey30",
           arrow = arrow(ends = "last", 
                         angle = 15, 
                         length = unit(.2,"cm")))+
  
  annotate("text", x = 0.5, y = 31, hjust = 0, vjust = 0.5,
           size = 3,
           color = "grey30",
           label = "Popular interest")+
  
  annotate("segment", x = -0.5, xend = -3.5, y = 30.3, yend = 30.3,
           color = "grey30",
           arrow = arrow(ends = "last", 
                         angle = 15, 
                         length = unit(.2,"cm")))+
  
  annotate("text", x = -0.5, y = 31, hjust =1, vjust = 0.5,
           size = 3,
           color = "grey30",
           label = "Scientific interest")+
  
ggimage::geom_phylopic(data = data.frame(x = 9, y = 18, kingdom = "Animalia",
                                           image = "f3309b41-d0d9-4b50-9e1a-325a0693cf5e"),  
                         aes(x = x, y = y, image = image),
                         size=.1, color = custom_color[1]) +
  
ggimage::geom_phylopic(data = data.frame(x = 9, y = 8, kingdom = "Fungi",
                                             image = "8cff2d66-6549-44d2-8304-d2dfecf53d78"),  
                           aes(x = x, y = y, image = image),
                           size=.1, color = custom_color[2]) + 
    
ggimage::geom_phylopic(data = data.frame(x = 9.5, y = 3, kingdom = "Plantae",
                                             image = "29762b5d-82b9-4fd5-908e-986b5340cadc"),  
                           aes(x = x, y = y, image = image),
                           size=.1, color = custom_color[3]) +  
      
  geom_vline(lty = 3, size = 0.5, col = "black", xintercept = 0) +
  scale_color_manual(values = custom_color)+
  scale_fill_manual(values = custom_color)+
  scale_y_discrete(drop=FALSE, expand=c(0.05,.05))+
  theme_classic() +
  custom_theme +
  theme(legend.position = "none")

# Merge
pdf(file = "Figures/Figure_1.pdf", width = 12, height = 5)
ggpubr::ggarrange(f1.panelA, f1.panelB, ncol = 2, nrow = 1, labels = c("A", "B"))
dev.off()

# Figure 2 ----------------------------------------------------------------

(plot.VPA <- db.VPA %>% 
   ggplot(aes(x=Model, y=Values, fill = factor(Type, levels = ord))) +
   geom_bar(stat="identity", color = "grey50")+
   labs(y = "Variance explained", x = " ") +
   scale_fill_manual("",values = c("grey20", "grey40", "darkorchid4", rev(color_models)))+
   theme_classic() + theme(axis.text = element_text(size = 9), 
                           axis.title = element_text(size = 10),
                           axis.line.x = element_line(color="grey10"), 
                           axis.line.y = element_line(color="grey10"),
                           legend.text = element_text(size = 9),
                           plot.margin = unit(c(1.2, .2, .2, 0.9), units = , "cm")
   ))

#Merging
pdf(file = "Figures/Figure_2.pdf", width = 13, height = 5)

ggpubr::ggarrange(M1.2.forest_plot, plot.VPA,
          ncol = 2, nrow = 1,  labels = c("A", "B"),
          widths = c(2, 1), heights = c(1, 1),
          common.legend = FALSE)

dev.off()

# Figure 3 ----------------------------------------------------------------

pdf(file = "Figures/Figure_3.pdf", width = 12, height = 9)
ggpubr::ggarrange(M.WOS.sub.forest_plot,
                  M.wiki.sub.forest_plot,
                  common.legend = FALSE,
                  hjust = -5,
                  align = "v",
                  labels = c("A", "B"),
                  ncol=1, nrow=2) 
dev.off()

# Figure S1 (Map) --------------------------------------------------------------

# Load world map
world <- ggplot2::map_data("world")

points <- db %>% 
  dplyr::select(centroid_long,centroid_lat) %>% 
  na.omit() %>% 
  sf::st_as_sf(coords=c("centroid_long", "centroid_lat"))

r <- raster::raster(points, ncols = 30, nrows = 30)
extent(r) <- c(-180,180,-90,90)

# Count the number of points on each pixel
map.total <- raster::rasterize(points, 
                               r, 
                               fun="count") %>% 
  as.data.frame(xy = TRUE) %>%
  ggplot2::ggplot() +
  geom_raster(aes(x=x, y=y, fill = layer), alpha = .9) +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "transparent", size = 0.1) +
  scale_fill_gradient("", low="grey70", high="blue", na.value="white") +
  ggthemes::theme_map() +
  theme(legend.position = c(0,0.1),
        legend.direction ="horizontal",
        legend.key.height = unit(.5, 'cm'), #change legend key height
        legend.key.width = unit(.5, 'cm'))

map.animal <- ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "transparent", size = 0.1) +
  geom_point(data = db[db$kingdom == "Animalia",],
             aes(centroid_long, centroid_lat), alpha = 0.7, shape =21, color = "black", fill = custom_color[1], size = 1.8) +
  #annotation_custom(grid::rasterGrob(animal_png), xmin = -140, xmax = -100, ymin = -10, ymax = -70)+
  ggimage::geom_phylopic(data = data.frame(x = -120, y = -45, kingdom = "Animalia",
                                           image = "f3309b41-d0d9-4b50-9e1a-325a0693cf5e"),  
                         aes(x = x, y = y, image = image),
                         size=.15, color = custom_color[1]) +
  
  ggthemes::theme_map() 

map.plantae <- ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "transparent", size = 0.1) +
  geom_point(data = db[db$kingdom == "Plantae",],
             aes(centroid_long, centroid_lat), alpha = 0.7, shape =21, color = "black", fill = custom_color[3], size = 1.8) +
  #annotation_custom(grid::rasterGrob(plant_png), xmin = -140, xmax = -100, ymin = -10, ymax = -70)+
  ggimage::geom_phylopic(data = data.frame(x = -120, y = -45, kingdom = "Plantae",
                                           image = "29762b5d-82b9-4fd5-908e-986b5340cadc"),  
                         aes(x = x, y = y, image = image),
                         size=.15, color = custom_color[3]) + 
  
  ggthemes::theme_map()

map.fungi <- ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "transparent", size = 0.1) +
  geom_point(data = db[db$kingdom == "Fungi",],
             aes(centroid_long, centroid_lat), alpha = 0.7, shape =21, color = "black", fill = custom_color[2], size = 1.8) +
  #annotation_custom(grid::rasterGrob(fungi_png), xmin = -140, xmax = -100, ymin = -10, ymax = -70)+
  ggimage::geom_phylopic(data = data.frame(x = -120, y = -45, kingdom = "Fungi",
                                           image = "8cff2d66-6549-44d2-8304-d2dfecf53d78"),  
                         aes(x = x, y = y, image = image),
                         size=.15, color = custom_color[2]) + 
  
  ggthemes::theme_map()

pdf(file = "Figures/Figure_S1.pdf", width = 9, height = 6)
ggpubr::ggarrange(map.total, map.animal, map.plantae, map.fungi,
                  common.legend = FALSE,
                  hjust = 0,
                  align = "hv",
                  labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2) 
dev.off()

# Figure S2 ----------------------------------------------------------------

(plotS2a <- db %>% 
   group_by(phylum) %>%
   mutate(median_wos = median(Total_wos, na.rm=T),
          n = n()) %>%
   ungroup() %>%
   arrange(desc(kingdom),phylum) %>%
   mutate(phylum = factor(phylum, levels = unique(.$phylum))) %>%
   ggplot(aes(x = Total_wos, y = phylum,
              fill = kingdom, color = kingdom)) +
   geom_point(position = position_jitter(width = 0.35), size = 1, alpha = 0.3) +
   geom_boxplot(width = .8, outlier.shape = NA, alpha = 0.2, col = "grey20") +
   labs(x = "N° papers in the Web of Science [log-scaled axis]", y = NULL) +
   
   # ggimage::geom_phylopic(aes(x = 1, y = phylum, image = image),
   #                        size = .2, color = "grey20") +
    scale_color_manual(values = custom_color)+
    scale_fill_manual(values = custom_color)+
    scale_x_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 10, 100, 1000, 10000)) + 
    theme_classic() +
    theme(legend.position = "none", axis.text.y = element_text(size = 12)))

plotS2b <- db %>% group_by(phylum) %>%
    mutate(median_wiki = median(total_wiki_pgviews, na.rm=T),
           n = n()) %>%
    ungroup() %>%
    arrange(desc(kingdom) ,phylum) %>% 
    mutate(phylum = factor(phylum, levels = unique(.$phylum))) %>% 
    ggplot(aes(x = total_wiki_pgviews, y = phylum, fill = kingdom, color = kingdom)) +
    geom_point(position = position_jitter(width = 0.35), size = 1, alpha = 0.3) +
    geom_boxplot(width = .8, outlier.shape = NA, alpha = 0.2, col = "grey20") +
    labs(x = "N° views in Wikipedia [log-scaled axis]", y = NULL) +
    scale_color_manual("",values = custom_color)+
    scale_fill_manual("",values = custom_color)+
    scale_x_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 1000, 100000, 50000000)) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 12))

pdf(file = "Figures/Figure_S2.pdf", width = 12, height = 7)
ggpubr::ggarrange(plotS2a,
                  plotS2b,
                  common.legend = FALSE,
                  hjust = -0.5,
                  align = "h",
                  labels = c("A", "B"),
                  ncol=2, nrow=1) 
dev.off()

# Figure S3 ----------------------------------------------------------------

pdf(file = "Figures/Figure_S3.pdf", width = 7.5, height = 4.5)
M0.forest_plot
dev.off()
