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

# Questions: VPA?

# Loading R packages ------------------------------------------------------

library("Amelia")
library("dplyr")
library("glmmTMB")
library("sjPlot")
library("tidylog")
library("tidyverse") 

# Functions ---------------------------------------------------------------

source("Functions/Functions.r") #custom functions

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
# wiki_mean_month_pgviews has the best distribution

# Correlations 
db %>% dplyr::select(Total_wos, total_wiki_pgviews, wiki_langs, wiki_mean_month_pgviews) %>% 
  GGally::ggpairs() + theme_classic()

## Variables X ##
# Distribution 
ggplot(db, aes(x = uniqueness_family)) + geom_dotplot(binaxis='x', stackdir='center',dotsize=0.1) + theme_bw()
ggplot(db, aes(x = uniqueness_genus)) + geom_dotplot(binaxis='x', stackdir='center',dotsize=0.1) + theme_bw()
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
nlevels(db$order) #Devilish!
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

# Homogenize distribution
db <- db %>% 
  dplyr::mutate(log_uniqueness_family = log(uniqueness_family+1),
                log_uniqueness_genus = log(uniqueness_genus+1),
                log_range_size = log(range_size+1),
                log_size_avg = log(size_avg+1),
                log_distance_hu = log(mean_divergence_time_Mya+1))

# scaling size and range size by group
scaled_size <- c()
for(i in 1:nlevels(db$phylum)) {
  db2 <- db[db$phylum == levels(db$phylum)[i], ]
  scaled_size <- append(scaled_size, scale(db2$log_size_avg)[,1]) }

scaled_range_size <- c()
for(i in 1:nlevels(db$phylum)) {
  db2 <- db[db$phylum == levels(db$phylum)[i], ]
  scaled_range_size <- append(scaled_range_size, scale(db2$log_range_size)[,1]) }

db <- data.frame(db, scaled_size, scaled_range_size) ; rm(db2, scaled_size, scaled_range_size, i)

db$scaled_uniqueness_family   <- scale(db$log_uniqueness_family, center = TRUE, scale = TRUE)
db$scaled_uniqueness_genus    <- scale(db$log_uniqueness_genus, center = TRUE, scale = TRUE)
db$scaled_log_distance_human  <- scale(db$log_distance_hu, center = TRUE, scale = TRUE)



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
                               # x = centroid_long,
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
nlevels(dbWOS$phylum)
nlevels(dbWOS$class)
nlevels(dbWOS$order)
nlevels(dbWOS$family)

# dbWOS$cor.str <- numFactor(dbWOS$x,dbWOS$y)
# dbWOS$group   <- factor(rep(1, nrow(dbWOS)))

model.formula.WOS <- as.formula(paste0("WOS ~",
                                   paste(colnames(dbWOS)[7:ncol(dbWOS)], collapse = " + "),
                                   #"+ exp(cor.str + 0 | group)",
                                   "+ (1 | phylum) + (1 | class) + (1 | order) + (1 | family)"))

# Fit the model -----------------------------------------------------------

# First model
M1 <- glmmTMB::glmmTMB(model.formula.WOS,
                       family = poisson, 
                       data = dbWOS)
                       
# Is the model overdispersed?
performance::check_overdispersion(M1) #yes!

# Negative binomial
M1_nbinom <- update(M1, family = nbinom1)

performance::check_zeroinflation(M1_nbinom) #yes!

# # Check for zero-inflation
# 
# Observed zeros: 1032
# Predicted zeros: 21
# Ratio: 0.02

# Zero-inflated
M1_zinb1 <- glmmTMB::glmmTMB(model.formula.WOS,
                                 ziformula = ~ 1,
                                 data = dbWOS,
                                 family = nbinom1)

M1_zinb2 <- update(M1_zinb1, family = nbinom2)

M1_hurdle <- update(M1_zinb2,
                       ziformula = ~ .,
                       family = truncated_nbinom2) # Hurdle model

# Comparing the models
AIC(M1, M1_nbinom, M1_zinb1, M1_zinb2, M1_hurdle)

# Validation
performance::check_model(M1_hurdle)    
performance::check_collinearity(M1_hurdle)

# 
# # Checking spatial autocorrelation 
# #
# # # Is there spatial autocorrelation?
# res   <- DHARMa::simulateResiduals(M2)
# 
# DHARMa::testSpatialAutocorrelation(res, 
#                                    jitter(dbWOS$x,0.0000000001),
#                                    dbWOS$y, plot = FALSE)

# Summary table
(par.M1 <- parameters::parameters(M1_hurdle))

table.M1 <- par.M1 %>% dplyr::select(Parameter,
                         Component,
                         Effects,
                         Beta = Coefficient,
                         SE,
                         CI_low,
                         CI_high,
                         p) %>% 
                         data.frame() %>% 
                         mutate_if(is.numeric, ~ round(.,3)) 

table.M1 <- table.M1[table.M1$Effects == "fixed",] %>% 
            dplyr::select(-c(Effects)) %>% 
            na.omit()

table.M1$Parameter <- as.factor(as.character(table.M1$Parameter))
table.M1$Component <- as.factor(as.character(table.M1$Component))


levels(table.M1$Component) <- c("Conditional", "Zero-inflated")

levels(table.M1$Parameter) <- c("Intercept",
                                "Color blue",
                                "Color red",
                                "Colorful",
                                "Common name [yes]",
                                "Domain [aquatic]",
                                "Domain [terrestrial]",
                                "Harmful to humans [yes]",
                                "Human use [yes]",
                                "IUCN [non-endangered]",
                                "IUCN [endangered]",
                                "Phylogenetic distance to humans",
                                "Range size",
                                "Organism size",
                                "Family uniqueness (N° species)",
                                "Genus uniqueness (N° species)")

# R^2
(M1.R2 <- my.r2(M1_hurdle))

# A general look
sjPlot::plot_model(M1_hurdle, sort.est = FALSE, se = TRUE,
                   vline.color ="grey70",
                   show.values = TRUE, value.offset = .3) + theme_bw()

sign <- ifelse(table.M1$p > 0.05, "", ifelse(table.M1$p>0.01," *", " **")) #Significance
col_p <- ifelse(table.M1$p > 0.05, "grey5", ifelse(table.M1$Beta>0,"orange", "blue") )

table.M1 %>% ggplot2::ggplot(aes(Parameter, Beta)) + facet_wrap(. ~ Component, nrow = 1, ncol = 2) +               
                geom_hline(lty = 3, size = 0.5, col = "grey50", yintercept = 0) +
                geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0, col = "grey10")+
                geom_point(size = 2, pch = 21, col = "grey10", fill = "grey20") +
                geom_text(
                label = paste0(round(table.M1$Beta, 3), sign, sep = "  "), vjust = - 1, size = 2) +
                labs(title = "Scientific interest [N° papers in the Web of Science]",
                     subtitle = paste0("[Sample size = ", nrow(dbWOS) ," observations]"),
                     y = expression(paste("Estimated beta" %+-% "95% Confidence interval")),
                       x = NULL) +
                theme_classic() + 
                coord_flip()  +
                geom_text(data = data.frame(x = 2, y = 2.4, Component = "Zero-inflated", 
                          label = paste0("R^2 ==",round(as.numeric(M1.R2[2]),2))), 
                          aes(x = x, y = y, label = label), 
                          size = 3, parse = TRUE)

            










# A general look
dbWOS %>% ggplot2::ggplot(aes(x = scaled_log_distance_human, y = WOS)) +
  geom_point(col = "grey10", fill = "grey30", size = 3, shape = 21, alpha = 0.3)+
  geom_smooth(method = "gam",  se = TRUE, 
              formula = y ~ s(x),
              method.args = list(family = poisson)) +
  theme_classic() 

# A general look
dbWOS %>% ggplot2::ggplot(aes(x = log_size_avg, y = WOS)) +
  geom_point(col = "grey10", fill = "grey30", size = 3, shape = 21, alpha = 0.3)+
  geom_smooth(method = "gam",  se = TRUE, 
              formula = y ~ s(x),
              method.args = list(family = poisson)) +
  theme_classic() 

# A general look
dbWOS %>% ggplot2::ggplot(aes(x = log_range_size, y = WOS)) +
  geom_point(col = "grey10", fill = "grey30", size = 3, shape = 21, alpha = 0.3)+
  geom_smooth(method = "gam",  se = TRUE, 
              formula = y ~ s(x),
              method.args = list(family = poisson)) +
  theme_classic() 

# A general look
dbWOS %>% ggplot2::ggplot(aes(x = centroid_lat, y = WOS)) +
  geom_point(col = "grey10", fill = "grey30", size = 3, shape = 21, alpha = 0.3)+
  geom_smooth(method = "gam",  se = TRUE, 
              formula = y ~ s(x),
              method.args = list(family = poisson)) +
  theme_classic() 

colnames(db)
######################
## Popular interets ##
######################

dbWIKI2 <- db %>% dplyr::select(wiki = wiki_mean_month_pgviews,
                                kingdom,
                                phylum,
                                class,
                                order,
                                family,
                                y = centroid_lat,
                                x = centroid_long,
                                harmful_to_human, #cultural
                                human_use, #cultural
                                common_name, #cultural
                                colorful,  #aesthethic
                                color_blu, #aesthethic
                                color_red, #aesthethic
                                IUCN_rec, #rarity
                                domain_rec,
                                starts_with("scaled_")) #aesthethic / rarity
# Scaling variables 
# dbWIKI2 <- dbWIKI2 %>% dplyr::select(-c(wiki)) %>%
#   mutate_if(is.numeric, ~ scale(., center = TRUE, scale = TRUE)) %>% 
#   cbind(wiki = dbWIKI2$wiki, .)

# Missing data
Amelia::missmap(dbWIKI2)
dbWIKI <- na.omit(dbWIKI2)

model.formula <- as.formula(paste0("wiki ~",
                                   paste(colnames(dbWIKI)[7:ncol(dbWIKI)], collapse = " + "),
                                   "+ (1 | phylum) + (1 | class) + (1 | order) + (1 | family)"))


(M3 <- glmmTMB::glmmTMB(model.formula,
                        family = poisson, data = dbWIKI))

# Is the model overdispersed?
performance::check_overdispersion(M3)  

(M4 <- glmmTMB::glmmTMB(model.formula,
                        family = nbinom1, data = dbWIKI)) 
#control=glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS"))

performance::check_model(M4)                            

summary(M4)

my.r2(M4)

sjPlot::plot_model(M4, sort.est = TRUE, se = TRUE,
                   vline.color ="grey70",
                   show.values = TRUE, value.offset = .3) + theme_bw()

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

pdf(file = "./Figures/Figure 1.pdf", width = 7, height = 5)

(plot_cor <- gridExtra::grid.arrange(density_WOS, blankPlot, cor_plot, density_WIKI, 
                          ncol=2, nrow=2, widths=c(3, 1.7), heights=c(1.7, 3)) )

dev.off()
