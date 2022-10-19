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
library("glmmTMB")
library("tidyverse") 

# Loading the database  ---------------------------------------------------

db <- read.csv2(file = "./Data/SampleTREE_181022.csv", sep = '\t', dec = ',', header = TRUE, as.is = FALSE)
db <- db[db$FOSSIL != "1",] ; db <- droplevels(db)

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
ggplot(db, aes(x = Total_wos)) + geom_histogram(bins = 100) + theme_bw()

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

# Correlations 
db %>% dplyr::select(uniqueness_family, uniqueness_genus, range_size, size_avg, centroid_lat) %>% 
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

# Homogenize distribution
db <- db %>% 
  dplyr::mutate(log_uniqueness_family = log(uniqueness_family+1),
                log_uniqueness_genus = log(uniqueness_genus+1),
                log_range_size = log(range_size+1),
                log_size_avg = log(size_avg+1))

# scaling size??
scale_size <- c()

for(i in 1:nlevels(db$phylum)) {
  db2 <- db[db$phylum == levels(db$phylum)[i], ]
  scale_size <- append(scale_size, scale(db2$size_avg)[,1]) }

db <- data.frame(db,scale_size)

# Assembling a final database ---------------------------------------------

#######################
## Research interets ##
#######################

dbWOS2 <- db %>% dplyr::select(WOS = Total_wos,
                            #wiki = wiki_mean_month_pgviews,
                            kingdom,
                            phylum,
                            class,
                            order,
                            family,
                            harmful_to_human,
                            human_use,
                            common_name,
                            colorful,
                            color_blu,
                            color_red,
                            IUCN_rec,
                            centroid_lat,
                            starts_with("log_"))

# Scaling variables 
dbWOS2 <- dbWOS2 %>% dplyr::select(-c(WOS)) %>%
  mutate_if(is.numeric, ~ scale(., center = TRUE, scale = TRUE)) %>% 
  cbind(WOS = dbWOS2$WOS, .)

# Missing data
Amelia::missmap(dbWOS2)
dbWOS <- na.omit(dbWOS2)

# Outliers

dbWOS <- dbWOS[dbWOS$WOS < 2000,]

# Setting formula ---------------------------------------------------------

nlevels(dbWOS$phylum)
nlevels(dbWOS$class)
nlevels(dbWOS$order)
nlevels(dbWOS$family)

model.formula <- as.formula(paste0("WOS ~",
                                   paste(colnames(dbWOS)[7:ncol(dbWOS)], collapse = " + "),
                                   "+ (1 | phylum) + (1 | class) + (1 | order) + (1 | family)"))



(M1 <- glmmTMB::glmmTMB(model.formula,
                        family=poisson, data=dbWOS))

# Is the model overdispersed?
performance::check_overdispersion(M1)  

(M2 <- glmmTMB::glmmTMB(model.formula,
               family = nbinom1, data = dbWOS)) #control=glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS"))

performance::check_model(M2)                            

summary(M2)

# A general look
dbWOS %>% ggplot2::ggplot(aes(x = centroid_lat, y = WOS)) +
  geom_point(col = "grey10", fill = "grey30", size = 3, shape = 21, alpha = 0.3)+
  geom_smooth(method = "gam",  se = TRUE, 
              formula = y ~ s(x),
              method.args = list(family = poisson)) +
  theme_classic() #Vague problem at the tail of the distribution, but overall OK!

colnames(dbWOS)

# A general look
dbWOS %>% ggplot2::ggplot(aes(x = log_size_avg, y = WOS)) +
  geom_point(col = "grey10", fill = "grey30", size = 3, shape = 21, alpha = 0.3)+
  geom_smooth(method = "gam",  se = TRUE, 
              formula = y ~ s(x),
              method.args = list(family = poisson)) +
  theme_classic() #Vague problem at the tail of the distribution, but overall OK!


# A general look
dbWOS %>% ggplot2::ggplot(aes(x = log_range_size, y = WOS)) +
  geom_point(col = "grey10", fill = "grey30", size = 3, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm",  se = TRUE, 
              formula = y ~ x,
              method.args = list(family = poisson)) +
  theme_classic() #Vague problem at the tail of the distribution, but overall OK!

