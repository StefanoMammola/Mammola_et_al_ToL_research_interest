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

library("dplyr")
library("tidyverse") 

# Loading the database  ---------------------------------------------------

db <- read.csv2(file = "./Data/SampleTREE_171022.csv", sep = '\t', dec = ',', header = TRUE, as.is = FALSE)

db <- db[db$FOSSIL != "1",]

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

# Variable Y -------------------------------------------------------

# Correlations

db %>% dplyr::select(Total_wos, total_wiki_pgviews, wiki_langs, wiki_mean_month_pgviews) %>% 
  GGally::ggpairs() + theme_classic()

# Y (distribution)
ggplot(db, aes(x = total_wiki_pgviews)) + geom_dotplot(binaxis='x', stackdir='center',dotsize=0.1) + theme_bw()
ggplot(db, aes(x = total_wiki_pgviews)) + geom_histogram(bins = 100) + theme_bw()
ggplot(db, aes(x = log(total_wiki_pgviews+1))) + geom_histogram() + theme_bw()

ggplot(db, aes(x = Total_wos)) + geom_dotplot(binaxis='x', stackdir='center',dotsize=0.1) + theme_bw()
ggplot(db, aes(x = Total_wos)) + geom_histogram(bins = 100) + theme_bw()
ggplot(db, aes(x = log(Total_wos+1))) + geom_histogram() + theme_bw()

# 
ggplot(db, aes(x = uniqueness_family)) + geom_dotplot(binaxis='x', stackdir='center',dotsize=0.1) + theme_bw()
ggplot(db, aes(x = uniqueness_genus)) + geom_dotplot(binaxis='x', stackdir='center',dotsize=0.1) + theme_bw()
ggplot(db, aes(x = range_size)) + geom_dotplot(binaxis='x', stackdir='center',dotsize=0.1) + theme_bw()


plot(log(db$range_size+1),db$Total_wos)

# scaling size
scale_size <- c()

for(i in 1:nlevels(db$phylum)) {
  db2 <- db[db$phylum == levels(db$phylum)[i], ]
  scale_size <- append(scale_size, scale(db2$size_avg)[,1]) }

db <- data.frame(db,scale_size)
plot(log(db$scale_size+1),db$Total_wos)
plot(db$mean_divergence_time_Mya,db$Total_wos)

plot(db$centroid_lat,db$total_wiki_pgviews)

# Relationships

ggplot(db, aes(x = phylum, y = log(Total_wos+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = phylum, y = log(total_wiki_pgviews+1))) + geom_boxplot() + theme_bw() + coord_flip()

ggplot(db, aes(x = harmful_to_human, y = log(Total_wos+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = human_use, y = log(Total_wos+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = model_organism, y = log(Total_wos+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = common_name, y = log(Total_wos+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = colorful, y = log(Total_wos+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = color_blu, y = log(Total_wos+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = color_red, y = log(Total_wos+1))) + geom_boxplot() + theme_bw() + coord_flip()

ggplot(db, aes(x = harmful_to_human, y = log(total_wiki_pgviews+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = human_use, y = log(total_wiki_pgviews+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = model_organism, y = log(total_wiki_pgviews+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = common_name, y = log(total_wiki_pgviews+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = colorful, y = log(total_wiki_pgviews+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = color_blu, y = log(total_wiki_pgviews+1))) + geom_boxplot() + theme_bw() + coord_flip()
ggplot(db, aes(x = color_red, y = log(total_wiki_pgviews+1))) + geom_boxplot() + theme_bw() + coord_flip()

db %>% dplyr::select(size_avg, mean_divergence_time_Mya, n_occurrences, range_size, range_size_AOO.y, range_size_MCP.y) %>% 
  GGally::ggpairs() + theme_classic()


