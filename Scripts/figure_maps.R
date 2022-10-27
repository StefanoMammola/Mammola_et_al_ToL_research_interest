# Load libraries
library(ggplot2)
library(tidyverse)

# Load world map
world <- map_data("world")

# Load dataset ####
sp <- read.csv("/media/tcanc/PHILIPS UFD/Articles_Tom/Tree_of_life_Mammola_()/Dataset/SampleTREE_191022.csv", sep = "\t")
head(sp); tail(sp); str(sp)

# filter
sp.1 <- sp[ ,c("kingdom", "phylum", "class", "order", "family", "name", "centroid_lat", "centroid_long")]
head(sp.1)

write.csv(sp.1, "/media/tcanc/PHILIPS UFD/Articles_Tom/Tree_of_life_Mammola_()/Dataset/SampleTREE_191022_NoNA.csv", row.names = FALSE)

# Check NA in long and lat
length(which(is.na(sp.1$centroid_lat)))

# Remove rows with NA
sp.1 <- sp.1[complete.cases(sp.1), ]
length(which(is.na(sp.1$centroid_lat))) # <- OK

# Convert to numeric
sp.1$centroid_long <- gsub(",", ".", sp.1$centroid_long)
sp.1$centroid_long <- as.numeric(sp.1$centroid_long)
sp.1$centroid_lat <- gsub(",", ".", sp.1$centroid_lat)
sp.1$centroid_lat <- as.numeric(sp.1$centroid_lat)


ggplot() +
  geom_map(data = world, map = world,
    aes(long, lat, map_id = region),
    color = "black", fill = "lightgray", size = 0.1) +
  geom_point(data = sp.1,
    aes(centroid_long, centroid_lat), alpha = 0.7, color = "purple", size = .5) +
  theme_bw() + 
  labs(title="Total occurrences")

# Animalia
sp.1.animalia <- sp.1[sp.1$kingdom == "Animalia", ]
unique(sp.1.animalia$kingdom)

ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.1) +
  geom_point(data = sp.1.animalia,
             aes(centroid_long, centroid_lat), alpha = 0.7, color = "#FFAE42", size = .3) +
  theme_bw() + 
  labs(title="Animalia")

# Plantae
sp.1.plantae <- sp.1[sp.1$kingdom == "Plantae", ]
unique(sp.1.plantae$kingdom)

ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.1) +
  geom_point(data = sp.1.plantae,
             aes(centroid_long, centroid_lat), alpha = 0.7, color = "green", size = .5) +
  theme_bw() + 
  labs(title="Plantae")

# Fungi
sp.1.fungi <- sp.1[sp.1$kingdom == "Fungi", ]
unique(sp.1.fungi$kingdom)

ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.1) +
  geom_point(data = sp.1.fungi,
             aes(centroid_long, centroid_lat), alpha = 0.7, color = "purple", size = .5) +
  theme_bw() + 
  labs(title="Fungi")

