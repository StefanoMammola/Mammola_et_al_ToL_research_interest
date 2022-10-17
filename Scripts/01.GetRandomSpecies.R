## ------------------------------------------------------------------------
## 'Species popularity and research interests across the Tree of Life' 
## ------------------------------------------------------------------------

##################
# Created by Stefano Mammola
# Updated 17/03/2022
##################

## Started: 03/10/2020, Helsinki: core loop & test sample.
## Final sampling: 15/10/2020, Helsinki

## R script to extract random species from GBIF backbone taxonomy

# ____________________________________________________

# Random stratified sampling of the eukaryote Tree of Life [Animalia, Fungi (Agaricomycetes), and Plantae (excluding unicellular Algae)].
# Based on the backbone taxonomy of GBIF: https://www.gbif.org/dataset/d7dddbf4-2cf0-4f39-9b2a-bb099caae36c

# Conditions:

# 1) Sample at the level of species within Orders.

# 2) For each Order, a fraction of 0.002 species are sampled.

# 3) To avoid having an excessively uneven number of species among Orders, thresholds are set:
#   - If the number of species in an Order is comprised between 10001 and 50000, 20 species are arbitrarily sampled.
#   - If the number of species in an Order is comprised between 50001 and 100000, 40 species are arbitrarily sampled.
#   - If the number of species in an Order is >100000, 60 species are arbitrarily sampled.

# 4) To incorporate a broader sample of tetrapods so as to reflect the typical knowledge bias ("Institutional vertebratism" sensu Leather), two conditions are further set:
#   - For each tetrapod Order, 20 species are arbitrarily sampled.
#   - However, if the number of species in a tetrapod Order is < 10, only 1 species is sampled.

# 5) Subspecies and varieties are excluded.

# 6) Fossils are excluded. 
#    Note, however, that a few fossil not properly labelled as "FOSSIL_SPECIMEN" in GBIF are not excluded by the loop. 
#    These need to be manually resampled.

# ____________________________________________________


# Loading R packages ------------------------------------------------------

library("data.table") #<---- to read tsv
library("httr")       #<---- to set internet settings
library("rgbif")      #<---- to interact with GBIF
library("sp")         #<---- miscellaneous functions
library("xlsx")       #<---- to save .xlsx files

# Loading the database  ---------------------------------------------------

# Backbone species database from GBIF
data <- as.data.frame(fread("taxon.tsv")) #<---- read GBIF backbone table 
# Note that this database is not provided due to large size (>2 Giga). 
# To be extracted from GBIF: https://www.gbif.org/dataset/d7dddbf4-2cf0-4f39-9b2a-bb099caae36c

head(data,3) #checking the first 3 rows
str(data)    #checking the structure of the database

# Preparing data ----------------------------------------------------------

# Selecting only accepted names
data$taxonomicStatus <- as.factor(data$taxonomicStatus)
data <- data[data$taxonomicStatus == "accepted", ]

# Removing subspecies & varieties
data <- data[data$taxonRank != "subspecies", ]
data <- data[data$taxonRank != "variety", ]

# Subsetting the database
db <- data.frame(taxonID = data$taxonID, 
                 kingdom = data$kingdom, 
                 phylum = data$phylum, 
                 class  = data$class, 
                 order = data$order,
                 family = data$family,
                 genus = data$genericName,
                 species = data$specificEpithet,
                 name = data$canonicalName,
                 author = data$scientificNameAuthorship)

# Cleaning the database and selecting only the eukaryotic Tree of Life ------------

db <- db[db$kingdom == "Animalia" | db$kingdom == "Fungi" & db$class == "Agaricomycetes" | db$kingdom == "Plantae", ] #take animals, plants, and Macroscopic fungi

# Removing divisions with mostly unicellular algae.
db <- db[db$phylum != "Chlorophyta",]
db <- db[db$phylum != "Glaucophyta",] 
db <- db[db$phylum != "Charophyta",] 

# Removing obvious fossil groups
db <- db[db$class != "Trilobita",] 

db <- db[db$order != "Protocoleoptera",]
db <- db[db$order != "Thecodontia",]
db <- db[db$order != "Bolosauria",]
db <- db[db$order != "Pterosauria",]
db <- db[db$order != "Dinosauria",]
db <- db[db$order != "Rhynchosauria",]
db <- db[db$order != "Pareiasaurida",]
db <- db[db$order != "Ammonoidea",]
db <- db[db$order != "Miomoptera",]

db <- db[db$family != "Megalonychidae",]
db <- db[db$family != "Megatheriidae",] 
db <- db[db$family != "Mylodontidae",]
db <- db[db$family != "Scelidotheriidae",]
db <- db[db$family != "Orophodontidae",]
db <- db[db$family != "Nothrotheriidae",]
db <- db[db$family != "Megalocnidae",]
db <- db[db$family != "Palaeotheriidae",]

# ... and selecting only living fossils

ginko <- db[db$name ==  "Ginkgo biloba",] ; db <- db[db$class !=  "Ginkgoopsida",] 
db <- rbind(db,ginko) ; rm(ginko)

latimeria <- db[db$name ==  "Latimeria chalumnae",] ; db <- db[db$family !=  "Latimeriidae",] 
db <- rbind(db,latimeria) ; rm(latimeria)

notiothauma <- db[db$name ==  "Notiothauma reedi",] ; db <- db[db$family !=  "Eomeropidae",] 
db <- rbind(db,notiothauma) ; rm(notiothauma)

jurodes <- db[db$name == "Jurodes ignoramus",] ; db <- db[db$family !=  "Jurodidae",] 
db <- rbind(db,jurodes) ; rm(jurodes)

syntexis <- db[db$name == "Syntexix libocedrii",] ; db <- db[db$family != "Anaxyelidae",]
db <- rbind(db,syntexis) ; rm(syntexis)

eleph <- db[db$name == "Elephas maximus" | db$name == "Loxodonta africana" | db$name == "Loxodonta cyclotis",] ;  db <- db[db$order != "Proboscidea" ,]

# Removing missing levels 
db <- db[db$phylum  != "",]
db <- db[db$class   != "",]
db <- db[db$order   != "",]
db <- db[db$family  != "",]
db <- db[db$genus   != "",]
db <- db[db$species != "",]
db <- db[db$name    != "",]
db <- db[db$author  != "",]

# Dropping unused levels of the factors
db$kingdom <- droplevels(db)$kingdom
db$phylum  <- droplevels(db)$phylum
db$class   <- droplevels(db)$class
db$order   <- droplevels(db)$order
db$family  <- droplevels(db)$family
db$genus   <- droplevels(db)$genus
db$species <- droplevels(db)$species

rm(data) #cleaning

# Sampling a random number of species -------------------------------------

#Run from here (Warning: it takes > 1 hour):
n_kingdom  <- length(unique(db$kingdom)) #How many kingdoms
SampleTree <- data.frame(db[1,], uniqueness_family = NA , uniqueness_genus = NA) #db to fill
fossil_sp <- c() #vector of fossil species to omit

# ________________________________ Looping the Kingdom __________________________________
for(i in 1 : n_kingdom ){
  
  message(paste("———————————— Selecting the kingdom ", as.character(unique(db$kingdom)[i]), " ———————————————————————— " ))
  db2 <- db[db$kingdom == as.character(unique(db$kingdom)[i]),  ]
  
  n_phyla <-  length(unique(db2$phylum)) #How many Phyla
  
  # ________________________________ Looping the phylum __________________________________
  for(j in 1 : n_phyla) {  
    
    message(paste("_ Selecting the phylum ", as.character(j) , " out of " , n_phyla ,": " ,as.character(unique(db2$phylum)[j])))
    db3 <- db2[db2$phylum == as.character(unique(db2$phylum)[j]),]
    
    n_class <- length(unique(db3$class)) #How many classes
    
    # ________________________________Looping the class_________________________________
    for(k in 1 : n_class) { 
      
      message(paste("_______ Selecting the class ", as.character(k) , " out of " , n_class ,": " ,as.character(unique(db3$class)[k])))
      db4 <- db3[db3$class == as.character(unique(db3$class)[k]),]
      
      n_ord <- length(unique(db4$order)) #How many orders
      
      # ________________________________Looping the Order________________________________
       for(l in 1 : n_ord) { 
        
        message(paste("___________ Selecting the order ", as.character(l) , " out of " , n_ord ,": " ,as.character(unique(db4$order)[l])))
        db5 <- db4[db4$order == as.character(unique(db4$order)[l]),]
        
        n_sp <- nrow(db5)
        
        ### ____Estimating number of species to sample____
        if(n_sp > 100000) { frac_to_sample = 60 } else if(n_sp > 50000) { frac_to_sample = 40 } else if(n_sp > 10000) { frac_to_sample = 20 }  else {  frac_to_sample = round(n_sp * 0.002, 1) }
        
        # Setting number of species to 20 for tetrapods to reflect the typical taxonomic bias
        if(unique(db4$class) == "Aves" )     { frac_to_sample = 20 }
        if(unique(db4$class) == "Amphibia" ) { frac_to_sample = 20 }
        if(unique(db4$class) == "Mammalia" ) { frac_to_sample = 20 }
        if(unique(db4$class) == "Reptilia" ) { frac_to_sample = 20 }
        if(frac_to_sample > n_sp)            { frac_to_sample = 1  } #but not for groups with less than 10 species
        
        # Always sample at least one species
        if(frac_to_sample == 0) { frac_to_sample = 1 }
        
        ### ____Sampling the Species____
        row_sequence  <- seq(from = 1, to = n_sp, by = 1)
        random_sample <- sample(row_sequence)[1 : frac_to_sample] 
        
        db6 <- db5[ random_sample, ] 
        
        ## ____Calculating taxonomic uniqueness and range size____
        uniqueness_genus  <- c() ; uniqueness_family <- c() 
        
        print("... Calculating taxonomic uniqueness ... ")
        
        for(n in 1: nrow(db6)) { 
          
          # ____taxonomic uniqueness____
          uniqueness_genus  <- append(uniqueness_genus,  nrow(db5[db5$genus ==  as.character( db6[n,]$genus   ), ]) )
          uniqueness_family <- append(uniqueness_family, nrow(db5[db5$family ==  as.character( db6[n,]$family ), ]) )
          
        }
        
        ## ____Storing the result____
        SampledTree <- data.frame(db6, 
                                  uniqueness_family, uniqueness_genus)
        
        SampleTree <- rbind(SampleTree, SampledTree)
        
      }
    }
  }
} #end

# Store results
SampleTree <- SampleTree[-1,] #removing the first row
SampleTree <- SampleTree[ !(SampleTree$name) %in% c(fossil_sp),] #removing fossils 

# Save
write.csv2(SampleTre, "./Data/SampleTree.csv", row.names = F)

#end