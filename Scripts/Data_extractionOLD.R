## ------------------------------------------------------------------------
## 'Species popularity and research interests across the Tree of Life' 
## ------------------------------------------------------------------------

## Started: 03 October 2020, Helsinki: core loop & test sample.
## Updated: 12 October 2020, Helsinki: created extract.IUCN & extract.GBIF functions.
## Final sampling: 15 October 2020, Helsinki

# Loading R packages ------------------------------------------------------

library("data.table") #<---- to read tsv
library("geodist")
library("httr")       #<---- to set internet settings
library("raster")     #<---- for operation with rasters
library("red")        #<---- for range size (MCP and AOO)
library("rredlist")   #<---- to interact with IUCN list
library("rgbif")      #<---- to interact with GBIF
library("sp")         #<---- miscellaneous functions
library("xlsx")       #<---- to save .xlsx files

# Useful functions --------------------------------------------------------

## Function to extract IUCN data avoiding errors in internet searches:
extract.IUCN <- function (sp_name, api_IUCN) {
  
  # check input data
  if(class(sp_name) != "character")
    stop("Species name should be a character")
  
  if(class(api_IUCN) != "character")
    stop("Api key should be a character")
  
  tryCatch(
    
    expr = {
      
      #check global database
      if(is.data.frame(rl_history(name = sp_name , key = api_IUCN, region = "global")$result) == TRUE){ 
        common_name  <- ifelse(is.data.frame(rl_common_names(name = sp_name, key = api_IUCN, region = "global")$result) == TRUE, 1, 0) 
        IUCN         <- rl_history(name = sp_name, key = api_IUCN, region = "global")$result[1,3] 
        IUCN_habitat <- ifelse(is.data.frame(rl_habitats(name = sp_name,  key = api_IUCN, region = "global")$result) == TRUE, paste(unique(rl_habitats(name = sp_name,  key = api_IUCN, region = "global")$result[,2]),collapse=" ; "), NA)
      } 
      #check also in the european base:
      else if(is.data.frame(rl_history(name = sp_name , key = api_IUCN, region = "europe")$result) == TRUE){ 
        common_name  <- ifelse(is.data.frame(rl_common_names(name = sp_name, key = api_IUCN, region = "europe")$result) == TRUE, 1, 0)
        IUCN         <- rl_history(name = sp_name, key = api_IUCN, region = "europe")$result[1,3]
        IUCN_habitat <- ifelse(is.data.frame(rl_habitats(name = sp_name,  key = api_IUCN, region = "europe")$result) == TRUE, paste(unique(rl_habitats(name = sp_name,  key = api_IUCN, region = "europe")$result[,2]),collapse=" ; "), NA)
      } 
      #if not evaluated:
      else { 
        common_name  <- NA
        IUCN         <- "Not Evaluated"
        IUCN_habitat <- NA
      }
      
      results_IUCN <- list(common_name, IUCN, IUCN_habitat)
      
    },
    
    #if Bad Gateway (HTTP 502) error happens, recursively call the same function
    error = function(e){ 
      
      print(" ... Bad Gateway (HTTP 502) error occurred ... ")
      Sys.sleep(10) #pause for 10 seconds
      print(" ... recalculating ...")
      
      results_IUCN2 <- extract.IUCN(sp_name = sp_name, api_IUCN = api_IUCN)
      
      results_IUCN <- list(common_name = results_IUCN2[[1]], IUCN = results_IUCN2[[2]], IUCN_habitat = results_IUCN2[[3]])
      
    }
    
  )
  
  return(results_IUCN)
  
}

## Function to extract GBIF data avoiding errors in internet searches:
extract.GBIF <- function (taxonID) {
  
  if(class(taxonID) != "integer")
    stop("TaxonID should be an integer")

  tryCatch(
  
  expr = {
    
    GBIF <- rgbif::occ_search(taxonKey = taxonID) #extracting GBIF coords
    
  },
  #if search fails, recursively call the same function...
  error = function(e){ 
    
    print(" ... Internet search error occurred ... ")
    Sys.sleep(10) #pause for 10 seconds
    print(" ... Trying again ... ")
    
    GBIF <- extract.GBIF(taxonID = taxonID)  #extracting GBIF coords
    
  }
)
  
 return(GBIF) 
  
}

# Setting working directory -----------------------------------------------

setwd("/Users/stefanomammola/Desktop/SAMPLING THE THREE OF LIFE/")  #<---- change with your directory

# Loading elevation data --------------------------------------------------

# Global elevation data at a resolution of 5 min.
# Source: WorldClim (https://www.worldclim.org/)
# Fick, S.E. and R.J. Hijmans, 2017. WorldClim 2: new 1km spatial resolution climate surfaces for global land areas. International Journal of Climatology 37 (12): 4302-4315.

elevation <- raster::raster("wc2.1_5m_elev.tif")
#plot(elevation)

# Loading the database  ---------------------------------------------------

# Backbone species database from GBIF
data <- as.data.frame(fread("taxon.tsv")) #<---- read GBIF backbone table

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

# ____________________________________________________

# Random stratified sampling of the eukaryote Tree of Life [Animalia, Fungi (Agaricomycetes), and Plantae (excluding unicellular Algae)].

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

# 5) To reduce the number of missing data (NA) in the database, species with less than 3 georeferenced localities in GBIF are resampled.
#    This is because 3 occurrences are needed to construct a minimum convex polygon to estimate range size.
#    However, if all the species in the Order have less than 3 occurrence, 1 species is arbitrarily kept.

# 6) Subspecies and varieties are excluded.

# 7) Fossils are excluded. 
#    Note, however, that a few fossil not properly labelled as "FOSSIL_SPECIMEN" in GBIF are not excluded by the loop. 
#    These need to be manually resampled.

# ____________________________________________________

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
        
        ### ____Checking if sampled species have enough occurrences____
        remove <- NULL
        
        for(m in 1 : nrow(db6)) {
          
          GBIF <- extract.GBIF(taxonID = db6$taxonID[m])
          coord <- unique(na.omit(cbind(GBIF$data$decimalLongitude, GBIF$data$decimalLatitude)))
          
          if(is.null(coord) == TRUE){
            remove <- append(remove,m)} else if(nrow(coord) < 3){
            remove <- append(remove,m)} else{Sys.sleep(0) }
            
        }
        
        if(is.null(remove) == FALSE) { 
          
          db7 <- db6[-remove, ] #species to keep
          if(nrow(db7) == 0){ db7 = db6[1,] } # If all species are removed, re-establish one
          db8 <- db5[-random_sample, ] # database to resample
          db8 <- db8[sample(nrow(db8)), ] #shuffle rows of db8
          
          print("... Resampling species with too few occurrences in GBIF ... ")
          
          while(nrow(db7) < nrow(db6)){ 
            
            for(y in 1:nrow(db8)){ 
              
              GBIF <- extract.GBIF(taxonID = db8$taxonID[y])
              coord <- unique(na.omit(cbind(GBIF$data$decimalLongitude, GBIF$data$decimalLatitude)))
              
              if(is.null(coord) == TRUE){
                Sys.sleep(0)} else if(nrow(coord) < 3){
                Sys.sleep(0)} else {
                db7 <- rbind(db7,db8[y,]) } 
              
              if (nrow(db7) == nrow(db6)) # If the number of species to sample is reached, break
                break 
              
            }
            
            if (y == nrow(db8))  # If loop is over, and not sufficient species are there, break
              break 
            
          } 
          
          db6 <- db7 ; rm(db7,db8)
          
        }
        
        ### ____Checking if sampled species have fossils____
        remove <- NULL
        
        for(p in 1 : nrow(db6)) {
          if(any(occ_search(taxonKey= db6[p,]$taxonID )$data$basisOfRecord == "FOSSIL_SPECIMEN"))
            remove <- append(remove,p) #omit fossil 
        }
        
        if(is.null(remove) == FALSE) { 
          
          db7 <- db6[-remove, ] #species to keep
          if(nrow(db7) == 0){ db7 = db6[1,] ; fossil_sp <- append(fossil_sp, paste(db6[1,]$genus,db6[1,]$species, sep = ' ')) } # If all species are removed, re-establish one
          db8 <- db5[-random_sample, ] # database to resample
          db8 <- db8[sample(nrow(db8)), ] #shuffle rows of db8
          
          print("... Resampling fossil species  ... ")
          
          while(nrow(db7) < nrow(db6)){ 
            
            for(y in 1:nrow(db8)){ 
              
              if(any(occ_search(taxonKey= db8[y,]$taxonID)$data$basisOfRecord != "FOSSIL_SPECIMEN")) { 
                db7 <- rbind(db7,db8[y,])  } #select new random species with at least 5 occurrence
              else {  NULL } # Skip
              
              if (nrow(db7) == nrow(db6)) # If the number of species to sample is reached, break
                break 
              
            }
            
            if (y == nrow(db8))  # If loop is over, and not sufficient species are there, break
              break 
            
          } 
          
          db6 <- db7 ; rm(db7,db8)
          
        }
        
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
}
#end

# Clean
rm(db2, db3, db4, db5, db6, GBIF, coord, SampledTree, i, j, k, l, m, n, p, n_kingdom, n_class, n_ord, n_phyla, n_sp, 
   random_sample, remove, uniqueness_family, uniqueness_genus, frac_to_sample) #clean

# Saving the results
SampleTree <- SampleTree[-1,] #removing the first row
SampleTree <- SampleTree[ !(SampleTree$name) %in% c(fossil_sp),] #removing fossils 

xlsx::write.xlsx(SampleTree,"SampleTREE.xlsx")
write.table(SampleTree,"SampleTREE.txt")

# Extracting GBIF data ----------------------------------------------------

# Read data again ---------------------------------------------------------

SampleTree <- read.csv("/Users/stefanomammola/Desktop/SAMPLING THE THREE OF LIFE/Analysis/sampleTREE.csv", sep='\t', dec=',', header=T, as.is = F)
str(SampleTree)

# Configuration to minimize problem with online searches ...
httr::set_config(config(ssl_verifypeer = 0L))

n_occurrences     <- c() ; range_size        <- c() 
elevation_min     <- c() ; elevation_max     <- c() ; elevation_range   <- c()
higherGeography   <- c() ; verbatimHabitat   <- c()

n_sp <- nrow(SampleTree)
#n_sp <- 1500

#Run from here: 
#(Warning: it takes approx. 1-2 hours)
for(i in 1501 : n_sp) {  
  
  message(paste("_ Selecting species ", as.character(i) , " out of " , as.character(n_sp), sep = ''))
  
# ____range size____
print("... Calculating range size ... ")  

# Extract GBIF data
GBIF <- extract.GBIF(taxonID = SampleTree[i,]$taxonID)

# omit coordinates with missing data
coord <- na.omit(cbind(GBIF$data$decimalLongitude, GBIF$data$decimalLatitude))

if(length(unique(coord)) > 1){
  n_occurrences   <- append(n_occurrences  , nrow(unique(coord))) #number of unique known localities
  centroid        <- c(mean(coord[,1]),mean(coord[,2]))
  range_size      <- append(range_size, mean(geodist::geodist(coord[,1:2],centroid, measure = "geodesic")))
  elev_r          <- range(extract(elevation,coord), na.rm = TRUE)
  elevation_min   <- append(elevation_min, elev_r[1]  )
  elevation_max   <- append(elevation_max, elev_r[2]  )
  elevation_range <- append(elevation_range, elev_r[2]-elev_r[1]  ) #elevation range
} else if (length(unique(coord)) < 2){
  n_occurrences   <- append(n_occurrences, 1)
  range_size      <- append(range_size, 0)
  elevation_min   <- append(elevation_min, NA)
  elevation_max   <- append(elevation_max, NA)
  elevation_range <- append(elevation_range,0)
} 

# ____Geography verbatim____
print("... Extracting geography ... ")  

biogeography <- unique(na.omit(GBIF$data$higherGeography))

if(is.null(biogeography) == TRUE)
  higherGeography <- append(higherGeography, NA)
else if(is.logical(biogeography) == TRUE)
  higherGeography <- append(higherGeography, NA)
else 
  higherGeography <- append(higherGeography, paste(biogeography,collapse=" ; "))

# ____Habitat verbatim____
print("... Extracting habitat ... ")  

habitat <- unique(na.omit(GBIF$data$habitat))

if(is.null(habitat) == TRUE)
  verbatimHabitat <- append(verbatimHabitat, NA)
else if(is.logical(habitat) == TRUE)
  verbatimHabitat <- append(verbatimHabitat, NA)
else
  verbatimHabitat <- append(verbatimHabitat, paste(habitat,collapse=" ; "))

} #end

# Store results
SampleTree2 <- data.frame(SampleTree[1501:n_sp,1:10], n_occurrences, round(range_size,0), 
                          elevation_min, elevation_max, elevation_range, higherGeography, verbatimHabitat)

# Clean
rm(coord, elevation, GBIF, i, biogeography, elev_r, n_occurrences, range_size_MCP, range_size_AOO, 
   elevation_min, elevation_max, elevation_range, higherGeography, verbatimHabitat) #clean

# Overwrite the results
write.csv(SampleTree2,"range_size_1501_3000.csv")

# Read data again ---------------------------------------------------------

SampleTree <- read.csv("/Users/stefanomammola/Desktop/SAMPLING THE THREE OF LIFE/Analysis/sampleTREE.csv", sep='\t', dec='.', header=T, as.is=F)

str(SampleTree)

# Extracting IUCN data ----------------------------------------------------

# Defining an API for IUCN
api_IUCN <- "1178e14bfcb77f2e7b7f0e402a91eb0fceff98e0f4dead2b3c6e9f5d833828d9" #change with your IUCN API key 

# Starting the loop to extract IUCN data
common_name <- c() ; IUCN <- c() ; IUCN_habitat <- c()

n_sp <- nrow(SampleTree)

#Run from here: 
#(Warning: it takes 1-2 hours)
for(i in 1 : n_sp) {  
  
  message(paste("________ Selecting species ", as.character(i) , " out of " , as.character(n_sp), sep = ''))
  
  #Selecting species name
  sp_name <- as.character(SampleTree[i,]$name)
  
  #Extracting IUCN data
  print("... Extracting IUCN data ... ") 
  
  results_IUCN <- extract.IUCN(sp_name = sp_name, api_IUCN = api_IUCN)
  
  # Store results of the run
  common_name  <- append(common_name, results_IUCN[[1]])
  IUCN         <- append(IUCN, results_IUCN[[2]])
  IUCN_habitat <- append(IUCN_habitat, results_IUCN[[3]])

} #end

# Store results
SampleTree <- data.frame(SampleTree, common_name, IUCN, IUCN_habitat)

# Clean
rm(api_IUCN, n_sp, i, common_name, IUCN, IUCN_habitat,sp_name,results_IUCN)

# Overwrite the results
xlsx::write.xlsx(SampleTree,"SampleTREE.xlsx")
write.table(SampleTree,"SampleTREE.txt")

# Test --------------------------------------------------------------------

abc <- read.csv(file = "/Users/stefanomammola/Desktop/SAMPLING THE THREE OF LIFE/SampleTREE.csv",sep='\t', dec='.',header=T)
str(abc)
dotchart(abc$Total)

table(abc$Total)


levels(abc$IUCN) <- c("DD","CR", "DD", "EN", "EX", "NE", "LC", "LC", "NT", "NT", "NE","VU")

par(mfrow=c(2,3))
hist(log(abc$WoS+1),breaks=70,ylab="Frequency",xlab="N° of papers (log)",main="Distribution of WoS papers")
boxplot(log(abc$WoS+1)~as.factor(abc$kingdom),ylab="N° of papers (log)",main="Animal bias?")
boxplot(log(abc$WoS+1)~as.factor(abc$IUCN),ylab="N° of papers (log)", main="Are IUCN species more studied?")
plot(abc$uniqueness_family,log(abc$WoS+1),ylab="N° of papers (log)",xlab="Taxonomic uniqueness (family)",main="Are unique species more studied?")
plot(log(abc$range_size_MCP+1),log(abc$WoS+1),ylab="N° of papers (log)",xlab="Range size (log)",main="Are wide ranging species more studied?")
plot(log(abc$elevation_range+1),log(abc$WoS+1),ylab="N° of papers (log)",xlab="Elevation range (log)",main="Are wide ranging species more studied?")

par(mfrow=c(2,3))
hist(log(abc$Biodiversity...Conservation+1),breaks=70,ylab="Frequency",xlab="N° of papers in conservation (log)",main="Distribution of WoS papers")
boxplot(log(abc$Biodiversity...Conservation+1)~as.factor(abc$kingdom),ylab="N° of papers in conservationers (log)",main="Animal bias?")
boxplot(log(abc$Biodiversity...Conservation+1)~as.factor(abc$IUCN),ylab="N° of papers in conservation (log)", main="Are IUCN species more studied?")
plot(abc$uniqueness_family,log(abc$Biodiversity...Conservation+1),ylab="N° of papers in conservation (log)",xlab="Taxonomic uniqueness (family)",main="Are unique species more studied?")
plot(log(abc$range_size_MCP+1),log(abc$Biodiversity...Conservation+1),ylab="N° of papers in conservation (log)",xlab="Range size (log)",main="Are wide ranging species more studied?")
plot(log(abc$elevation_range+1),log(abc$Biodiversity...Conservation+1),ylab="N° of papers in conservation (log)",xlab="Elevation range (log)",main="Are wide ranging species more studied?")


db <- data.frame(Var1=c(1,5,6,7,8,10,2,4,6,5),Var2=c(rep(1,10)))

library("dplyr")
sampled_db <- dplyr::sample_n(db,1)

while( sum(sampled_db$Var1) < 32){ 
  
  sampled_db <- rbind(sampled_db, dplyr::sample_n(db,1) )

  }


