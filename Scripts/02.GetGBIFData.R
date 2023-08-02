## ------------------------------------------------------------------------
## 'Drivers of species knowledge across the Tree of Life' 
## ------------------------------------------------------------------------

##################
# Created by Ricardo Correia
# Updated 14/03/2022
##################

# Script to extract data from GBIF

# Loading R packages ------------------------------------------------------

library("countrycode")
library("dplyr")
library("httr") 
library("red") 
library("rgbif")
library("tidyverse")
library("xlsx")       

# Read data ---------------------------------------------------------------

SampleTree <- read.csv("./Data/SampleTREE.csv", sep = ',', dec = '.', header=T, as.is = F)

# Configuration to minimize problem with online searches ...
httr::set_config(config(ssl_verifypeer = 0L))

# Get country codes
ctr <- codelist$iso2c[!is.na(codelist$iso2c)]

# Extracting GBIF data ----------------------------------------------------

# Store results
SampleTree2 <- data.frame(taxonID = SampleTree$taxonID,
                          n_occurrences = NA,
                          n_geo_occurrences = NA,
                          n_sample_occurrences = NA,
                          centroid_lat = NA,
                          centroid_long = NA,
                          range_size = NA,
                          range_size_MCP = NA,
                          range_size_AOO = NA,
                          higherGeography = NA,
                          verbatimHabitat = NA)

#Run from here: 
#(Warning: it takes a few hours)
for(i in 1:nrow(SampleTree)) {
  
  message(paste("_ Selecting species ", i , " out of " , nrow(SampleTree), sep = ''))
  
  if(is.na(SampleTree2$taxonID[i])){
    
    next()
    
  } else {
    
    print(paste0("... Counting occurrences for species ", SampleTree$name[i], " ... "))
    
    SampleTree2$n_occurrences[i] <- occ_count(taxonKey = SampleTree2$taxonID[i])
    SampleTree2$n_geo_occurrences[i] <- occ_count(taxonKey = SampleTree2$taxonID[i],
                                                  georeferenced = T)
    
    print(paste0("... ", SampleTree2$n_geo_occurrences[i], " occurrences found for species ", SampleTree$name[i] ,"... "))
    
    if(SampleTree2$n_geo_occurrences[i] <= 1){
      
      next()
      
    } else if(SampleTree2$n_geo_occurrences[i] > 1 & SampleTree2$n_geo_occurrences[i] < 99000){
           
      while(T){
        occ_data <- try(occ_data(taxonKey = SampleTree2$taxonID[i],
                                 hasCoordinate = T,
                                 hasGeospatialIssue = F,
                                 limit = 99000),
                        silent = T)
                
        if(inherits(occ_data, "try-error")){
          message("restart loop")
          Sys.sleep(60)
        } else {
          break
        }
      }
      
      occ_data <- occ_data$data
      print(paste0("... Extracted ", nrow(occ_data), " occurrences for ", SampleTree$name[i], " ..."))
      
    } else {
      occ_data <- data.frame()
      for(j in 1:length(ctr)){
        while(T){
          occ_data_ctr <- try(occ_data(taxonKey = SampleTree2$taxonID[i],
                                       country = ctr[j],
                                       hasCoordinate = T,
                                       hasGeospatialIssue = F,
                                       limit = 99000),
                              silent = T)
          if(inherits(occ_data_ctr, "try-error")) {
            message("restart loop")
            Sys.sleep(60)
          } else {
            break
          }
        }
        
        print(paste0("... Extracted ", nrow(occ_data_ctr$data)," observation for ", SampleTree$name[i], " from country ", ctr[j], " ..."))
        
        if(!is.null(occ_data_ctr$data) && nrow(occ_data_ctr$data) != 0){
          occ_data <- bind_rows(occ_data, occ_data_ctr$data)
          rm(occ_data_ctr)
          Sys.sleep(30)
        }
        
      }
      
      print(paste0("... Extracted ", nrow(occ_data), " occurrences for ", SampleTree$name[i], " ..."))
      
    }
    
    print("... Extracting range ... ")
    coord <- na.omit(cbind(occ_data$decimalLongitude, occ_data$decimalLatitude))
    
    SampleTree2$n_sample_occurrences[i] <- nrow(occ_data)
    SampleTree2$centroid_lat[i] <- mean(coord[,2])
    SampleTree2$centroid_long[i] <- mean(coord[,1])
    SampleTree2$range_size[i] <- mean(geodist::geodist(coord[,1:2],c(SampleTree2$centroid_long[i], SampleTree2$centroid_lat[i]), measure = "geodesic"))
    coord <- data.frame(coord)
    SampleTree2$range_size_MCP[i] <- red::eoo(coord) #minimum convex polygon
    SampleTree2$range_size_AOO[i] <- red::aoo(coord) #area of occupancy
    
    print("... Extracting geography ... ")  
    
    biogeography <- unique(na.omit(occ_data$higherGeography))
    if(is.null(biogeography) == TRUE){
      SampleTree2$higherGeography[i] <- NA
    } else if(is.logical(biogeography) == TRUE) {
      SampleTree2$higherGeography[i] <- NA
    } else {
      SampleTree2$higherGeography[i] <- paste(biogeography,collapse=" ; ")
    }
    
    print("... Extracting habitat ... ")  
    
    habitat <- unique(na.omit(occ_data$habitat))
    if(is.null(habitat) == TRUE){
      SampleTree2$verbatimHabitat[i] <- NA
    } else if(is.logical(habitat) == TRUE) {
      SampleTree2$verbatimHabitat[i] <- NA
    } else {
      SampleTree2$verbatimHabitat[i] <- paste(habitat,collapse=" ; ")
    }
    
    rm(occ_data)  
  } 

  saveRDS(SampleTree2, "./SampleTREE2.rds")
  s3sync(data_folder, "gbif", direction = c("download","upload"), region ="")

}  #end

# Store results
SampleTree <- left_join(SampleTree, SampleTree2, by = "taxonID")

# Saving
write.csv2(SampleTre, "./Data/SampleTree.csv", row.names = F)

#end