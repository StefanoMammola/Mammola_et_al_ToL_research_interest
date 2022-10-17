## ------------------------------------------------------------------------
## 'Species popularity and research interests across the Tree of Life' 
## ------------------------------------------------------------------------

##################
# Created by Stefano Mammola
# Updated 17/03/2022
##################

# Script to extract data from IUCN

# Loading R packages ------------------------------------------------------

library("httr")       #<---- to set internet settings
library("rredlist")   #<---- to interact with IUCN list
library("sp")         #<---- miscellaneous functions
library("xlsx")       #<---- to save .xlsx files

# Function to extract IUCN data avoiding errors in internet searches:
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

# Read data ---------------------------------------------------------------

SampleTree <- read.csv("./Data/SampleTREE.csv", sep = ',', dec = '.', header=T, as.is = F)

# Extracting IUCN data ----------------------------------------------------

# Defining an API for IUCN
api_IUCN <- " " #------> Add your IUCN API key 

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

# Saving
write.csv2(SampleTree, "./Data/SampleTree.csv", row.names = F)

#end
