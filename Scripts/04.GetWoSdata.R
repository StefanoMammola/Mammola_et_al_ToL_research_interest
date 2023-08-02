## ------------------------------------------------------------------------
## 'Drivers of species knowledge across the Tree of Life' 
## ------------------------------------------------------------------------

##################
# Created by Ricardo Correia
# Updated 19/11/2022
##################

# Code to extract data from the Web of Science

# Loading R packages ------------------------------------------------------

library("wosr")
library("readxl")
library("tm")
library("tidyverse")
library("lubridate")

# Read data ---------------------------------------------------------------

sp_list <- as.data.frame(read_xlsx("./Data/SampleTREE.xlsx"))
head(sp_list)

# Open research areas
ra_list <- as.data.frame(read_xlsx("./Data/sel_reareas.xlsx"))
head(ra_list)

# Selected research areas
sel <- readRDS("./Data/sel_rareas.rds")
sel <- ra_list[sel,]

# Data frame to save selected search results
wos_results_sel <- matrix(nrow = nrow(sp_list),
                      ncol = nrow(sel)+2)
wos_results_sel <- as.data.frame(wos_results_sel)
colnames(wos_results_sel)[1:2] <- c("Species", "Total")
colnames(wos_results_sel)[3:ncol(wos_results_sel)] <- sel[,1]
wos_results_sel[,1] <- sp_list[,6]
head(wos_results_sel)

# Sampling function
extract.WoS <- function (q, s) {
  
  tryCatch(
    
    expr = {
      
      # extracting Wos records
      q_sp <- query_wos(query = q,
                        sid = s,
                        editions = c("SCI", "SSCI", "AHCI", "ISTP", "ISSHP", "ESCI", "BSCI", "BHCI")) 
      
    },
    
    #if search fails, recursively call the same function...
    error = function(e){
      
      print(" ... Internet search error occurred ... ")
      Sys.sleep(30) #pause for 30 seconds
      sid <- auth(username = NULL, password = NULL)
      print(" ... Trying again ... ")
      
      #extracting WoS records
      q_sp <- extract.WoS(q = q, s = s)
      
    }
  )
  
  return(q_sp)
}

# Extract data ------------------------------------------------------------

sid <- auth(username = NULL, password = NULL)
# Sample WoS
for (i in 1:nrow(sp_list)){
  
  if(i %in% seq(1,3200, 100)){
    sid <- auth(username = NULL, password = NULL)
  }
  
  # Set general species query
  sp_query <- paste0("(TS=(\"", sp_list[i,6], "\"))")
  q_sp <- extract.WoS(q = sp_query, s = sid)
  
  # Save results of general species query
  wos_results_sel[i,2] <- q_sp$rec_cnt
  
  if(q_sp$rec_cnt > 0){
    for (j in 1:nrow(sel)){
      sp_query_ra <- paste0("(TS=(\"", sp_list[i,6], "\") AND SU=(\"", sel[j,1],"\"))")
      q_sp_ra <- extract.WoS(q = sp_query_ra, s = sid)
      
      # Save results of general species query
      wos_results_sel[i,j+2] <- q_sp_ra$rec_cnt
      
    }  
  } else {
    wos_results_sel[i,3:ncol(wos_results_sel)] <- 0
  }
  
  print(paste0("Finished searches for species ", i, " out of ", nrow(sp_list)))
  Sys.sleep(1)
  
} #end

# Store results
SampleTree <- left_join(sp_list, wos_results_sel, by = "taxonID")

# Saving
write.csv2(SampleTre, "./Data/SampleTree.csv", row.names = F)

#end