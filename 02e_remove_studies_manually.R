# ---
# title: "MoveTraits Database"
# author: "Anne Hertel"
# date: "March 2025"
# ---

### in this script: 
## open access studies where owners have not approved reuse are set to "exclude = yes" in reference table

library(lubridate);library(metafor);library(tidyverse);library(amt);library(Hmisc)
library(adehabitatHR); library(move2); library(epitools); library(suncalc); library(purrr); library(bit64)

## ----Import movement data per individual-------------------------------------------------------------
pathTOfolder <- "./DATA/MoveTraitsData/"

metadata <- readRDS(paste0(pathTOfolder,"/referenceTableStudies_ALL_excludedColumn.rds"))
metadata[metadata$MBid %in% 1541820092,"excluded"] <- "yes" # exclude Cagan study

metadata[metadata$fileName %in%  c("1120749252_3069649656.rds",
                                   "1120749252_3069649678.rds",
                                   "1120749252_3069649685.rds",
                                   "2950149_2950167.rds") ,"excluded"] <- "yes" # exclude coaties

## ---- Remove studies from trabnslocated, semi-domesticated or experimental individuals -------------------------------------------------------------

agreements <- readRDS("./DATA/MoveTraitsData/data_agreements/data_agreements.rds")

rem <- agreements |> 
  filter(peculiarity_experimental == "Yes" | peculiarity_semi_domesticated == "Yes") |> 
  dplyr::select(Study_Name,Study_id, peculiarity_experimental,peculiarity_semi_domesticated)
#View(rem)
metadata[metadata$MBid %in% rem$Study_id ,"excluded"] <- "yes"

# check study names for peculiarities eg "animal hit by car"
#View(agreements[agreements$Study_id %nin% rem$Study_id,])
metadata[metadata$MBid %in% c(2991437203,8086049754,7865931038,8086475605,8086521811,5839913205) ,"excluded"] <- "yes"

## ---- Studies to remove - manual cleaning -------------------------------------------------------------
## ---- REVISE AFTER WE HAVE IMPROVED CLEANING ROUTINE -------------------------------------------------------------

myremove <- read.csv(paste0(pathTOfolder,"/manual_remove.csv"))
myremove <- paste(myremove$study_id,myremove$animal_id,sep="_")
myremove <- paste0(myremove,".rds")

metadata[metadata$fileName %in%  myremove ,"excluded"] <- "yes" 

saveRDS(metadata, paste0(pathTOfolder,"/referenceTableStudies_ALL_excludedColumn_excludedStudies.rds"))

