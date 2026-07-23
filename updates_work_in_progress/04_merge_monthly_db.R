# ---
# title: "MoveTraits Database"
# author: "Anne Hertel"
# date: "March 2025"
# ---

### in this script: 

# merge the databases created from movebank data and tucker data
# summarize the data at the species level
# save the database with repeated within-individual trait measures

## remove duplicates and outliers and studies without permission
## order columns in a logical way and rename if necessary


library(lubridate);library(metafor);library(tidyverse);library(amt);library(Hmisc)
library(adehabitatHR); library(move2); library(epitools); library(suncalc); library(purrr); library(bit64)

#-------------------------------------------------------------------------------------------
# ## ----Monthly level Database-------------------------------------------------------------
#-------------------------------------------------------------------------------------------

## ----Import movement data per individual-------------------------------------------------------------
pathTOfolder <- "./updates_work_in_progress/DATA/Movebank/"

#dir for individual monthly summaries
pthtraitsum <- paste0(pathTOfolder,"6.MB_indv_monthly_traitsum/")

flsTS <- list.files(pthtraitsum, full.names = T)

# Read and combine all files while keeping all columns
db.movebank <- flsTS %>%
  lapply(readRDS) %>%        # Read each file into a list of data frames
  bind_rows()                 # Combine them into one data frame

# tr <- readRDS(flsTS[8])
# colnames(tr)

# View the combined data
print(db.movebank)

dim(db.movebank)
colnames(db.movebank)

# Remove rows with all NA values
db.movebank.1 <- 
  db.movebank %>%
  filter(!if_all(c("n1h":"di.05"), is.na))

## ----Merge individual level information data-------------------------------------------------------------

## Merge meta data and exclude duplicates and studies without permission
pathTOfolder2 <- "/Users/ahertel/Documents/Work/Study_MoveTraits/database v 0.0/MoveTraitsDatabase_Git/MoveTraitsDatabase_202510/DATA/MoveTraitsData"
metadata <- readRDS(paste0(pathTOfolder2,"/referenceTableStudies_ALL_excludedColumn_excludedStudies.rds"))

metadata <- 
  metadata |>
  mutate(individual_id = sapply(str_split(fileName, "_|\\."), function(x) x[2])) |>
  dplyr::select(MBid,individual_id,species,sex,animal_mass,
                animal_life_stage,manipulation_type,median_timelag_mins,
                tracking_duration_days, tracking_start_date, tracking_end_date,excluded) |> 
  rename(study_id = MBid) |> 
  mutate(study_individual = paste(study_id,individual_id,sep="_"))

metadata <- metadata |> 
  filter(!duplicated(study_individual)) |> 
  dplyr::select(study_individual,species,sex,animal_mass,
                animal_life_stage,median_timelag_mins,
                tracking_duration_days, tracking_start_date, tracking_end_date,excluded) 

#library(bit64)
#metadata$study_id <- as.integer64(metadata$study_id)

db.movebank.2 <- db.movebank.1 |>
  mutate(study_individual = paste(study_id,individual_id,sep="_")) |> 
  left_join(metadata, by = "study_individual") |> 
  filter(excluded == "no") |> 
  dplyr::select(-excluded)

## ----Merge study level information data-------------------------------------------------------------

metadata2 <- readRDS(paste0(pathTOfolder2,"/full_table_all_studies.rds"))
colnames(metadata2)[7] <- "study_id"

db.movebank.3 <- db.movebank.2 |>
  left_join(metadata2[,c("study_id","contact_person_name","license_type","citation")], by = "study_id") 

## ----Rename species name-------------------------------------------------------------

#  Martes pennanti == Pekania pennanti
db.movebank.3[db.movebank.3$species %in% c("Pekania pennanti"),"species"] <-"Martes pennanti"

## ----Merge common name and exclude reptiles, fish etc.-------------------------------------------------------------

commonname <- read.csv("/Users/ahertel/Documents/Work/Study_MoveTraits/database v 0.0/MoveTraitsDatabase_Git/MoveTraits_Git/DATA/SpeciesList_commonname.csv")
commonname <-
  commonname |> 
  rename("species" = "Species")

db.movebank.4 <- db.movebank.3 |> 
  left_join(commonname, by = "species") 

# species that were removed from database
table(db.movebank.4[db.movebank.4$include == "no",c("common_name")])

db.movebank.4 <- db.movebank.4 |> 
  filter(include == "yes") |> 
  dplyr::select(-include,-study_individual) |> 
  mutate(source = "movebank.june2026")

#nrow(db.movebank.4) == nrow(db.movebank.3)  

# 4574 individuals
dim(db.movebank.4)
# 243 studies
length(unique(db.movebank.4$study_id))

## ----Merge Tucker data-------------------------------------------------------------

db.tucker <- readRDS("updates_work_in_progress/DATA/Tucker/tucker_monthly.sum_20260521.rds")

db.tucker <-
  db.tucker |> 
  filter(individual_id != 8) |> 
  filter(!duplicated(individual_id)) |> 
  mutate(animal_mass = animal_mass*1000) |>  # body mass in tucker in kg in movebank in grams
  mutate(source = "Tucker2023")

# Remove rows with all NA values
db.tucker <- 
  db.tucker %>%
  filter(!if_all(c("n1h":"di.05"), is.na))

# remove two studies that are duplicated in movebank! (see script 04b)
# This step needs to be automized in the future
db.tucker <- db.tucker |> 
  filter(species != "Connochaetes taurinus") |> 
  filter(!(species == "Cervus canadensis" & contact_person_name == "Mark Hebblewhite"))  

# remove czech red deer
db.tucker <- db.tucker |> 
  filter(contact_person_name != "Miloš Ježek") 

# all column names lower case  
#names(db.tucker) <- tolower(names(db.tucker)) 

# add taxon - tucker
db.tucker <- db.tucker |> 
  left_join(commonname, by = "species") |> 
  dplyr::select(-include)

# dimensions of movebank and Tucker files
dim(db.tucker)
dim(db.movebank.4)

colnames(db.movebank.4)[colnames(db.movebank.4) %nin% colnames(db.tucker)]
colnames(db.tucker)[colnames(db.tucker) %nin% colnames(db.movebank.4)]

## bind database
monthly.db <- plyr::rbind.fill(db.movebank.4,db.tucker)
dim(monthly.db)

## ----Save monthly database-------------------------------------------------------------
monthly.db <- monthly.db |> 
  dplyr::select("study_id","individual_id",
                "species","common_name",
                "month","year",
                "grid.id.100km","grid.id.10km",
                "n1h":"d1.05",
                "median_timelag_mins","tracking_duration_days","tracking_start_date","tracking_end_date",
                "contact_person_name","license_type","citation","source")

saveRDS(MoveTrait.v0.1, file=paste0(pthdb,"MoveTrait.v0.1_individual.sum_20260519.rds"))

