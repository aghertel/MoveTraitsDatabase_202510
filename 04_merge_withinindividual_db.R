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


pathTOfolder <- "./DATA/MoveTraitsData/Movebank/"
pathTOfolder2 <- "./DATA/MoveTraitsData/"

#-------------------------------------------------------------------------------------------
## ----Save within-individual level Database-------------------------------------------------------------
#-------------------------------------------------------------------------------------------

#dir for individual underlying traits
pthtrait <- paste0(pathTOfolder,"7.MB_indv_trait/")

flsTS <- list.files(pthtrait, full.names = T)

# Read and combine all files while keeping all columns
db.movebank <- flsTS %>%
  lapply(readRDS) %>%        
  bind_rows() 

# # sort columns
# db.movebank <- db.movebank[,c(1:82,85,89,86,92,93,83,84,90,87,91,94,95,88)]

## ----Merge individual level information data-------------------------------------------------------------

## Merge meta data and exclude duplicates
metadata <- readRDS(paste0(pathTOfolder2,"/referenceTableStudies_ALL_excludedColumn_excludedStudies.rds"))

metadata <- 
  metadata |>
  mutate(individual_id = sapply(str_split(fileName, "_|\\."), function(x) x[2])) |>
  dplyr::select(MBid,individual_id,individual_local_identifier,species,sex,animal_mass,
                animal_life_stage,manipulation_type,median_timelag_mins,
                tracking_duration_days, tracking_start_date, tracking_end_date,excluded) |> 
  rename(study_id = MBid) |> 
  mutate(study_individual = paste(study_id,individual_id,sep="_"))

metadata <- metadata |> 
  filter(!duplicated(study_individual)) |> 
  dplyr::select(study_individual,individual_local_identifier,species,sex,animal_mass,
                animal_life_stage,median_timelag_mins,
                tracking_duration_days, tracking_start_date, tracking_end_date,excluded) 

db.movebank.2 <- db.movebank |>
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

commonname <- read.csv(
  paste0(pathTOfolder2,"/SpeciesList_commonname.csv"))

commonname <-
  commonname |>
  dplyr::select("Species","common_name","class","movement.mode","include") |> 
  rename("species" = "Species") 

db.movebank.4 <- db.movebank.3 |> 
  left_join(commonname, by = "species") 

# species that were removed from database
table(db.movebank.4[db.movebank.4$include == "no",c("common_name")])

db.movebank.4 <- db.movebank.4 |> 
  filter(include == "yes") |> 
  dplyr::select(-include,-study_individual) |> 
  mutate(source = "movebank.august2026")

#nrow(db.movebank.4) == nrow(db.movebank.3)  

# 14372 individuals
dim(db.movebank.4)
# 417 studies
length(unique(db.movebank.4$study_id))

## ----Merge Tucker data-------------------------------------------------------------

db.tucker <- readRDS("DATA/MoveTraitsData/Tucker/tucker_withinindividual_20260521.rds")

db.tucker <-
  db.tucker |> 
  filter(individual_id != 8) |> 
  mutate(animal_mass = animal_mass*1000) |> 
  mutate(source = "Tucker2023")

# remove duplicate wildebeest and elk
db.tucker <- db.tucker |> 
  filter(species != "Connochaetes taurinus") |> 
  filter(!(species == "Cervus canadensis" & contact_person_name == "Mark Hebblewhite"))  

# remove czech red deer
db.tucker <- db.tucker |> 
  filter(contact_person_name != "Miloš Ježek") 
  
names(db.tucker) <- tolower(names(db.tucker)) 
names(db.movebank.4) <- tolower(names(db.movebank.4)) 

# add taxon - tucker
db.tucker <- db.tucker |> 
  left_join(commonname, by = "species") |> 
  dplyr::select(-include)

db.tucker |> group_by(class) |>  tally()

# dimensions of movebank and Tucker files
dim(db.tucker)
dim(db.movebank.4)

colnames(db.movebank.4)[colnames(db.movebank.4) %nin% colnames(db.tucker)]
colnames(db.tucker)[colnames(db.tucker) %nin% colnames(db.movebank.4)]

## bind database
withinindividual.db <- plyr::rbind.fill(db.movebank.4,db.tucker)
dim(withinindividual.db)

withinindividual.db <- withinindividual.db |> 
  dplyr::select("study_id","individual_id","individual_local_identifier",
                "species","common_name","class","movement.mode",
                "sex","animal_mass","animal_life_stage","source",
                "d1h","d24h",
                "dmax24h","dmax1m","dmax12m",
                "mcp24h","mcp1m","mcp12m",
                "iou24h","iou1m","iou12m","di",
                "median_timelag_mins","tracking_duration_days",
                "tracking_start_date","tracking_end_date","contact_person_name",
                "license_type","citation")

## ----Save within-individual level Database-------------------------------------------------------------

# 8220 bird ind., 7986 mammal ind. - 16206 ind total
withinindividual.db |> tally()
withinindividual.db |> group_by(class) |>  tally()
withinindividual.db |> group_by(source) |>  tally()

pthdb <- paste0(pathTOfolder,"8.MoveTraits_db/")
saveRDS(withinindividual.db, file="./DATA/MoveTraitsData/8.MoveTraits_db/MoveTrait.v0.1_withinindividual_20260807.rds")
