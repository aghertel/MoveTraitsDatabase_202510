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
# ## ----Individual level Database-------------------------------------------------------------
#-------------------------------------------------------------------------------------------

## ----Import movement data per individual-------------------------------------------------------------
pathTOfolder <- "./updates_work_in_progress/DATA/"

#dir for individual summaries
pthtraitsum <- paste0(pathTOfolder,"5.MB_indv_traitsum/")

flsTS <- list.files(pthtraitsum, full.names = T)

# Read and combine all files while keeping all columns
db.movebank <- flsTS %>%
  lapply(readRDS) %>%        # Read each file into a list of data frames
  bind_rows()                 # Combine them into one data frame

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

db.tucker <- readRDS("./DATA/Tucker/MoveTraitsDB.v0.1_Tucker.rds")

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
names(db.tucker) <- tolower(names(db.tucker)) 

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
MoveTrait.v0.1 <- plyr::rbind.fill(db.movebank.4,db.tucker)
dim(MoveTrait.v0.1)

MoveTrait.v0.1 <- MoveTrait.v0.1 |> 
  dplyr::select("study_id","individual_id",
                "species","common_name","class","movement.mode",
                "sex","animal_mass","animal_life_stage","source",
                "grid.id.100km":"di.05","median_timelag_mins","tracking_duration_days",
                "tracking_start_date","tracking_end_date","contact_person_name",
                "license_type","citation")

## ----Save individual level Database-------------------------------------------------------------

# final recode of species labels 
MoveTrait.v0.1 <- MoveTrait.v0.1 |> 
  mutate(species = fct_recode(species, "Ovis canadensis" = "Ovis canadensis californiana")) |> 
  mutate(species = fct_recode(species, "Ovis canadensis" = "Ovis canadensis canadensis")) |> 
  mutate(species = fct_recode(species, "Ovis canadensis" = "Ovis canadensis nelsoni"))|> 
  mutate(species = fct_recode(species, "Loxodonta africana" = "Elephantidae"))|> 
  mutate(species = fct_recode(species, "Cervus elaphus" = "Cervus canadensis"))

# 108 bird sp., 55 mammal sp.
MoveTrait.v0.1 |> filter(!duplicated(species)) |> group_by(class) |>  tally()
# 3660 bird ind., 2691 mammal ind. - 6351 ind total
MoveTrait.v0.1 |> tally()
MoveTrait.v0.1 |> group_by(class) |>  tally()
#1777 tucker, 4574 movebank
MoveTrait.v0.1 |> group_by(source) |>  tally()

dir.create(paste0(pathTOfolder,"8.MoveTraits_db"))
pthdb <- paste0(pathTOfolder,"8.MoveTraits_db/")
saveRDS(MoveTrait.v0.1, file=paste0(pthdb,"MoveTrait.v0.1_individual.sum_20260519.rds"))


#-------------------------------------------------------------------------------------------
# ## ----Species level Database-------------------------------------------------------------
#-------------------------------------------------------------------------------------------

cv_fun <- function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)

MoveTrait.v0.1.sp <- MoveTrait.v0.1 |>
  mutate(common_name = recode(common_name, "reindeer" = "reindeer/caribou")) |>
  mutate(common_name = recode(common_name, "elk" = "red deer/elk")) |>
  mutate(common_name = recode(common_name, "red deer" = "red deer/elk"))

MoveTrait.v0.1.sp2 <-
  MoveTrait.v0.1.sp |>
  group_by(species) |>
  mutate(species = unique(species),
          common_name = unique(common_name),
          class = unique(class),
          movement.mode = unique(movement.mode),
          grid.id.100km = paste(unique(trimws(unlist(strsplit(na.omit(grid.id.100km), ";")))),collapse = ";"),
          grid.id.10km = paste(unique(trimws(unlist(strsplit(na.omit(grid.id.10km), ";")))),collapse = ";"),
         across(
           ends_with(".mean"),
           list(
             species_mean = ~ mean(.x, na.rm = TRUE),
             species_sd   = ~ sd(.x, na.rm = TRUE),
             species_cv   = ~ cv_fun(.x),
             species_n    = ~ sum(!is.na(.x))),
           .names = "{gsub('.mean', '', .col, fixed = TRUE)}_{.fn}"),
            contact_person_name = paste(unique(contact_person_name), collapse = ", ")) |>
  ungroup() |>
  dplyr::select(species,common_name,class,movement.mode,grid.id.100km,grid.id.10km,
                contact_person_name,ends_with("_mean"),ends_with("_mean"),
                ends_with("_cv"),ends_with("_n")) |> 
  distinct(species, .keep_all = TRUE) 


saveRDS(MoveTrait.v0.1.sp2, file=paste0(pthdb,"MoveTrait.v0.1_species.sum_20260519.rds"))

