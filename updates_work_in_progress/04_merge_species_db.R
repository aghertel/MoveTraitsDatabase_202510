# ---
# title: "MoveTraits Database"
# author: "Anne Hertel"
# date: "May 2026"
# ---

### in this script: 

# merge the databases created from movebank data and tucker data
# summarize the individual data at the species level
# save the database with repeated within-individual trait measures

## remove duplicates and outliers and studies without permission
## order columns in a logical way and rename if necessary
library(lubridate);library(metafor);library(tidyverse);library(amt);library(Hmisc)
library(adehabitatHR); library(move2); library(epitools); library(suncalc); library(purrr); library(bit64)

pathTOfolder <- "./updates_work_in_progress/DATA/"


#-------------------------------------------------------------------------------------------
# ## ----Species level Database-------------------------------------------------------------
#-------------------------------------------------------------------------------------------

MoveTrait.v0.1 <- readRDS(file=paste0(pthdb,"MoveTrait.v0.1_individual.sum_20260519.rds"))

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
