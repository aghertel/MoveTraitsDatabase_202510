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

#-------------------------------------------------------------------------------------------
# ## ----Species level Database-------------------------------------------------------------
#-------------------------------------------------------------------------------------------

MoveTrait.v0.1 <- readRDS(file="./DATA/MoveTraitsData/8.MoveTraits_db/MoveTrait.v0.1_individual.sum_20260807.rds")

cv_fun <- function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)

MoveTrait.v0.1.sp <- MoveTrait.v0.1 |>
  mutate(common_name = recode(common_name, "caribou" = "reindeer/caribou")) |>
  mutate(common_name = recode(common_name, "reindeer" = "reindeer/caribou")) |>
  mutate(common_name = recode(common_name, "elk" = "red deer/elk")) |>
  mutate(common_name = recode(common_name, "red deer" = "red deer/elk")) |> 
  mutate(common_name = recode(common_name, "african bush elephant" = "african elephant")) |> 
  mutate(common_name = recode(common_name, "sierra nevada bighorn sheep" = "bighorn sheep")) |> 
  mutate(common_name = recode(common_name, "rocky mountain bighorn sheep" = "bighorn sheep"))|> 
  mutate(common_name = recode(common_name, "desert bighorn sheep" = "bighorn sheep"))

MoveTrait.v0.1.sp2 <-
  MoveTrait.v0.1.sp |>
  group_by(species) |>
  mutate(species = unique(species),
         common_name = unique(common_name),
         class = unique(class),
         movement.mode = unique(movement.mode),
         n.ind = n(),
         grid.id.100km = paste(unique(trimws(unlist(strsplit(na.omit(grid.id.100km), ";")))),collapse = ";"),
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
  dplyr::select(species,common_name,n.ind,class,movement.mode,grid.id.100km,grid.id.10km,
                contact_person_name,starts_with("d1h_"),starts_with("d24h_"),
                starts_with("dmax24h_"),starts_with("dmax1m_"),starts_with("dmax12m_"),
                starts_with("mcp24h_"),starts_with("mcp1m_"),starts_with("mcp12m_"),
                starts_with("iou24h_"),starts_with("iou1m_"),starts_with("iou12m_"),
                starts_with("di_")) |> 
  distinct(species, .keep_all = TRUE) 

saveRDS(MoveTrait.v0.1.sp2, file="./DATA/MoveTraitsData/8.MoveTraits_db/MoveTrait.v0.1_species.sum_20260807.rds")

