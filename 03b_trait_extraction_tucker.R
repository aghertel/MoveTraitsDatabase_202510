# ---
# title: "MoveTraits Database"
# author: "Anne G. Hertel"
# date: "20/5/2026"
# ---

library(lubridate);library(metafor);library(tidyverse);library(amt);
library(adehabitatHR); library(move2); library(epitools); library(suncalc); library(purrr); library(bit64)
library(dggridR);library(mapview);library(ISOweek);library(sf)

## ----Import functions for traits-------------------------------------------------------------

source(here::here('trait_scripts', 'd1h.R'))
source(here::here('trait_scripts', 'd24h.R'))
source(here::here('trait_scripts', 'dmax24h.R'))
source(here::here('trait_scripts', 'dmax1m.R'))
source(here::here('trait_scripts', 'dmax12m.R'))
source(here::here('trait_scripts', 'mcp24h.R'))
source(here::here('trait_scripts', 'mcp1m.R'))
source(here::here('trait_scripts', 'mcp12m.R'))
source(here::here('trait_scripts', 'iou24h.R'))
source(here::here('trait_scripts', 'iou1m.R'))
source(here::here('trait_scripts', 'iou12m.R'))
source(here::here('trait_scripts', 'di.R'))

# amt function to make one track per individual with hourly steps
make_and_resample_track <- function(x, sampling.interval = 1, tolerance = 15){
  require(amt)
  make_track(x, 
             Longitude, 
             Latitude, 
             TimestampUTC, 
             individual_id = individual_id, 
             species = species, 
             contact_person_name = contact_person_name, 
             animal_mass = animal_mass,
             year = year, 
             month = month, 
             hour = hour, 
             crs = 4326) %>%
    track_resample(., rate = hours(sampling.interval),
                   tolerance = minutes(tolerance)) %>%
    return()
}

# function to determine the hour of the day with the most hourly relocations
# this hour will be used to resample the track to daily locations
get_preferred_hour <- function(trk) {
  tab <- tabulate(hour(trk$t_) + 1L, nbins = 24L)
  if (length(unique(tab)) == 1L) return(12L)
  which.max(tab) - 1L
}


daily_from_hourly <- function(trk, tolerance_mins = 60) {
  dat <- as.data.frame(trk)
  tz <- attr(dat$t_, "tzone")
  pref_hour <- get_preferred_hour(trk)

  # Vectorized: compute each row's target time based on its date
  dates <- as.Date(dat$t_, tz = tz)
  target_times <- as.POSIXct(
    paste(dates, sprintf("%02d:00:00", pref_hour)),
    tz = tz
  )
  diffs_mins <- abs(as.numeric(difftime(dat$t_, target_times, units = "mins")))

  # Filter to within-tolerance rows, then keep closest per day
  daily_df <- dat |>
    mutate(.date = dates, .diff = diffs_mins) |>
    filter(.diff <= tolerance_mins) |>
    slice_min(.diff, n = 1, by = .date, with_ties = FALSE) |>
    select(-.date, -.diff)

  if (nrow(daily_df) == 0) return(NULL)

  amt::make_track(
    daily_df,
    x_, y_, t_,
    individual_id = individual_id,
    species = species,
    contact_person_name = contact_person_name,
    animal_mass = animal_mass,
    year = year,
    month = month,
    hour = hour,
    crs = 4326
  )
}


# Aim

# We here provide code for a first version of the MoveTraits database. MoveTraits uses animal movement data collected from GPS sensors to summarize a suite of movement metrics on the individual level.

# The workflow is as follow:
# 1. resampling of the raw GPS relocation data to regular time intervals, here: 1hour, 24 hour
# 2. using resampled data to quantify movement metrics
#   a) using these resampled data to build regular movement trajectories to quantify step lengths of successive locations (at the hourly or 24 hourly rate)
#   b) using resampled data to quantify maximum displacement within a set time interval from all pairwise distance comparisons
#   c)
# 3. summarize each metric to obtain one value per individual (mean, median, cv, 5 & 95 %ile)
# 4. create a database with metrics summarized at the indivdual level AND provide the raw metrics without spatial information
# 
# We here compile a first version from open access data obtained from Tucker et al. 2023 "Behavioral responses of terrestrial mammals to COVID-19 lockdowns" (https://zenodo.org/records/7704108, file "").
# 
# We only used data from 2019 (i.e., not from 2020 during COVID lockdowns).

## Load a prepare raw spatial data


## ----load data--------------------------------------------------------------------------------------------------------
movedata <- readRDS("./DATA/Tucker/Tucker_Road_Spatial.rds")

#spatial grid
dggs.100    <- dgconstruct(projection = "ISEA", area = 10000, resround='nearest')
dggs.10     <- dgconstruct(projection = "ISEA", area = 100, resround='nearest')
dggs.1      <- dgconstruct(projection = "ISEA", area = 1, resround='nearest')

## ----select columns---------------------------------------------------------------------------------------------------
movedata <- movedata[,c("Species","ID","TimestampUTC","Latitude", "Longitude","BodyMass_kg","ContactPerson")]
colnames(movedata) <- c("species","individual_id","TimestampUTC","Latitude","Longitude", "animal_mass","contact_person_name")

#' Add time information to aggregate later; set timezone to UTC
movedata <-
  movedata %>% 
  mutate(
    year = lubridate::year(TimestampUTC),
    month = lubridate::month(TimestampUTC),
    hour = lubridate::hour(TimestampUTC))

#' Keep only data from 2019, i.e. not during Covid lockdowns
movedata <- movedata %>% 
  filter(year < 2020) 

## ----Resample data to 1hr-------------------------------------------------------------
#Resample data to 1h time scales using amt

animlocs.1hourly <- movedata %>%
  group_by(individual_id) %>%
  group_split() %>%
  map(~ make_and_resample_track(., sampling.interval = 1,
                                 tolerance = 15)) %>%
  setNames(sort(unique(movedata$individual_id)))

## ----Resample data to 24 hrs-------------------------------------------------------------
#Resample data to 24h time scales using amt

animlocs.daily <- animlocs.1hourly  |>
  map(~ daily_from_hourly(.x, tolerance_mins = 60))

## ----Calculate movement metrics-------------------------------------------------------------

#1h displacement----
d1h <- purrr::map_dfr(animlocs.1hourly,calc_d1h,dggs.10,dggs.1,.id = "track_name")
sum.ind.d1h <- f_sum.ind.d1h(d1h) 
#sum.monthly.ind.d1h <- f_sum.monthly.ind.d1h(d1h)
d1h_split <- split(d1h, d1h$individual_id)
sum.monthly.ind.d1h <- lapply(d1h_split, f_sum.monthly.ind.d1h)
sum.monthly.ind.d1h <- do.call(rbind, sum.monthly.ind.d1h)

#24hr displacement distance----
d24h <- purrr::map_dfr(animlocs.daily,calc_d24h,dggs.10,dggs.1,.id = "track_name")
sum.ind.d24h <- f_sum.ind.d24h(d24h)
d24h_split <- split(d24h, d24h$individual_id)
sum.monthly.ind.d24h <- lapply(d24h_split, f_sum.monthly.ind.d24h)
sum.monthly.ind.d24h <- do.call(rbind, sum.monthly.ind.d24h)

#Maximum 24hr displacement distance----
dmax24h <- purrr::map_dfr(animlocs.1hourly,calc_dmax24h,dggs.10,dggs.1,.id = "track_name")
sum.ind.dmax24h <- f_sum.ind.dmax24h(dmax24h)
dmax24h_split <- split(dmax24h, dmax24h$individual_id)
sum.monthly.ind.dmax24h <- lapply(dmax24h_split, f_sum.monthly.ind.dmax24h)
sum.monthly.ind.dmax24h <- do.call(rbind, sum.monthly.ind.dmax24h)

#Maximum 1month displacement distance----
dmax1m <- purrr::map_dfr(animlocs.daily,calc_dmax1m,dggs.10,dggs.1,.id = "track_name")
sum.ind.dmax1m <- f_sum.ind.dmax1m(dmax1m)
sum.monthly.ind.dmax1m <- f_sum.monthly.ind.dmax1m(dmax1m)

#Daily MCP----
mcp24h <- purrr::map_dfr(animlocs.1hourly,calc_mcp24h,dggs.10,dggs.1,.id = "track_name")
sum.ind.mcp24h <- f_sum.ind.mcp24h(mcp24h)
mcp24h_split <- split(mcp24h, mcp24h$individual_id)
sum.monthly.ind.mcp24h <- lapply(mcp24h_split, f_sum.monthly.ind.mcp24h)
sum.monthly.ind.mcp24h <- do.call(rbind, sum.monthly.ind.mcp24h)

#Monthly MCP----
mcp1m <- purrr::map_dfr(animlocs.daily,calc_mcp1m,dggs.10,dggs.1,.id = "track_name")
sum.ind.mcp1m <- f_sum.ind.mcp1m(mcp1m)
sum.monthly.ind.mcp1m <- f_sum.monthly.ind.mcp1m(mcp1m)

#Daily IOU----
iou24h <- calc_iou24h(mcp24h,d1h, dggs.10, dggs.1)
sum.ind.iou24h <- f_sum.ind.iou24h(iou24h)
iou24h_split <- split(iou24h, iou24h$individual_id)
sum.monthly.ind.iou24h <- lapply(iou24h_split, f_sum.monthly.ind.iou24h)
sum.monthly.ind.iou24h <- do.call(rbind, sum.monthly.ind.iou24h)

#Monthly IOU----
iou1m <- calc_iou1m(mcp1m,d24h, dggs.10, dggs.1)
sum.ind.iou1m <- f_sum.ind.iou1m(iou1m)
sum.monthly.ind.iou1m <- f_sum.monthly.ind.iou1m(iou1m) 

#Diurnality Index----
di <- calc_di(d1h, dggs.10, dggs.1)
sum.ind.di <- f_sum.ind.di(di)
di_split <- split(di, di$individual_id)
sum.monthly.ind.di <- lapply(di_split, f_sum.monthly.ind.di)
sum.monthly.ind.di <- do.call(rbind, sum.monthly.ind.di)

## ----Build database with individual summary values-------------------------------------------------------------
individual.traits.tucker <- 
  full_join(sum.ind.d1h, full_join(sum.ind.d24h, full_join(sum.ind.dmax24h, 
             full_join(sum.ind.dmax1m, full_join(sum.ind.mcp24h, full_join(sum.ind.mcp1m, 
             full_join(sum.ind.iou24h, full_join(sum.ind.iou1m, sum.ind.di))))))))

movedata2 <- movedata %>% 
  dplyr::select(individual_id,species,animal_mass,contact_person_name) %>% 
  dplyr::filter(!duplicated(individual_id))

individual.traits.tucker <- merge(movedata2, individual.traits.tucker, by = "individual_id", all.y=T)

#gridded coordinates for individual summaries
get_grid_ids <- function(trk, dggs_100) {
  cell_info.100 <- dgGEO_to_SEQNUM(dggs_100, trk$x_, trk$y_)
  unique(cell_info.100$seqnum)
}

grid.id.100km <- purrr::map(animlocs.1hourly, get_grid_ids, dggs = dggs.100)
grid.id.10km <- purrr::map(animlocs.1hourly, get_grid_ids, dggs = dggs.10)

grid.df <- data.frame(
  individual_id = names(grid.id.100km),
  grid.id.100km = sapply(grid.id.100km, paste, collapse = ";"),
  grid.id.10km = sapply(grid.id.10km, paste, collapse = ";"),
  row.names = NULL)

individual.traits.tucker <- individual.traits.tucker  |> 
  left_join(grid.df, by = "individual_id") |> 
  droplevels()

saveRDS(individual.traits.tucker,"updates_work_in_progress/DATA/Tucker/tucker_individual.sum_20260521.rds")

## ----Build database with individual monthly summary values-------------------------------------------------------------
dfs <- list(
  sum.monthly.ind.d1h,
  sum.monthly.ind.d24h,
  sum.monthly.ind.dmax24h,
  sum.monthly.ind.dmax1m,
  sum.monthly.ind.mcp24h,
  sum.monthly.ind.mcp1m,
  sum.monthly.ind.iou24h,
  sum.monthly.ind.iou1m,
  sum.monthly.ind.di)

# check for duplicates in monthly traits
# dup_check <- purrr::imap_dfr(dfs, function(df, nm) {
#   df %>%
#     dplyr::mutate(individual_id = as.character(individual_id)) %>%
#     dplyr::count(individual_id, month, year, name = "n") %>%
#     dplyr::filter(n > 1) %>%
#     dplyr::mutate(table = nm)
# })
# dup_check

individual.monthly.traits.tucker <-
  purrr::reduce(
  dfs,
  dplyr::full_join,
  by = c("individual_id", "month", "year"))|>
  filter(if_any(everything(), ~ !is.na(.)))

individual.monthly.traits.tucker <- merge(movedata2, individual.monthly.traits.tucker, by = "individual_id", all.y=T)

#gridded coordinates for individual summaries
get_monthly_grids <- function(trk, dggs) {
  id_monthly <- trk |>
    mutate(
      month = lubridate::month(t_),
      year = lubridate::year(t_),
      month_year = paste(year, month, sep = "_")
    )
  
  cell_info <- dgGEO_to_SEQNUM(dggs, id_monthly$x_, id_monthly$y_)
  id_monthly$grid.id <- cell_info$seqnum
  
  id_monthly |>
    distinct(individual_id, month_year, grid.id, .keep_all = TRUE) |>
    group_by(individual_id, month_year, year, month) |>
    summarise(
      grid.id = paste(unique(grid.id), collapse = ";"),
      .groups = "drop"
    )
}

grid.id.100km <- purrr::map(animlocs.1hourly, get_monthly_grids, dggs = dggs.100)
grid.id.100km <- purrr::map_dfr(grid.id.100km, ~ .x, .id = "individual_id")
colnames(grid.id.100km)[5] <- "grid.id.100km"

grid.id.10km <- purrr::map(animlocs.1hourly, get_monthly_grids, dggs = dggs.10)
grid.id.10km <- purrr::map_dfr(grid.id.10km, ~ .x, .id = "individual_id")
colnames(grid.id.10km)[5] <- "grid.id.10km"

individual.monthly.traits.tucker <- individual.monthly.traits.tucker  |> 
  left_join(grid.df, by = "individual_id") |> 
  droplevels()

saveRDS(individual.monthly.traits.tucker,"updates_work_in_progress/DATA/Tucker/tucker_monthly.sum_20260521.rds")

## ----Build database with repeated trait measures-------------------------------------------------------------
repeated.traits.tucker <- 
  movedata2  |> 
  tidyr::nest(data = -individual_id)|> 
  dplyr::select(-data) |>  # just a quick trick to keep things flowing.
  left_join(movedata2[!duplicated(movedata2$individual_id),], 
            by = c("individual_id" = "individual_id"))  |>  
  # 1 hourly
  left_join(d1h %>%
              data.frame %>% 
              tidyr::nest(data = -individual_id), 
            by = c("individual_id" = "individual_id")) %>% 
  dplyr::rename(d1h = data) |>  
  
  # 24 hourly
  left_join(d24h %>%
              tibble() |> 
              tidyr::nest(data = -individual_id), 
            by = c("individual_id" = "individual_id")) %>% 
  dplyr::rename(d24h = data) |>  
  
  # Dmax24
  left_join(dmax24h %>%
              tidyr::nest(data = -individual_id), 
            by = c("individual_id" = "individual_id")) %>% 
  dplyr::rename(dmax24h = data)  |>  
  
  # Dmax1m
  left_join(dmax1m %>%
              tidyr::nest(data = -individual_id), 
            by = c("individual_id" = "individual_id")) %>% 
  dplyr::rename(dmax1m = data)  |>  
  
  # mcp.daily
  left_join(mcp24h %>%
              tidyr::nest(data = -individual_id), 
            by = c("individual_id" = "individual_id")) %>% 
  dplyr::rename(mcp24h = data) |>  
  
  # mcp.monthly
  left_join(mcp1m %>%
              tidyr::nest(data = -individual_id), 
            by = c("individual_id" = "individual_id")) %>% 
  dplyr::rename(mcp1m = data)  |>  
  
  # Intensity of Use 24h
  left_join(iou24h %>%
              tibble() |> 
              tidyr::nest(data = -individual_id), 
            by = c("individual_id" = "individual_id")) %>% 
  dplyr::rename(iou24h = data)  |>  
  
  # Intensity of Use 1m
  left_join(iou1m %>%
              tibble() |> 
              tidyr::nest(data = -individual_id), 
            by = c("individual_id" = "individual_id")) %>% 
  dplyr::rename(iou1m = data)  |>  
  
  # Diurnality Index
  left_join(di %>%
              tidyr::nest(data = -individual_id), 
            by = c("individual_id" = "individual_id")) %>% 
  dplyr::rename(di = data) 

## ---------------------------------------------------------------------------------------------------------------------
saveRDS(repeated.traits.tucker,"updates_work_in_progress/DATA/Tucker/tucker_withinindividual_20260521.rds")  
