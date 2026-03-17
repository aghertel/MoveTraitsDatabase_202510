# ---
# title: "MoveTraits Database"
# author: "Anne G. Hertel"
# date: "17/3/2026"
# ---

library(lubridate);library(metafor);library(tidyverse);library(amt);
library(adehabitatHR); library(move2); library(epitools); library(suncalc); library(purrr); library(bit64)
library(dggridR);library(mapview);library(ISOweek)

## ----Import functions for traits-------------------------------------------------------------

source(here::here('updates_work_in_progress/trait_scripts', 'd1h.R'))
source(here::here('updates_work_in_progress/trait_scripts', 'd24h.R'))
source(here::here('updates_work_in_progress/trait_scripts', 'dmax24h.R'))
source(here::here('updates_work_in_progress/trait_scripts', 'dmax7d.R'))
source(here::here('updates_work_in_progress/trait_scripts', 'dmax12m.R'))
source(here::here('updates_work_in_progress/trait_scripts', 'mcp24h.R'))
source(here::here('updates_work_in_progress/trait_scripts', 'mcp7d.R'))
source(here::here('updates_work_in_progress/trait_scripts', 'mcp1m.R'))
source(here::here('updates_work_in_progress/trait_scripts', 'mcp12m.R'))
source(here::here('updates_work_in_progress/trait_scripts', 'iou24h.R'))
source(here::here('updates_work_in_progress/trait_scripts', 'iou1m.R'))
source(here::here('updates_work_in_progress/trait_scripts', 'iou12m.R'))
source(here::here('updates_work_in_progress/trait_scripts', 'di.R'))

## ----Import movement data per individual-------------------------------------------------------------
pathTOfolder <- "./DATA/MoveTraitsData/"
pthamt1h <- paste0(pathTOfolder,"4.MB_indv_amt_1h/")

#dir for individual summaries
pathTOfolder2 <- "./updates_work_in_progress/DATA/"
dir.create(paste0(pathTOfolder2,"5.MB_indv_traitsum"))
pthtraitsum <- paste0(pathTOfolder2,"5.MB_indv_traitsum/")

#dir for individual monthly summaries
pathTOfolder2 <- "./updates_work_in_progress/DATA/"
dir.create(paste0(pathTOfolder2,"6.MB_indv_monthly_traitsum"))
pthtraitsummonthly <- paste0(pathTOfolder2,"6.MB_indv_monthly_traitsum/")

#dir for individual underlying traits
dir.create(paste0(pathTOfolder2,"7.MB_indv_trait"))
pthtrait <- paste0(pathTOfolder2,"7.MB_indv_trait/")

#spatial grid
dggs.100    <- dgconstruct(projection = "ISEA", area = 10000, resround='nearest')
dggs.10     <- dgconstruct(projection = "ISEA", area = 100, resround='nearest')
dggs.1      <- dgconstruct(projection = "ISEA", area = 1, resround='nearest')

# list files
flsMV <- list.files(pthamt1h, full.names = F)
done <- list.files(pthtraitsum, full.names = F)
flsMV <- flsMV[!flsMV%in%done]

#filter excluded studies
referenceTableStudies <- readRDS(paste0(pathTOfolder,"referenceTableStudies_ALL_excludedColumn.rds"))
referenceTableStudiesUsed <- referenceTableStudies[referenceTableStudies$excluded=="no",]

flsMV <- flsMV[flsMV %in% referenceTableStudiesUsed$fileName]

# load data of one individual and apply all trait extraction functions its tracking data 
#lapply(flsMV, function(indPth)
#  {

# example ciconia ciconia - 10596067_10666985.rds
# example anas platyrhynchos - #446579_7945602.rds
# example panthera leo - 3809257699_3809647332.rds

  animlocs.1hourly <- readRDS(paste0(pthamt1h,"3809257699_3809647332.rds")) 
  
## ----Resample data-------------------------------------------------------------
#Resample data to 24h, 7 week time scales using amt

animlocs.daily <- animlocs.1hourly |>
  track_resample(rate = hours(24),
               tolerance = minutes(60))

animlocs.weekly <- animlocs.1hourly |>
  track_resample(rate = hours(24*7),
                 tolerance = minutes(60*24))

## ----Trait extraction-------------------------------------------------------------

#1h displacement----
d1h <- calc_d1h(animlocs.1hourly, dggs.10, dggs.1)
sum.ind.d1h <- f_sum.ind.d1h(d1h)
sum.monthly.ind.d1h <- f_sum.monthly.ind.d1h(d1h)

#24hr displacement distance----
d24h <- calc_d24h(animlocs.daily, dggs.10, dggs.1)
sum.ind.d24h <- f_sum.ind.d24h(d24h)
sum.monthly.ind.d24h <- f_sum.monthly.ind.d24h(d24h)

#Maximum 24hr displacement----
dmax24h <- calc_dmax24h(animlocs.1hourly, dggs.10, dggs.1)
sum.ind.dmax24h <- f_sum.ind.dmax24h(dmax24h)
sum.monthly.ind.dmax24h <- f_sum.monthly.ind.dmax24h(dmax24h)

#Maximum 7day displacement distance----
dmax7d <- calc_dmax7d(animlocs.daily, dggs.10, dggs.1)
sum.ind.dmax7d <- f_sum.ind.dmax7d(dmax7d)
sum.monthly.ind.dmax7d <- f_sum.monthly.ind.dmax7d(dmax7d)

#Maximum annual displacement distance----
dmax12m <- calc_dmax12m(animlocs.weekly, dggs.10, dggs.1)
sum.ind.dmax12m <- f_sum.ind.dmax12m(dmax12m)

#Daily MCP----
mcp24h <- calc_mcp24h(animlocs.1hourly, dggs.10, dggs.1)
sum.ind.mcp24h <- f_sum.ind.mcp24h(mcp24h)
sum.monthly.ind.mcp24h <- f_sum.monthly.ind.mcp24h(mcp24h)

#Weekly MCP----
mcp7d <- calc_mcp7d(animlocs.1hourly, dggs.10, dggs.1)
sum.ind.mcp7d <- f_sum.ind.mcp7d(mcp7d)
sum.monthly.ind.mcp7d <- f_sum.monthly.ind.mcp7d(mcp7d)

#Monthly MCP----
mcp1m <- calc_mcp1m(animlocs.daily, dggs.10, dggs.1)
sum.ind.mcp1m <- f_sum.ind.mcp1m(mcp1m)
sum.monthly.ind.mcp1m <- mcp1m[,1:4] |> 
  rename(mcp1m = area) |> 
  mutate(year = as.numeric(substr(year_month, 1, 4))) |> 
  dplyr::select(individual_id,month,year,mcp1m)

#Annual MCP----
mcp12m <- calc_mcp12m(animlocs.weekly, dggs.10, dggs.1)
sum.ind.mcp12m <- f_sum.ind.mcp12m(mcp12m)

#Daily IOU----
iou24h <- calc_iou24h(mcp24h,d1h, dggs.10, dggs.1)
sum.ind.iou24h <- f_sum.ind.iou24h(iou24h)
sum.monthly.ind.iou24h <- f_sum.monthly.ind.iou24h(iou24h)

#Monthly IOU----
iou1m <- calc_iou1m(mcp1m,d24h, dggs.10, dggs.1)
sum.ind.iou1m <- f_sum.ind.iou1m(iou1m)
sum.monthly.ind.iou1m <- iou1m[,c(1,2,3,5)] 

#Annual IOU----
iou12m <- calc_iou12m(mcp12m,d24h, dggs.10, dggs.1)
sum.ind.iou12m <- f_sum.ind.iou12m(iou12m)

#Diurnality Index----
di <- calc_di(d1h, dggs.10, dggs.1)
sum.ind.di <- f_sum.ind.di(di)
sum.monthly.ind.di <- f_sum.monthly.ind.di(di)

## ----Individual summary database--------------------
#cbind individual summaries of all traits 
MoveTrait.ind.sum <- bind_cols(sum.ind.d1h, 
                          bind_cols(sum.ind.d24h[,2:length(sum.ind.d24h)], 
                          bind_cols(sum.ind.dmax24h[,2:length(sum.ind.dmax24h)], 
                          bind_cols(sum.ind.dmax7d[,2:length(sum.ind.dmax7d)], 
                          bind_cols(sum.ind.dmax12m[,2:length(sum.ind.dmax12m)], 
                          bind_cols(sum.ind.mcp24h[,2:length(sum.ind.mcp24h)], 
                          bind_cols(sum.ind.mcp7d[,2:length(sum.ind.mcp7d)], 
                          bind_cols(sum.ind.mcp1m[,2:length(sum.ind.mcp1m)], 
                          bind_cols(sum.ind.mcp12m[,2:length(sum.ind.mcp12m)], 
                          bind_cols(sum.ind.iou24h[,2:length(sum.ind.iou24h)],
                          bind_cols(sum.ind.iou1m[,2:length(sum.ind.iou1m)],
                          bind_cols(sum.ind.iou12m[,2:length(sum.ind.iou12m)], 
                                    sum.ind.di[,2:length(sum.ind.di)])))))))))))) %>% 
  filter(if_any(everything(), ~ !is.na(.)))

library(bit64)
movedata2 <- animlocs.1hourly %>% 
  dplyr::select(study_id,individual_id) %>% 
  mutate(individual_id = as.character(individual_id)) %>% 
  distinct()

#gridded coordinates for individual summaries
cell_info.100 <- dgGEO_to_SEQNUM(dggs.100, animlocs.1hourly$x_, animlocs.1hourly$y_)
grid.id.100km <- unique(cell_info.100$seqnum)
lon.100km <- dgSEQNUM_to_GEO(dggs.100, grid.id.100km)$lon_deg
lat.100km <- dgSEQNUM_to_GEO(dggs.100, grid.id.100km)$lat_deg

cell_info.10 <- dgGEO_to_SEQNUM(dggs.10, animlocs.1hourly$x_, animlocs.1hourly$y_)
grid.id.10km <- unique(cell_info.10$seqnum)
lon.10km <- dgSEQNUM_to_GEO(dggs.10, grid.id.10km)$lon_deg
lat.10km <- dgSEQNUM_to_GEO(dggs.10, grid.id.10km)$lat_deg

movedata2$grid.id.100km <- paste(grid.id.100km, collapse = ";")
movedata2$grid.id.10km <- paste(grid.id.10km, collapse = ";")

#join database
MoveTrait.ind.sum <- 
  if (is.null(MoveTrait.ind.sum) | nrow(MoveTrait.ind.sum) == 0) NULL else {
    movedata2 %>% 
    left_join(MoveTrait.ind.sum, by = "individual_id") %>% 
    droplevels() }

#Save database with summaries
saveRDS(MoveTrait.ind.sum, file=paste0(pthtraitsum,"3809257699_3809647332.rds"))

## ----Individual monthly summary database--------------------

MoveTrait.monthly.ind.sum <- 
sum.monthly.ind.d1h |>  
  full_join(sum.monthly.ind.d24h, 
            by = join_by(individual_id, month, year)) |> 
  full_join(sum.monthly.ind.dmax24h, 
            by = join_by(individual_id, month, year)) |> 
  full_join(sum.monthly.ind.dmax7d, 
            by = join_by(individual_id, month, year)) |> 
  full_join(sum.monthly.ind.mcp24h, 
            by = join_by(individual_id, month, year)) |> 
  full_join(sum.monthly.ind.mcp7d, 
            by = join_by(individual_id, month, year)) |> 
  full_join(sum.monthly.ind.mcp1m, 
            by = join_by(individual_id, month, year)) |> 
  full_join(sum.monthly.ind.iou24h, 
            by = join_by(individual_id, month, year)) |> 
  full_join(sum.monthly.ind.iou1m, 
            by = join_by(individual_id, month, year)) |> 
  full_join(sum.monthly.ind.di, 
            by = join_by(individual_id, month, year)) |>  
  filter(if_any(everything(), ~ !is.na(.)))
  
# #gridded coordinates for monthly individual summaries
id_monthly <- animlocs.1hourly |>
  mutate(month = month(t_),
         year = year(t_),
         month_year = paste(month,year,sep="_")) |> 
  mutate(individual_id = as.character(individual_id))

cell_info.100 <- dgGEO_to_SEQNUM(dggs.100, animlocs.1hourly$x_, animlocs.1hourly$y_)
id_monthly$grid.id.100km <- cell_info.100$seqnum
id_monthly_100km <- id_monthly |>
  distinct(month_year,grid.id.100km,.keep_all = TRUE)|> 
  group_by(individual_id,study_id,month, year) |>               
  summarise(grid.id.100km = paste(grid.id.100km, collapse = ";"),
            .groups = "drop")

cell_info.10 <- dgGEO_to_SEQNUM(dggs.10, animlocs.1hourly$x_, animlocs.1hourly$y_)
id_monthly$grid.id.10km <- cell_info.10$seqnum
id_monthly_10km <- id_monthly |>
  distinct(month_year,grid.id.10km,.keep_all = TRUE) |> 
  group_by(individual_id,study_id,month, year) |>               
  summarise(grid.id.10km = paste(grid.id.10km, collapse = ";"),
            .groups = "drop")

#join database
MoveTrait.monthly.ind.sum <- 
  if (is.null(MoveTrait.monthly.ind.sum) | nrow(MoveTrait.monthly.ind.sum) == 0) NULL else {
    id_monthly_100km  |>  
      left_join(id_monthly_10km, by =  join_by(individual_id, study_id, month, year)) %>% 
      left_join(MoveTrait.monthly.ind.sum, by =  join_by(individual_id, month, year)) %>% 
      droplevels()  }

#Save database with individual monthly summaries
saveRDS(MoveTrait.monthly.ind.sum, file=paste0(pthtraitsummonthly,"3809257699_3809647332.rds"))

## ----Build full database including raw metrics data--------------------

MoveTrait.v0.1_spatial <- 
  if (is.null(MoveTrait.v0.1)) NULL else {
  MoveTrait.v0.1 %>%
  tidyr::nest(data = -individual_id) %>% 
  dplyr::select(-data) %>% # just a quick trick to keep things flowing.
  left_join(., MoveTrait.v0.1[!duplicated(MoveTrait.v0.1$individual_id)], 
            by = c("individual_id" = "individual_id")) %>% 
  
  # 1 hourly
  {if (!is.null(animlocs.1hourly_sl)) left_join(.,  animlocs.1hourly_sl %>%
                                    data.frame %>%
                                    mutate(individual_id = as.character(individual_id)) %>%
                                    dplyr::select("individual_id","t_","d1h","x_","y_",
                                                  "grid.id.10km","lon.10km","lat.10km",
                                                  "grid.id.1km","lon.1km","lat.1km" ) %>% 
                tidyr::nest(Displ.1h = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 
  
  # 24 hourly
  {if (!is.null(animlocs.daily_sl)) left_join(.,  animlocs.daily_sl %>%
                                    mutate(individual_id = as.character(individual_id)) %>%
                                    dplyr::select("individual_id","t_","d24h","x_","y_",
                                                  "grid.id.10km","lon.10km","lat.10km",
                                                  "grid.id.1km","lon.1km","lat.1km" ) %>% 
                tidyr::nest(Displ.24h = -individual_id), 
            by = c("individual_id" = "individual_id")) else .}  %>% 
  
  # Dmax24
  {if (!is.null(dmax24)) left_join(.,  dmax24 %>%
                                    mutate(individual_id = as.character(individual_id)) %>%
                                    dplyr::select("individual_id","ymd","dmax24h","mean.x", "mean.y",
                                                  "grid.id.10km","lon.10km","lat.10km",
                                                  "grid.id.1km","lon.1km","lat.1km" ) %>% 
            tidyr::nest(MaxDispl.24h = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 
  
  # Dmax7d
  {if (!is.null(dmax7d)) left_join(.,  dmax7d %>%
                                     mutate(individual_id = as.character(individual_id)) %>%
                                     dplyr::select("individual_id","week","year_week","dmax7d","mean.x", "mean.y",
                                                   "grid.id.10km","lon.10km","lat.10km",
                                                   "grid.id.1km","lon.1km","lat.1km" ) %>% 
                tidyr::nest(MaxDispl.7d = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 

  # Dmax12m
  {if (!is.null(dmax12m)) left_join(.,  dmax12m %>%
                                      mutate(individual_id = as.character(individual_id)) %>%
                                      dplyr::select("individual_id","year","dmax12m","mean.x", "mean.y",
                                                    "grid.id.10km","lon.10km","lat.10km",
                                                    "grid.id.1km","lon.1km","lat.1km" ) %>% 
                tidyr::nest(MaxDispl.12m = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 
  
  # mcp.daily
  {if (!is.null(mcp.daily)) left_join(.,  mcp.daily %>%
                                        mutate(individual_id = as.character(individual_id)) %>%
                                        dplyr::select("individual_id","ymd","area","mean.x", "mean.y",
                                                      "grid.id.10km","lon.10km","lat.10km",
                                                      "grid.id.1km","lon.1km","lat.1km" ) %>% 
                tidyr::nest(Mcp.24h = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 
  
  # mcp.weekly
  {if (!is.null(mcp.weekly)) left_join(.,  mcp.weekly %>%
                                         mutate(individual_id = as.character(individual_id)) %>%
                                         dplyr::select("individual_id","week","year_week","area","mean.x", "mean.y",
                                                       "grid.id.10km","lon.10km","lat.10km",
                                                       "grid.id.1km","lon.1km","lat.1km" ) %>% 
                tidyr::nest(Mcp.7d = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>%
  
  # mcp.monthly
  {if (!is.null(mcp.monthly)) left_join(.,  mcp.monthly %>%
                                          mutate(individual_id = as.character(individual_id)) %>%
                                          dplyr::select("individual_id","month","year_month","area","mean.x", "mean.y",
                                                        "grid.id.10km","lon.10km","lat.10km",
                                                        "grid.id.1km","lon.1km","lat.1km" ) %>% 
                tidyr::nest(Mcp.1m = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 

  # mcp.annual
  {if (!is.null(mcp.annual)) left_join(.,  mcp.annual %>%
                                         mutate(individual_id = as.character(individual_id)) %>%
                                         dplyr::select("individual_id","year","area","mean.x", "mean.y",
                                                       "grid.id.10km","lon.10km","lat.10km",
                                                       "grid.id.1km","lon.1km","lat.1km" ) %>% 
                tidyr::nest(Mcp.12m = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 

  # Intensity of Use 24h
  {if (!is.null(df.IoU24h)) left_join(.,  df.IoU24h %>%
              ungroup() |>
                mutate(individual_id = as.character(individual_id)) %>%
                dplyr::select("individual_id","ymd","iou24h","mean.x", "mean.y",
                              "grid.id.10km","lon.10km","lat.10km",
                              "grid.id.1km","lon.1km","lat.1km") %>% 
                tidyr::nest(IoU.24h = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 

  # Intensity of Use 1m
  {if (!is.null(df.IoU1m)) left_join(.,  df.IoU1m %>%
              ungroup() |> 
                mutate(individual_id = as.character(individual_id)) %>%
                dplyr::select("individual_id","month","year","iou1m","mean.x", "mean.y",
                              "grid.id.10km","lon.10km","lat.10km",
                              "grid.id.1km","lon.1km","lat.1km" ) %>% 
                tidyr::nest(IoU.1m = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 

  # Intensity of Use 12m
  {if (!is.null(df.IoU12m)) left_join(.,  df.IoU12m %>%
              ungroup() |>
                mutate(individual_id = as.character(individual_id)) %>%
                dplyr::select("individual_id","year","iou12m","mean.x", "mean.y",
                              "grid.id.10km","lon.10km","lat.10km",
                              "grid.id.1km","lon.1km","lat.1km" ) %>% 
                tidyr::nest(IoU.12m = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 

  # Diurnality Index
  {if (!is.null(DI)) left_join(.,  DI %>%
                                 mutate(individual_id = as.character(individual_id)) %>%
                                 dplyr::select("individual_id","ymd","diurnality","mean.x","mean.y",
                                               "grid.id.10km","lon.10km","lat.10km",
                                               "grid.id.1km","lon.1km","lat.1km" ) %>% 
                tidyr::nest(Diurnality = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} }

    
## ----Save full database including raw metrics data--------------------
saveRDS(MoveTrait.v0.1_spatial, file=paste0(pthtrait,"3809257699_3809647332.rds"))


rm(animlocs.1hourly_sl);rm(animlocs.daily_sl);rm(dmax7d);rm(dmax12m);rm(mcp.daily);
rm(mcp.weekly);rm(mcp.monthly);rm(mcp.annual);rm(df.IoU24h);rm(df.IoU1m);
rm(df.IoU12m);rm(DI)

#})

