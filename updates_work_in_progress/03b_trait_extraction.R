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
flsMV
lapply(flsMV, function(indPth)
  {
  #indPth<-keep[1]
  animlocs.1hourly <- readRDS(file.path(pthamt1h, kp))
  
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
sum.monthly.ind.mcp1m <- f_sum.monthly.ind.mcp1m(mcp1m)

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
sum.monthly.ind.iou1m <- f_sum.monthly.ind.iou1m(iou1m) 

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
saveRDS(MoveTrait.ind.sum, file=paste0(pthtraitsum,indPth))

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
saveRDS(MoveTrait.monthly.ind.sum, file=paste0(pthtraitsummonthly,indPth))

## ----Build full database including raw metrics data--------------------

MoveTrait.repeats <- 
  if (is.null(MoveTrait.ind.sum)) NULL else {
    MoveTrait.ind.sum %>%
  dplyr::select(individual_id,study_id) |> 
  tidyr::nest(data = -individual_id) %>% 
  dplyr::select(-data) %>% # just a quick trick to keep things flowing.
  left_join(., MoveTrait.ind.sum[!duplicated(MoveTrait.ind.sum$individual_id),1:2], 
            by = c("individual_id" = "individual_id")) %>% 
  
  # 1 hourly
  {if (!is.null(d1h)) left_join(.,  d1h %>%
                                    data.frame %>%
                                    mutate(individual_id = as.character(individual_id)) %>%
                tidyr::nest(d1h = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 
  
  # 24 hourly
  {if (!is.null(d24h)) left_join(.,  d24h %>%
                                    mutate(individual_id = as.character(individual_id)) %>%
                tidyr::nest(d24h = -individual_id), 
            by = c("individual_id" = "individual_id")) else .}  %>% 
  
  # Dmax24
  {if (!is.null(dmax24h)) left_join(.,  dmax24h %>%
                                    mutate(individual_id = as.character(individual_id)) %>%
            tidyr::nest(dmax24h = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 
  
  # Dmax7d
  {if (!is.null(dmax7d)) left_join(.,  dmax7d %>%
                                     mutate(individual_id = as.character(individual_id)) %>%
                tidyr::nest(dmax7d = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 

  # Dmax12m
  {if (!is.null(dmax12m)) left_join(.,  dmax12m %>%
                                      mutate(individual_id = as.character(individual_id)) %>%
                tidyr::nest(dmax12m = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 
  
  # mcp.daily
  {if (!is.null(mcp24h)) left_join(.,  mcp24h %>%
                                        mutate(individual_id = as.character(individual_id)) %>%
                tidyr::nest(mcp24h = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 
  
  # mcp.weekly
  {if (!is.null(mcp7d)) left_join(.,  mcp7d %>%
                                         mutate(individual_id = as.character(individual_id)) %>%
                tidyr::nest(mcp7d = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>%
  
  # mcp.monthly
  {if (!is.null(mcp1m)) left_join(.,  mcp1m %>%
                                          mutate(individual_id = as.character(individual_id)) %>%
                tidyr::nest(mcp1m = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 

  # mcp.annual
  {if (!is.null(mcp12m)) left_join(.,  mcp12m %>%
                                         mutate(individual_id = as.character(individual_id)) %>%
                tidyr::nest(mcp12m = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 

  # Intensity of Use 24h
  {if (!is.null(iou24h)) left_join(.,  iou24h %>%
              ungroup() |>
                mutate(individual_id = as.character(individual_id)) %>%
                tidyr::nest(iou24h = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 

  # Intensity of Use 1m
  {if (!is.null(iou1m)) left_join(.,  iou1m %>%
              ungroup() |> 
                mutate(individual_id = as.character(individual_id)) %>%
                tidyr::nest(iou1m = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 

  # Intensity of Use 12m
  {if (!is.null(iou12m)) left_join(.,  iou12m %>%
              ungroup() |>
                mutate(individual_id = as.character(individual_id)) %>%
                tidyr::nest(iou12m = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} %>% 

  # Diurnality Index
  {if (!is.null(di)) left_join(.,  di %>%
                                 mutate(individual_id = as.character(individual_id)) %>%
                tidyr::nest(di = -individual_id), 
            by = c("individual_id" = "individual_id"))  else .} }

    
## ----Save full database including raw metrics data--------------------
saveRDS(MoveTrait.repeats, file=paste0(pthtrait,indPth))

# rm(d1h);rm(d24h);rm(dmax7d);rm(dmax12m);rm(mcp24h);
# rm(mcp7d);rm(mcp1m);rm(mcp12m);rm(iou24h);rm(iou1m);
# rm(iou12m);rm(di);rm(animlocs.1hourly);rm(animlocs.daily);rm(animlocs.weekly)

})

