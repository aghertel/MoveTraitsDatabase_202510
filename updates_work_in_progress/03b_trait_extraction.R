# ---
# title: "MoveTraits Database"
# author: "Anne G. Hertel"
# date: "28/1/2025"
# ---

library(lubridate);library(metafor);library(tidyverse);library(amt);
library(adehabitatHR); library(move2); library(epitools); library(suncalc); library(purrr); library(bit64)
library(dggridR);library(mapview);library(ISOweek);library(lubridate)


## ----Import functions for traits-------------------------------------------------------------
source(here::here('trait_scripts', 'd1h.R'))
source(here::here('trait_scripts', 'dmax24h.R'))

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

## ----Displacement traits-------------------------------------------------------------

## ----1h displacement-------------------------------------------------------------
d1h <- calc_d1h(animlocs.1hourly, dggs.10, dggs.1)
sum.ind.d1h <- f_sum.ind.d1h(d1h)
sum.monthly.ind.d1h <- f_sum.monthly.ind.d1h(d1h)

## ----Maximum 24hr displacement-------------------------------------------------------------
dmax24h <- calc_dmax24h(animlocs.1hourly, dggs.10, dggs.1)
sum.ind.dmax24h <- f_sum.ind.dmax24h(dmax24h)
sum.monthly.ind.dmax24h <- f_sum.monthly.ind.dmax24h(dmax24h)

## ----24hr displacement distance-------------------------------------------------------------
d24h <- calc_d24h(animlocs.daily, dggs.10, dggs.1)
sum.ind.d24h <- f_sum.ind.d24h(d24h)
sum.monthly.ind.d24h <- f_sum.monthly.ind.d24h(d24h)

## ----Maximum 7day displacement distance-------------------------------------------------------------
dmax7d <- calc_dmax7d(animlocs.daily, dggs.10, dggs.1)
sum.ind.dmax7d <- f_sum.ind.dmax7d(dmax7d)
sum.monthly.ind.dmax7d <- f_sum.monthly.ind.dmax7d(dmax7d)

## ----Maximum annual displacement distance-------------------------------------------------------------
dmax12m <- calc_dmax12m(animlocs.weekly, dggs.10, dggs.1)
sum.ind.dmax12m <- f_sum.ind.dmax12m(dmax12m)
#sum.monthly.ind.dmax12m<- f_sum.monthly.ind.dmax12m(dmax12m)
## ???

## ----Range size traits-------------------------------------------------------------
## ----Daily MCP-------------------------------------------------------------
mcp24h <- calc_mcp24h(animlocs.1hourly, dggs.10, dggs.1)
sum.ind.mcp24h <- f_sum.ind.mcp24h(mcp24h)
sum.monthly.ind.mcp24h <- f_sum.monthly.ind.mcp24h(mcp24h)

## ----Weekly MCP-------------------------------------------------------------
mcp7d <- calc_mcp7d(animlocs.1hourly, dggs.10, dggs.1)
sum.ind.mcp7d <- f_sum.ind.mcp7d(mcp7d)
sum.monthly.ind.mcp7d <- f_sum.monthly.ind.mcp24h(mcp7d)

## ----Monthly MCP-------------------------------------------------------------
mcp1m <- calc_mcp1m(animlocs.daily, dggs.10, dggs.1)
sum.ind.mcp1m <- f_sum.ind.mcp1m(mcp1m)
sum.monthly.ind.mcp1m <- mcp1m[,1:4] |> 
  rename(mcp1m = area) |> 
  mutate(year = as.numeric(substr(year_month, 1, 4))) |> 
  dplyr::select(individual_id,month,year,mcp1m)

## ----Annual MCP-------------------------------------------------------------
mcp12m <- calc_mcp12m(animlocs.weekly, dggs.10, dggs.1)
sum.ind.mcp12m <- f_sum.ind.mcp12m(mcp12m)

## ----Intensity of use-------------------------------------------------------------

## ----Daily IOU-------------------------------------------------------------

mcp.daily <- 
  if(is.null(mcp.daily)) NULL else {
    mcp.daily|>
      mutate(id_ymd = paste(mcp.daily$individual_id,mcp.daily$ymd,sep="."))
    }

df.IoU24h <- 
  if(is.null(animlocs.1hourly_sl) | is.null(mcp.daily)) NULL else {
  animlocs.1hourly_sl |>  
  group_by(ymd) %>% 
  mutate(cumsumD1h = sum(d1h,na.rm=T),
         mean.x = mean(x_),
         mean.y = mean(y_)) %>% 
  dplyr::select(individual_id,ymd,cumsumD1h,mean.x,mean.y) %>% distinct() %>%
  left_join(mcp.daily[,c("ymd","area")],by = "ymd") %>% 
  mutate(iou24h = cumsumD1h/sqrt(area)) %>% 
  filter(!is.na(iou24h))
    }

if(is.null(df.IoU24h)) NULL else {
  # Spatial annotation 10km
  cell_info <- dgGEO_to_SEQNUM(dggs.10, df.IoU24h$mean.x, df.IoU24h$mean.y)
  df.IoU24h$grid.id.10km <- cell_info$seqnum
  centers <- dgSEQNUM_to_GEO(dggs.10, df.IoU24h$grid.id.10km)
  df.IoU24h$lon.10km <- centers$lon_deg
  df.IoU24h$lat.10km <- centers$lat_deg
  rm(cell_info);rm(centers)
  
  # Spatial annotation 1km
  cell_info <- dgGEO_to_SEQNUM(dggs.1, df.IoU24h$mean.x, df.IoU24h$mean.y)
  df.IoU24h$grid.id.1km <- cell_info$seqnum
  centers <- dgSEQNUM_to_GEO(dggs.1, df.IoU24h$grid.id.1km)
  df.IoU24h$lon.1km <- centers$lon_deg
  df.IoU24h$lat.1km <- centers$lat_deg
  rm(cell_info);rm(centers)
}

## ----Monthly IOU-------------------------------------------------------------
mcp.monthly <- 
  if(is.null(mcp.monthly)) NULL else {
    mcp.monthly|>
      mutate(id_month = paste(mcp.monthly$individual_id,mcp.monthly$year_month,sep="."))
  }

df.IoU1m <- 
  if(is.null(animlocs.daily_sl) | is.null(mcp.monthly)) NULL else {
  animlocs.daily_sl  %>%  
  mutate(year_month = paste(year,month,sep="_"),
         id_month = paste(individual_id,year_month,sep="."))  |>  group_by(id_month)  |>  
  mutate(cumsum.d24h = sum(d24h,na.rm=T),
         mean.x = mean(x_),
         mean.y = mean(y_))  |>  
  dplyr::select(individual_id,id_month,month,year,year_month,cumsum.d24h,mean.x,mean.y)  |>
  distinct() |> 
  left_join(mcp.monthly[,c("id_month","area")],by = "id_month")  |>  
  mutate(iou1m = cumsum.d24h/sqrt(area))  |>  
  filter(!is.na(iou1m)) |> ungroup() |> 
  dplyr::select(individual_id,month,year,year_month,mean.x,mean.y,cumsum.d24h,area,iou1m)  
  }

if(is.null(df.IoU1m)) NULL else {
  # Spatial annotation 10km
  cell_info <- dgGEO_to_SEQNUM(dggs.10, df.IoU1m$mean.x, df.IoU1m$mean.y)
  df.IoU1m$grid.id.10km <- cell_info$seqnum
  centers <- dgSEQNUM_to_GEO(dggs.10, df.IoU1m$grid.id.10km)
  df.IoU1m$lon.10km <- centers$lon_deg
  df.IoU1m$lat.10km <- centers$lat_deg
  rm(cell_info);rm(centers)
  
  # Spatial annotation 1km
  cell_info <- dgGEO_to_SEQNUM(dggs.1, df.IoU1m$mean.x, df.IoU1m$mean.y)
  df.IoU1m$grid.id.1km <- cell_info$seqnum
  centers <- dgSEQNUM_to_GEO(dggs.1, df.IoU1m$grid.id.1km)
  df.IoU1m$lon.1km <- centers$lon_deg
  df.IoU1m$lat.1km <- centers$lat_deg
  rm(cell_info);rm(centers)
}
  
## ----Annual IOU-------------------------------------------------------------
mcp.annual <- 
  if(is.null(mcp.annual)) NULL else {
    mcp.annual|>
      mutate(id_year = paste(mcp.annual$individual_id,mcp.annual$year,sep="."))
  }

df.IoU12m <- 
  if(is.null(animlocs.daily_sl) | is.null(mcp.annual)) NULL else {
  animlocs.daily_sl  |>   
  mutate(id_year = paste(individual_id,year,sep="."))  |>  group_by(id_year) |>  
  mutate(cumsum.d24h = sum(d24h,na.rm=T),
         mean.x = mean(x_),
         mean.y = mean(y_)) |>  
  dplyr::select(id_year,year,individual_id,cumsum.d24h,mean.x,mean.y) |>  distinct() |> 
  left_join(mcp.annual[,c("id_year","area")],by = "id_year") |>  
  mutate(iou12m = cumsum.d24h/sqrt(area)) |> 
  filter(!is.na(iou12m)) 
  }

if(is.null(df.IoU12m)) NULL else {
  # Spatial annotation 10km
  cell_info <- dgGEO_to_SEQNUM(dggs.10, df.IoU12m$mean.x, df.IoU12m$mean.y)
  df.IoU12m$grid.id.10km <- cell_info$seqnum
  centers <- dgSEQNUM_to_GEO(dggs.10, df.IoU12m$grid.id.10km)
  df.IoU12m$lon.10km <- centers$lon_deg
  df.IoU12m$lat.10km <- centers$lat_deg
  rm(cell_info);rm(centers)
  
  # Spatial annotation 1km
  cell_info <- dgGEO_to_SEQNUM(dggs.1, df.IoU12m$mean.x, df.IoU12m$mean.y)
  df.IoU12m$grid.id.1km <- cell_info$seqnum
  centers <- dgSEQNUM_to_GEO(dggs.1, df.IoU12m$grid.id.1km)
  df.IoU12m$lon.1km <- centers$lon_deg
  df.IoU12m$lat.1km <- centers$lat_deg
  rm(cell_info);rm(centers)
}

## ----Diurnality Index-------------------------------------------------------------
#' Estimating proportiona daily diurnality, corrected for daylight changes,
#' using movement data. The diurnality index (from here on diurnality)
#' is based on Hoogenboom, Daan, Dallinga, and Schoenmakers (1984):
#' 
#' [(AD/DD) - (AN/DN)] / [(AD/DD) + (AN/DN)]
#' 
#' where AD and AN are the sums of the movement during the
#' day and night, respectively, and DD and DN are the durations of the
#' day and night, respectively. The diurnality index varies between -1
#' (night active) and 1 (day active).

# animlocs.1hourly_sl.t <- animlocs.1hourly_sl
# animlocs.1hourly_sl.t <- NULL

f.diurn <- function(a, b, c, d) {((a / b) - (c / d)) / ((a / b) + (c / d))}


DI <- 
  if(is.null(animlocs.1hourly_sl) ) NULL else {
    
diurn.ind <- animlocs.1hourly_sl |> 
  mutate(lat=y_,lon=x_,date=as.Date(t_)) |> 
  dplyr::select(lat,lon,date,individual_id,t_,ymd,d1h) 

sun_all <- getSunlightTimes(data = diurn.ind, 
                            tz = "UTC",
                            keep = c("sunrise", "sunset"))

diurn.ind$sunrise <- sun_all$sunrise
diurn.ind$sunset<- sun_all$sunset

diurn.ind$daytime <- 
  ifelse(diurn.ind$t_ < diurn.ind$sunrise |
         diurn.ind$t_ > diurn.ind$sunset, "night","day")

      mean.coord <- diurn.ind |> 
        mutate(id.day = paste(individual_id,ymd,sep=".")) |> group_by(id.day) |> 
        mutate(mean.x = mean(lon, na.rm=T),
               mean.y = mean(lat, na.rm=T)) |> 
        dplyr::select(id.day, mean.x, mean.y) |> distinct()

DI <- diurn.ind %>%
  data.frame %>% 
      group_by(individual_id, ymd)%>%
      summarise(dist.sum.day = sum(d1h[daytime=="day"]),
                dist.sum.night = sum(d1h[daytime=="night"]),
                dist.length.day = sum(daytime=="day"),
                dist.length.night = sum(daytime=="night"),
                total.daylength = dist.length.day + dist.length.night) |> 
  data.frame()

DI$diurnality <- f.diurn(DI$dist.sum.day, DI$dist.length.day, DI$dist.sum.night, DI$dist.length.night)

DI <- DI |> 
  mutate(id.day = paste(individual_id,ymd,sep=".")) |> 
  left_join(mean.coord,by = "id.day") |> 
  filter(total.daylength > 19)  |> 
  filter(!is.na(diurnality)) 

DI <- if(nrow(DI)==0) NULL else DI
}

if(is.null(DI)) NULL else {
  # Spatial annotation 10km
  cell_info <- dgGEO_to_SEQNUM(dggs.10, DI$mean.x, DI$mean.y)
  DI$grid.id.10km <- cell_info$seqnum
  centers <- dgSEQNUM_to_GEO(dggs.10, DI$grid.id.10km)
  DI$lon.10km <- centers$lon_deg
  DI$lat.10km <- centers$lat_deg
  rm(cell_info);rm(centers)
  
  # Spatial annotation 1km
  cell_info <- dgGEO_to_SEQNUM(dggs.1, DI$mean.x, DI$mean.y)
  DI$grid.id.1km <- cell_info$seqnum
  centers <- dgSEQNUM_to_GEO(dggs.1, DI$grid.id.1km)
  DI$lon.1km <- centers$lon_deg
  DI$lat.1km <- centers$lat_deg
  rm(cell_info);rm(centers)
}

### !!! ###





## ----function to summarize IoU24h--------------------
FSumIOU24h<-function(x)
  {
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n.iou24h.days = NA,
                       iou24h.mean = NA,iou24h.median = NA,iou24h.cv = NA,iou24h.95 = NA,iou24h.05 = NA)
  } else {
    
  individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
  
  # Get sample size per indivindividual_idual
  n.iou24h.days<-as.numeric(with(x, tapply(x$ymd,individual_id, length)))
  
  # 24hr Displacement
  iou24h.mean<-as.numeric(with(x, tapply(x$iou24h+0.001,individual_id, mean, na.rm=T)))
  iou24h.median<-as.numeric(with(x, tapply(x$iou24h+0.001,individual_id, median, na.rm=T)))
  iou24h.cv<-as.numeric(with(x, tapply(x$iou24h+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
  iou24h.95<-as.numeric(with(x, tapply(x$iou24h+0.001,individual_id, quantile,.95, na.rm=T)))
  iou24h.05<-as.numeric(with(x, tapply(x$iou24h+0.001,individual_id, quantile,.05, na.rm=T)))
  
  # build dataframe
  dats<-data.frame(individual_id,n.iou24h.days,
                   iou24h.mean,iou24h.median,iou24h.cv,iou24h.95,iou24h.05)
  
  return(dats)
  }}

IOU24h <- FSumIOU24h(df.IoU24h)

## ----function to summarize IoU1m--------------------
FSumIOU1m<-function(x)
  {
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n.iou1m.month = NA,
                       iou1m.mean = NA,iou1m.median = NA,iou1m.cv = NA,iou1m.95 = NA,iou1m.05 = NA)
  } else {
    
  individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
  
  # Get sample size per indivindividual_idual
  n.iou1m.month<-as.numeric(with(x, tapply(x$year_month,individual_id, length)))
  
  # 24hr Displacement
  iou1m.mean<-as.numeric(with(x, tapply(x$iou1m+0.001,individual_id, mean, na.rm=T)))
  iou1m.median<-as.numeric(with(x, tapply(x$iou1m+0.001,individual_id, median, na.rm=T)))
  iou1m.cv<-as.numeric(with(x, tapply(x$iou1m+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
  iou1m.95<-as.numeric(with(x, tapply(x$iou1m+0.001,individual_id, quantile,.95, na.rm=T)))
  iou1m.05<-as.numeric(with(x, tapply(x$iou1m+0.001,individual_id, quantile,.05, na.rm=T)))
  
  # build dataframe
  dats<-data.frame(individual_id,n.iou1m.month,
                   iou1m.mean,iou1m.median,iou1m.cv,iou1m.95,iou1m.05)
  
  return(dats)
}}

IOU1m <- FSumIOU1m(df.IoU1m)

## ----function to summarize IoU12m--------------------

FSumIOU12m<-function(x)
  {
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n.iou12m.year = NA,
                       iou12m.mean = NA,iou12m.median = NA,iou12m.cv = NA,iou12m.95 = NA,iou12m.05 = NA)
  } else {
    
  individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
  
  # Get sample size per indivindividual_idual
  n.iou12m.year<-as.numeric(with(x, tapply(x$year,individual_id, length)))
  
  # 24hr Displacement
  iou12m.mean<-as.numeric(with(x, tapply(x$iou12m+0.001,individual_id, mean, na.rm=T)))
  iou12m.median<-as.numeric(with(x, tapply(x$iou12m+0.001,individual_id, median, na.rm=T)))
  iou12m.cv<-as.numeric(with(x, tapply(x$iou12m+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
  iou12m.95<-as.numeric(with(x, tapply(x$iou12m+0.001,individual_id, quantile,.95, na.rm=T)))
  iou12m.05<-as.numeric(with(x, tapply(x$iou12m+0.001,individual_id, quantile,.05, na.rm=T)))
  
  # build dataframe
  dats<-data.frame(individual_id,n.iou12m.year,
                   iou12m.mean,iou12m.median,iou12m.cv,iou12m.95,iou12m.05)
  
  return(dats)
  }}

IOU12m <- FSumIOU12m(df.IoU12m)

## ----function to summarize Diurnality Index--------------------
FSumDI<-function(x)
  {
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n.di.days = NA,
                       di.mean = NA,di.median = NA,di.cv = NA,di.95 = NA,di.05 = NA)
  } else {
    
  individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
  
  # How many days per indivindividual_idual
  n.di.days<-as.numeric(with(x, tapply(x$ymd,individual_id, length)))
  
  # Diurnality
  di.mean<-as.numeric(with(x, tapply(x$diurnality,individual_id, mean, na.rm=T)))
  di.median<-as.numeric(with(x, tapply(x$diurnality,individual_id, median, na.rm=T)))
  di.cv<-as.numeric(with(x, tapply(x$diurnality,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
  di.95<-as.numeric(with(x, tapply(x$diurnality,individual_id, quantile,.95, na.rm=T)))
  di.05<-as.numeric(with(x, tapply(x$diurnality,individual_id, quantile,.05, na.rm=T)))
  
  # build dataframe
  dats<-data.frame(individual_id,n.di.days,
                   di.mean,di.median,di.cv,di.95,di.05)
  
  return(dats)
  }
  }

DI.12 <- FSumDI(DI)

## ----Build database with summary values--------------------
#cbind individual summaries of all traits 
MoveTrait.v0.1 <- bind_cols(Displ1h, 
                          bind_cols(Displ24h[,2:length(Displ24h)], 
                          bind_cols(MaxDispl24h[,2:length(MaxDispl24h)], 
                          bind_cols(MaxDispl7d[,2:length(MaxDispl7d)], 
                          bind_cols(MaxDispl12m[,2:length(MaxDispl12m)], 
                          bind_cols(MCP24h[,2:length(MCP24h)], 
                          bind_cols(MCP7d[,2:length(MCP7d)], 
                          bind_cols(MCP1m[,2:length(MCP1m)], 
                          bind_cols(MCP12m[,2:length(MCP12m)], 
                          bind_cols(IOU24h[,2:length(IOU24h)],
                          bind_cols(IOU1m[,2:length(IOU1m)],
                          bind_cols(IOU12m[,2:length(IOU12m)], DI.12[,2:length(DI.12)])))))))))))) %>% 
  filter(if_any(everything(), ~ !is.na(.)))

#mean coordinates for individual summaries
library(bit64)
movedata2 <- animlocs.1hourly %>% 
  mutate(study_id = as.integer64(study_id),
         individual_id = individual_id) |> 
  group_by(individual_id) %>% 
  mutate(mean.longitude = mean(x_,na.rm=T),
         mean.latitude = mean(y_,na.rm=T)) %>% 
  dplyr::select(study_id,individual_id,mean.longitude,mean.latitude) %>% 
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

# join database
MoveTrait.v0.1 <- 
  if (is.null(MoveTrait.v0.1) | nrow(MoveTrait.v0.1) == 0) NULL else {
    movedata2 %>% 
    left_join(MoveTrait.v0.1, by = "individual_id") %>% 
    droplevels() }

## ----Save database with summaries--------------------
saveRDS(MoveTrait.v0.1, file=paste0(pthtraitsum,"3809257699_3809647332.rds"))

# #coordinates for monthly individual summaries
# id_monthly <- animlocs.1hourly |> 
#   mutate(month = month(t_),
#          year = year(t_),
#          month_year = paste(month,year,sep="_"))
# 
# cell_info.100 <- dgGEO_to_SEQNUM(dggs.100, animlocs.1hourly$x_, animlocs.1hourly$y_)
# id_monthly$grid.id.100km <- cell_info.100$seqnum
# id_monthly_100km <- id_monthly |> 
#   distinct(month_year,grid.id.100km,.keep_all = TRUE)
# 
# cell_info.10 <- dgGEO_to_SEQNUM(dggs.10, animlocs.1hourly$x_, animlocs.1hourly$y_)
# id_monthly$grid.id.10km <- cell_info.10$seqnum
# id_monthly_10km <- id_monthly |> 
#   distinct(month_year,grid.id.10km,.keep_all = TRUE)

#pthtraitsummonthly

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
                              "grid.id.1km","lon.1km","lat.1km" ) %>% 
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

