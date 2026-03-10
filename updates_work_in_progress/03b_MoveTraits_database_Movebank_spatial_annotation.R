# ---
# title: "MoveTraits Database"
# author: "Anne G. Hertel"
# date: "28/1/2025"
# ---

library(lubridate);library(metafor);library(tidyverse);library(amt);
library(adehabitatHR); library(move2); library(epitools); library(suncalc); library(purrr); library(bit64)
library(dggridR);library(mapview)
## ----Import movement data per individual-------------------------------------------------------------
pathTOfolder <- "./StandardGrid/DATA/"

pthamt1h <- paste0(pathTOfolder,"4.MB_indv_amt_1h/")

#dir for individual summaries
dir.create(paste0(pathTOfolder,"5.MB_indv_traitsum"))
pthtraitsum <- paste0(pathTOfolder,"5.MB_indv_traitsum/")
#dir for individual underlying traits
dir.create(paste0(pathTOfolder,"6.MB_indv_trait"))
pthtrait <- paste0(pathTOfolder,"6.MB_indv_trait/")

#spatial grid
dggs.100    <- dgconstruct(projection = "ISEA", area = 10000, resround='nearest')
dggs.10     <- dgconstruct(projection = "ISEA", area = 100, resround='nearest')
dggs.1      <- dgconstruct(projection = "ISEA", area = 1, resround='nearest')

# list files
flsMV <- list.files(pthamt1h, full.names = F)
done <- list.files(pthtraitsum, full.names = F)
flsMV <- flsMV[!flsMV%in%done]

#filter excluded studies
referenceTableStudies <- readRDS(paste0(pathTOfolder,"/referenceTableStudies_ALL_excludedColumn.rds"))
referenceTableStudiesUsed <- referenceTableStudies[referenceTableStudies$excluded=="no",]

flsMV <- flsMV[flsMV %in% referenceTableStudiesUsed$fileName]

#lapply(flsMV, function(indPth)
#  {
animlocs.1hourly <- readRDS(paste0(pthamt1h,"3809257699_3809647332.rds")) #  #446579_7945602.rds

  #coordinates for individual summaries
  cell_info.100 <- dgGEO_to_SEQNUM(dggs.100, animlocs.1hourly$x_, animlocs.1hourly$y_)
  grid.id.100km <- unique(cell_info.100$seqnum)
  lon.100km <- dgSEQNUM_to_GEO(dggs.100, grid.id.100km)$lon_deg
  lat.100km <- dgSEQNUM_to_GEO(dggs.100, grid.id.100km)$lat_deg

  cell_info.10 <- dgGEO_to_SEQNUM(dggs.10, animlocs.1hourly$x_, animlocs.1hourly$y_)
  grid.id.10km <- unique(cell_info.10$seqnum)
  lon.10km <- dgSEQNUM_to_GEO(dggs.10, grid.id.10km)$lon_deg
  lat.10km <- dgSEQNUM_to_GEO(dggs.10, grid.id.10km)$lat_deg

#Visualization
# pts_sf <- sf::st_as_sf(animlocs.1hourly, coords = c("x_", "y_"), crs = 4326)
# intersect_cells_10 <- dgcellstogrid(dggs.10, grid.id.10km)
# intersect_cells_100 <- dgcellstogrid(dggs.100, grid.id.100km)
# mapview(intersect_cells_100, map.types = c("OpenStreetMap.DE"), col.regions = "pink")+
#   mapview(intersect_cells_10, alpha = 1, col.regions = "purple")+
#   mapview(pts_sf, col.regions = "black",cex = 0.2)

  
  #coordinates for monthly individual summaries
  id_monthly <- animlocs.1hourly |> 
    mutate(month = month(t_),
           year = year(t_),
           month_year = paste(month,year,sep="_"))
  
  cell_info.100 <- dgGEO_to_SEQNUM(dggs.100, animlocs.1hourly$x_, animlocs.1hourly$y_)
  id_monthly$grid.id.100km <- cell_info.100$seqnum
  id_monthly_100km <- id_monthly |> 
    distinct(month_year,grid.id.100km,.keep_all = TRUE)
  
  cell_info.10 <- dgGEO_to_SEQNUM(dggs.10, animlocs.1hourly$x_, animlocs.1hourly$y_)
  id_monthly$grid.id.10km <- cell_info.10$seqnum
  id_monthly_10km <- id_monthly |> 
    distinct(month_year,grid.id.10km,.keep_all = TRUE)
  
  #Visualization
  # intersect_cells_100 <- dgcellstogrid(dggs.100, id_monthly$grid.id.100km)
  # intersect_cells_100 <- merge(intersect_cells_100,id_monthly[,c("grid.id.100km","month_year")],by.x="seqnum",by.y="grid.id.100km")
  # intersect_cells_100$month_year <- factor(intersect_cells_100$month_year)
  # 
  # intersect_cells_10 <- dgcellstogrid(dggs.10, id_monthly$grid.id.10km)
  # intersect_cells_10 <- merge(intersect_cells_10,id_monthly[,c("grid.id.10km","month_year")],by.x="seqnum",by.y="grid.id.10km")
  # intersect_cells_10$month_year <- factor(intersect_cells_10$month_year)
  
  #mapview(intersect_cells_10, map.types = c("OpenStreetMap.DE"), zcol = "month_year")
  # 
  # ggplot() +
  #   geom_sf(data = intersect_cells_100, fill = NA, colour = "black") +
  #   theme_classic() +
  #   facet_wrap(~month_year) +
  #   theme(axis.text = element_blank())
  # 
  # ggsave("StandardGrid/visualization/panthera leo/monthly_sum_100km.png")
  
## ----Resample data-------------------------------------------------------------
#Resample data to 24h, 7 week time scales using amt

animlocs.daily <- animlocs.1hourly |>
  track_resample(rate = hours(24),
               tolerance = minutes(60))

animlocs.weekly <- animlocs.1hourly |>
  track_resample(rate = hours(24*7),
                 tolerance = minutes(60*24))


data_resampled <- animlocs.1hourly %>%
  tidyr::nest(data = -individual_id) %>% 
  dplyr::select(-data) %>% # just a quick trick to keep things flowing.
  left_join(., animlocs.1hourly[!duplicated(animlocs.1hourly$individual_id),
                        c("study_id","individual_id")], 
            by = c("individual_id" = "individual_id")) %>% 
  # 1 hourly
  left_join(.,  animlocs.1hourly %>%
              mutate(id = individual_id) %>% 
              tidyr::nest(data = -individual_id), 
            by = c("individual_id" = "individual_id")) %>% 
  dplyr::rename(animlocs.1hourly = data) %>% 
  # daily 
  left_join(.,  animlocs.daily %>%
            mutate(id = individual_id) %>% 
            tidyr::nest(data = -individual_id),
          by = c("individual_id" = "individual_id")) %>%
  dplyr::rename(animlocs.daily = data) %>%
  # weekly 
  left_join(.,  animlocs.weekly %>%
            mutate(id = individual_id) %>% 
            tidyr::nest(data = -individual_id),
          by = c("individual_id" = "individual_id")) %>%
  dplyr::rename(animlocs.weekly = data)

## ----Movement metrics-------------------------------------------------------------
## ----1h displacement-------------------------------------------------------------
animlocs.1hourly_sl <- 
  animlocs.1hourly %>% 
  mutate(d1h = step_lengths(.)*100000) |> 
  mutate(time_diff = as.numeric(difftime(lead(t_),t_, units = "mins"))) |> 
  filter(time_diff >= 45 & time_diff <= 75) |> 
  dplyr::select(individual_id, t_, d1h, x_, y_) |> 
  mutate(ymd = as.character(format(as.Date(t_), "%Y-%m-%d"))) |> 
  filter(!is.na(d1h))

animlocs.1hourly_sl <-  if(nrow(animlocs.1hourly_sl)<167) NULL else {animlocs.1hourly_sl}

if(is.null(animlocs.1hourly_sl)) NULL else {
# Spatial annotation 10km
cell_info <- dgGEO_to_SEQNUM(dggs.10, animlocs.1hourly_sl$x_, animlocs.1hourly_sl$y_)
animlocs.1hourly_sl$grid.id.10km <- cell_info$seqnum
centers <- dgSEQNUM_to_GEO(dggs.10, animlocs.1hourly_sl$grid.id.10km)
animlocs.1hourly_sl$lon.10km <- centers$lon_deg
animlocs.1hourly_sl$lat.10km <- centers$lat_deg
rm(cell_info);rm(centers)

# Spatial annotation 1km
cell_info <- dgGEO_to_SEQNUM(dggs.1, animlocs.1hourly_sl$x_, animlocs.1hourly_sl$y_)
animlocs.1hourly_sl$grid.id.1km <- cell_info$seqnum
centers <- dgSEQNUM_to_GEO(dggs.1, animlocs.1hourly_sl$grid.id.1km)
animlocs.1hourly_sl$lon.1km <- centers$lon_deg
animlocs.1hourly_sl$lat.1km <- centers$lat_deg
rm(cell_info);rm(centers)
}

## ----Maximum 24hr displacement-------------------------------------------------------------
locs1h <- animlocs.1hourly %>% 
  mutate(ymd = as.character(format(as.Date(t_), "%Y-%m-%d"))) |> 
  mutate(id.ymd = paste(individual_id,ymd,sep="_")) |> 
  group_by(id.ymd) %>% 
  filter(n() >= 12) %>%
  ungroup() 

mean.coord <- locs1h |> group_by(id.ymd) |> 
  mutate(mean.x = mean(x_, na.rm=T),
         mean.y = mean(y_, na.rm=T)) |> 
  dplyr::select(id.ymd, mean.x, mean.y) |> distinct()

locs1h.sf <- sf::st_as_sf(locs1h,
                          coords = c("x_", "y_"),
                          crs = 4326)

moveObjSplitTime <- split(locs1h.sf, locs1h$id.ymd)
maxNetDispL <- lapply(moveObjSplitTime, function(x){max(sf::st_distance(x))})
maxNetDisp <- do.call("rbind",maxNetDispL)

rm(moveObjSplitTime);rm(maxNetDispL);rm(locs1h)

dmax24 <- 
  if(is.null(maxNetDisp)) NULL else {
    data.frame(keyName=row.names(maxNetDisp), dmax24h=maxNetDisp[,1], row.names=NULL) |> 
      mutate(ymd = str_split(keyName, "_", simplify = TRUE)[,2],
             individual_id = str_split(keyName, "_", simplify = TRUE)[,1]) |> 
  filter(!is.na(dmax24h)) |>  
  group_by(individual_id) |>  
  filter(n() >= 7) %>%
  ungroup() |> 
  mutate(id.ymd = paste(individual_id,ymd,sep="_")) |> 
  left_join(mean.coord, by = "id.ymd") |> 
  dplyr::select(individual_id,ymd,dmax24h, mean.x, mean.y) }

dmax24 <- if(is.null(dmax24)) NULL else 
  if(nrow(dmax24) == 0) NULL else dmax24 

if(is.null(dmax24)) NULL else {
# Spatial annotation 10km
cell_info <- dgGEO_to_SEQNUM(dggs.10, mean.coord$mean.x, mean.coord$mean.y)
dmax24$grid.id.10km <- cell_info$seqnum
centers <- dgSEQNUM_to_GEO(dggs.10, dmax24$grid.id.10km)
dmax24$lon.10km <- centers$lon_deg
dmax24$lat.10km <- centers$lat_deg
rm(cell_info);rm(centers)

# Spatial annotation 1km
cell_info <- dgGEO_to_SEQNUM(dggs.1, mean.coord$mean.x, mean.coord$mean.y)
dmax24$grid.id.1km <- cell_info$seqnum
centers <- dgSEQNUM_to_GEO(dggs.1, dmax24$grid.id.1km)
dmax24$lon.1km <- centers$lon_deg
dmax24$lat.1km <- centers$lat_deg
rm(cell_info);rm(centers)
}

# # # Visualization
# pts_sf <- locs1h.sf[73:96,]
# pts_sf2 <- mean.coord |> filter(id.ymd %in% "3809647332_2009-08-01") |> sf::st_as_sf(coords = c("mean.x", "mean.y"))
# tmp1 <- dmax24 |> filter(ymd %in% "2009-08-01") |> dplyr::select("grid.id.10km")
# tmp2 <- dmax24 |> filter(ymd %in% "2009-08-01") |> dplyr::select("grid.id.1km")
# intersect_cells_10 <- dgcellstogrid(dggs.10, tmp1[[1]])
# intersect_cells_1 <- dgcellstogrid(dggs.1, tmp2[[1]])
# mapview(intersect_cells_10, map.types = c("OpenStreetMap.DE"), col.regions = "red")+
#  mapview(intersect_cells_1, alpha = 1, col.regions = "darkred")+
#  mapview(pts_sf, col.regions = "black")+
#  mapview(pts_sf2, col.regions = "yellow")

rm(mean.coord);rm(locs1h.sf)

## ----24hr displacement distance-------------------------------------------------------------
animlocs.daily_sl <- animlocs.daily %>% 
  mutate(d24h = step_lengths(.)*100000) |> 
  mutate(time_diff = as.numeric(difftime(lead(t_),t_, units = "hours"))) |> 
  filter(time_diff >= 20 & time_diff <= 28) |> 
  dplyr::select(individual_id, t_, d24h, x_, y_) |> 
  mutate(ymd = as.character(format(as.Date(t_), "%Y-%m-%d")),
         month = lubridate::month(t_),
         year = lubridate::year(t_))|> 
  filter(!is.na(d24h))

animlocs.daily_sl <-  if(nrow(animlocs.daily_sl)<30) NULL else {animlocs.daily_sl}

if(is.null(animlocs.daily_sl)) NULL else {
# Spatial annotation 10km
cell_info <- dgGEO_to_SEQNUM(dggs.10, animlocs.daily_sl$x_, animlocs.daily_sl$y_)
animlocs.daily_sl$grid.id.10km <- cell_info$seqnum
centers <- dgSEQNUM_to_GEO(dggs.10, animlocs.daily_sl$grid.id.10km)
animlocs.daily_sl$lon.10km <- centers$lon_deg
animlocs.daily_sl$lat.10km <- centers$lat_deg
rm(cell_info);rm(centers)

# Spatial annotation 1km
cell_info <- dgGEO_to_SEQNUM(dggs.1, animlocs.daily_sl$x_, animlocs.daily_sl$y_)
animlocs.daily_sl$grid.id.1km <- cell_info$seqnum
centers <- dgSEQNUM_to_GEO(dggs.1, animlocs.daily_sl$grid.id.1km)
animlocs.daily_sl$lon.1km <- centers$lon_deg
animlocs.daily_sl$lat.1km <- centers$lat_deg
rm(cell_info);rm(centers)
}

 # # Visualization
 #  pts_sf <- sf::st_as_sf(animlocs.daily_sl,
 #                           coords = c("x_", "y_"),
 #                           crs = 4326)
 # pts_sf <- pts_sf[5,]
 # tmp1 <- animlocs.daily_sl |> filter(ymd %in% "2009-08-01") |> dplyr::select("grid.id.10km")
 # tmp2 <- animlocs.daily_sl |> filter(ymd %in% "2009-08-01") |> dplyr::select("grid.id.1km")
 # intersect_cells_10 <- dgcellstogrid(dggs.10, tmp1[[1]])
 # intersect_cells_1 <- dgcellstogrid(dggs.1, tmp2[[1]])
 # mapview(intersect_cells_10, map.types = c("OpenStreetMap.DE"), col.regions = "red")+
 #  mapview(intersect_cells_1, alpha = 1, col.regions = "darkred")+
 #  mapview(pts_sf, col.regions = "black")
 
## ----Maximum 7day displacement distance-------------------------------------------------------------
locs24h <- flatten(data_resampled[,"animlocs.daily"]) |> bind_rows() |> tibble() |> 
  mutate(week = as.numeric(strftime(t_,format="%W")), year = as.numeric(strftime(t_,format="%Y")),
         year_week = paste(year,week, sep="_")) |> 
  mutate(id.year_week = paste(id,year_week,sep="_")) |>  
  group_by(id,year_week) |> filter(n() >= 5) |> ungroup()

mean.coord <- locs24h |> group_by(id.year_week) |> 
  mutate(mean.x = mean(x_, na.rm=T),
         mean.y = mean(y_, na.rm=T)) |> 
  dplyr::select(id.year_week, mean.x, mean.y) |> distinct()

locs24h.sf <- sf::st_as_sf(locs24h,
                          coords = c("x_", "y_"),
                          crs = 4326)

moveObjSplitTime <- split(locs24h.sf, locs24h.sf$id.year_week)
maxNetDispL <- lapply(moveObjSplitTime, function(x){max(sf::st_distance(x))})
maxNetDisp <- do.call("rbind",maxNetDispL)

rm(moveObjSplitTime);rm(maxNetDispL)

dmax7d <- 
  if(is.null(maxNetDisp)) NULL else {
  data.frame(keyName=row.names(maxNetDisp), dmax7d=maxNetDisp[,1], row.names=NULL) |> 
  mutate(week = str_split(keyName, "_", simplify = TRUE)[,3],
         year_week = paste(str_split(keyName, "_", simplify = TRUE)[,2],
                           str_split(keyName, "_", simplify = TRUE)[,3], sep ="_"),
         individual_id = str_split(keyName, "_", simplify = TRUE)[,1]) |> 
  filter(!is.na(dmax7d))|> 
  group_by(individual_id) |>  filter(n() >= 10) |> ungroup() |>  
  mutate(id.year_week = paste(individual_id,year_week,sep="_")) |> 
  left_join(mean.coord, by = "id.year_week") |> 
  dplyr::select(individual_id,week, year_week,dmax7d,mean.x, mean.y) }

dmax7d <- if(is.null(dmax7d)) NULL else 
          if(nrow(dmax7d) == 0) NULL else dmax7d 
  
if(is.null(dmax7d)) NULL else {
# Spatial annotation 10km
cell_info <- dgGEO_to_SEQNUM(dggs.10, mean.coord$mean.x, mean.coord$mean.y)
dmax7d$grid.id.10km <- cell_info$seqnum
centers <- dgSEQNUM_to_GEO(dggs.10, dmax7d$grid.id.10km)
dmax7d$lon.10km <- centers$lon_deg
dmax7d$lat.10km <- centers$lat_deg
rm(cell_info);rm(centers)

# Spatial annotation 1km
cell_info <- dgGEO_to_SEQNUM(dggs.1, mean.coord$mean.x, mean.coord$mean.y)
dmax7d$grid.id.1km <- cell_info$seqnum
centers <- dgSEQNUM_to_GEO(dggs.1, dmax7d$grid.id.1km)
dmax7d$lon.1km <- centers$lon_deg
dmax7d$lat.1km <- centers$lat_deg
rm(cell_info);rm(centers)
}

# # Visualization
# pts_sf <- locs24h.sf |> group_by(id.year_week) |> filter(id.year_week %in% "3809647332_2009_31")
# pts_sf2 <- mean.coord |> filter(id.year_week %in% "3809647332_2009_31") |> sf::st_as_sf(coords = c("mean.x", "mean.y"))
# tmp1 <- dmax7d |> filter(year_week %in% "2009_31") |> dplyr::select("grid.id.10km")
# tmp2 <- dmax7d |> filter(year_week %in% "2009_31") |> dplyr::select("grid.id.1km")
# intersect_cells_10 <- dgcellstogrid(dggs.10, tmp1[[1]])
# intersect_cells_1 <- dgcellstogrid(dggs.1, tmp2[[1]])
# mapview(intersect_cells_10, map.types = c("OpenStreetMap.DE"), col.regions = "red")+
#  mapview(intersect_cells_1, alpha = 1, col.regions = "darkred")+
#  mapview(pts_sf, col.regions = "black")+
#  mapview(pts_sf2, col.regions = "yellow")

rm(mean.coord);rm(locs24h.sf)

## ----Maximum annual displacement distance-------------------------------------------------------------
#' Based on daily (weekly) relocations we calculated the maximum annual displacement 
#' distance from all pairwise comparisons. We included only individuals with at least 36 weeks (9 months) of data

locs7d <- flatten(data_resampled[,"animlocs.weekly"]) %>% bind_rows() %>% tibble() %>%
  mutate(year = as.numeric(strftime(t_,format="%Y")),
         id.year = paste(id,year,sep="_")) |>   
  group_by(id.year)  |> filter(n() >= 36) |> ungroup()

mean.coord <- locs7d |> group_by(id.year) |> 
  mutate(mean.x = mean(x_, na.rm=T),
         mean.y = mean(y_, na.rm=T)) |> 
  dplyr::select(id.year, mean.x, mean.y) |> distinct()

locs7d.sf <- sf::st_as_sf(locs7d,
                           coords = c("x_", "y_"),
                           crs = 4326)

moveObjSplitTime <- split(locs7d.sf, locs7d.sf$id.year)
maxNetDispL <- lapply(moveObjSplitTime, function(x){max(sf::st_distance(x))})
maxNetDisp <- do.call("rbind",maxNetDispL)

rm(moveObjSplitTime);rm(maxNetDispL);rm(locs7d)

dmax12m <- 
  if(is.null(maxNetDisp)) NULL else {
  data.frame(keyName=row.names(maxNetDisp), dmax12m=maxNetDisp[,1], row.names=NULL) |> 
  filter(!is.na(dmax12m)) |> 
  mutate(year = str_split(keyName, "_", simplify = TRUE)[,2],
         individual_id = str_split(keyName, "_", simplify = TRUE)[,1]) |> 
  mutate(id.year = paste(individual_id,year,sep="_")) |>  left_join(mean.coord, by = "id.year") |> 
  dplyr::select(individual_id,year,dmax12m,mean.x, mean.y) }

dmax12m <- if(is.null(dmax12m)) NULL else 
  if(nrow(dmax12m) == 0) NULL else dmax12m 

if(is.null(dmax12m)) NULL else {
# Spatial annotation 10km
cell_info <- dgGEO_to_SEQNUM(dggs.10, mean.coord$mean.x, mean.coord$mean.y)
dmax12m$grid.id.10km <- cell_info$seqnum
centers <- dgSEQNUM_to_GEO(dggs.10, dmax12m$grid.id.10km)
dmax12m$lon.10km <- centers$lon_deg
dmax12m$lat.10km <- centers$lat_deg
rm(cell_info);rm(centers)

# Spatial annotation 1km
cell_info <- dgGEO_to_SEQNUM(dggs.1, mean.coord$mean.x, mean.coord$mean.y)
dmax12m$grid.id.1km <- cell_info$seqnum
centers <- dgSEQNUM_to_GEO(dggs.1, dmax12m$grid.id.1km)
dmax12m$lon.1km <- centers$lon_deg
dmax12m$lat.1km <- centers$lat_deg
rm(cell_info);rm(centers)
}

# # # Visualization
# pts_sf <- locs7d.sf |> filter(year %in% 2011)
#  pts_sf2 <- mean.coord |> filter(id.year %in% "3809647332_2011") |> sf::st_as_sf(coords = c("mean.x", "mean.y"))
#  tmp1 <- dmax12m |> filter(year %in% 2011) |> dplyr::select("grid.id.10km")
#  tmp2 <- dmax12m |> filter(year %in% 2011) |> dplyr::select("grid.id.1km")
#  intersect_cells_10 <- dgcellstogrid(dggs.10, tmp1[[1]])
#  intersect_cells_1 <- dgcellstogrid(dggs.1, tmp2[[1]])
# mapview(intersect_cells_10, map.types = c("OpenStreetMap.DE"), col.regions = "red")+
#  mapview(intersect_cells_1, alpha = 1, col.regions = "darkred")+
#  mapview(pts_sf, col.regions = "black")+
#  mapview(pts_sf2, col.regions = "yellow")

rm(mean.coord);rm(locs7d.sf)

## ----Range sizes - MCPs-------------------------------------------------------------
## ----Daily MCP-------------------------------------------------------------
#' Daily range used based on hourly relocations including only individuals with at least 12 locations on a given day.
dat.mcp.daily <- flatten(data_resampled[,"animlocs.1hourly"]) %>% bind_rows() %>% 
  tibble() %>% mutate(ymd = as.character(format(as.Date(t_), "%Y-%m-%d")))  %>%
  filter(!is.na(x_)) %>% filter(!is.na(y_)) %>% 
  mutate(id.day = paste(id,ymd,sep=".")) %>% group_by(id.day) %>% filter(n() > 12) %>% ungroup() %>% 
  dplyr::select(x_,y_,id.day) 

      mean.coord <- dat.mcp.daily |> group_by(id.day) |> 
        mutate(mean.x = mean(x_, na.rm=T),
               mean.y = mean(y_, na.rm=T)) |> 
        dplyr::select(id.day, mean.x, mean.y) |> distinct()

      tryCatch({
  coordinates(dat.mcp.daily) <- c("x_","y_")
  proj4string(dat.mcp.daily) <- CRS("EPSG:4326")
  
  # Bonne equal area projection - https://spatialreference.org/ref/esri/54024/
  # "+proj=bonne +lat_1=60 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"
  # Mollweide projection
  dat.mcp.daily <- spTransform(dat.mcp.daily,sp::CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"))
      }, error = function(e) {NA})
      
mcp.daily <- 
  if(nrow(dat.mcp.daily)==0) NULL else {
  mcp(dat.mcp.daily, percent = 95, unout = c( "m2")) %>% data.frame() %>% 
  left_join(mean.coord, by = c("id" = "id.day")) |> 
  mutate(ymd = str_split(id, '[.]', simplify = TRUE)[,2],
         individual_id = str_split(id, '[.]', simplify = TRUE)[,1]) %>% 
  #dplyr::select(individual_id,ymd,area, mean.x, mean.y)|> 
  filter(!is.na(area)) }

if(is.null(mcp.daily)) NULL else {
  # Spatial annotation 10km
  mean.coord$grid.id.10km <- dgGEO_to_SEQNUM(dggs.10, mean.coord$mean.x, mean.coord$mean.y)$seqnum
  mcp.daily <- mcp.daily |> left_join(mean.coord[,c("id.day","grid.id.10km")], by = c("id" = "id.day"))
  centers <- dgSEQNUM_to_GEO(dggs.10, mcp.daily$grid.id.10km)
  mcp.daily$lon.10km <- centers$lon_deg
  mcp.daily$lat.10km <- centers$lat_deg
  rm(centers)
  
  # Spatial annotation 1km
  mean.coord$grid.id.1km <- dgGEO_to_SEQNUM(dggs.1, mean.coord$mean.x, mean.coord$mean.y)$seqnum
  mcp.daily <- mcp.daily |> left_join(mean.coord[,c("id.day","grid.id.1km")], by = c("id" = "id.day"))
  centers <- dgSEQNUM_to_GEO(dggs.1, mcp.daily$grid.id.1km)
  mcp.daily$lon.1km <- centers$lon_deg
  mcp.daily$lat.1km <- centers$lat_deg
  rm(centers)
  
  mcp.daily <- mcp.daily |> 
  dplyr::select(individual_id,ymd,area, mean.x, mean.y,
                grid.id.10km, lon.10km,  lat.10km, 
                grid.id.1km,  lon.1km,   lat.1km) 
}

# # #Visualization
# pts_sf <- st_as_sf(dat.mcp.daily) |> st_transform(crs=4326) |>  filter(id.day %in% "3809647332.2009-08-01")
# pts_sf2 <- mean.coord |> filter( id.day %in% "3809647332.2009-08-01") |> sf::st_as_sf(coords = c("mean.x", "mean.y"))
# tmp1 <- mcp.daily |> filter(ymd %in% "2009-08-01") |> dplyr::select("grid.id.10km")
# tmp2 <- mcp.daily |> filter(ymd %in% "2009-08-01") |> dplyr::select("grid.id.1km")
# intersect_cells_10 <- dgcellstogrid(dggs.10, tmp1[[1]])
# intersect_cells_1 <- dgcellstogrid(dggs.1, tmp2[[1]])
# poly <- mcp(dat.mcp.daily, percent = 95, unout = c( "m2"))[4,]
# mapview(intersect_cells_10, map.types = c("OpenStreetMap.DE"), col.regions = "red")+
#  mapview(intersect_cells_1, alpha = 1, col.regions = "darkred")+
#  mapview(pts_sf2, col.regions = "yellow")+
#   mapview(poly, col.regions = "gray")+
#   mapview(pts_sf, cex= 0.8, col.regions = "black")

rm(mean.coord);rm(dat.mcp.daily)

## ----Weekly MCP-------------------------------------------------------------
dat.mcp.weekly <- flatten(data_resampled[,"animlocs.1hourly"]) %>% bind_rows() %>%  tibble() %>% 
 mutate(week = as.numeric(strftime(t_,format="%W")), 
        year = as.numeric(strftime(t_,format="%Y")),
        year_week = paste(year,week, sep="_")) %>% 
  mutate(id.week = paste(id,year_week,sep=".")) %>% 
  filter(!is.na(x_)) %>% filter(!is.na(y_)) %>%
  group_by(id.week) %>% filter(n() > 84) %>% ungroup() %>% 
  dplyr::select(x_,y_,id.week) 

      mean.coord <- dat.mcp.weekly |> group_by(id.week) |> 
        mutate(mean.x = mean(x_, na.rm=T),
               mean.y = mean(y_, na.rm=T)) |> 
        dplyr::select(id.week, mean.x, mean.y) |> distinct()

      tryCatch({
coordinates(dat.mcp.weekly) <- c("x_","y_")
proj4string(dat.mcp.weekly) <- CRS("EPSG:4326")
dat.mcp.weekly <- spTransform(dat.mcp.weekly,
                              sp::CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"))
      }, error = function(e) {NA})
      
mcp.weekly <- 
  if(nrow(dat.mcp.weekly)==0) NULL else {
  mcp(dat.mcp.weekly, percent = 95, unout = c("m2")) %>% data.frame() %>%
  left_join(mean.coord, by = c("id" = "id.week")) |> 
  mutate(week = as.numeric(stringr::str_extract(id, "(\\d+$)")),
         year_week = stringr::str_extract(id, "[^.]*$"),
         individual_id = str_extract(id, "[^.]+")) |> 
  filter(!is.na(area)) }

if(is.null(mcp.weekly)) NULL else {
  # Spatial annotation 10km
  mean.coord$grid.id.10km <- dgGEO_to_SEQNUM(dggs.10, mean.coord$mean.x, mean.coord$mean.y)$seqnum
  mcp.weekly <- mcp.weekly |> left_join(mean.coord[,c("id.week","grid.id.10km")], by = c("id" = "id.week"))
  centers <- dgSEQNUM_to_GEO(dggs.10, mcp.weekly$grid.id.10km)
  mcp.weekly$lon.10km <- centers$lon_deg
  mcp.weekly$lat.10km <- centers$lat_deg
  rm(centers)
  
  # Spatial annotation 1km
  mean.coord$grid.id.1km <- dgGEO_to_SEQNUM(dggs.1, mean.coord$mean.x, mean.coord$mean.y)$seqnum
  mcp.weekly <- mcp.weekly |> left_join(mean.coord[,c("id.week","grid.id.1km")], by = c("id" = "id.week"))
  centers <- dgSEQNUM_to_GEO(dggs.1, mcp.weekly$grid.id.1km)
  mcp.weekly$lon.1km <- centers$lon_deg
  mcp.weekly$lat.1km <- centers$lat_deg
  rm(centers)
  
  mcp.weekly <- mcp.weekly |> 
    dplyr::select(individual_id, week, year_week, area, mean.x, mean.y, 
                  grid.id.10km, lon.10km,  lat.10km, 
                  grid.id.1km,  lon.1km,   lat.1km)
}

#Visualization
# pts_sf <- st_as_sf(dat.mcp.weekly) |> st_transform(crs=4326) |>  filter(id.week %in% "3809647332.2009_31")
# pts_sf2 <- mean.coord |> filter(id.week %in% "3809647332.2009_31") |> sf::st_as_sf(coords = c("mean.x", "mean.y"))
# tmp1 <- mcp.weekly |> filter(year_week %in% "2009_31") |> dplyr::select("grid.id.10km")
# tmp2 <- mcp.weekly |> filter(year_week %in% "2009_31") |> dplyr::select("grid.id.1km")
# intersect_cells_10 <- dgcellstogrid(dggs.10, tmp1[[1]])
# intersect_cells_1 <- dgcellstogrid(dggs.1, tmp2[[1]])
# poly <- mcp(dat.mcp.weekly, percent = 95, unout = c("m2"))[2,]
# mapview(intersect_cells_10, map.types = c("OpenStreetMap.DE"), col.regions = "red")+
#   mapview(intersect_cells_1, alpha = 1, col.regions = "darkred")+
#   mapview(pts_sf2, col.regions = "yellow")+
#   mapview(poly, col.regions = "gray")+
#   mapview(pts_sf, cex= 0.2, col.regions = "black")

rm(dat.mcp.weekly);rm(mean.coord)

## ----Monthly MCP-------------------------------------------------------------
dat.mcp.monthly <- flatten(data_resampled[,"animlocs.daily"]) %>% bind_rows() %>% 
  tibble() %>% 
  mutate(month = as.numeric(strftime(t_,format="%m")), 
         year = as.numeric(strftime(t_,format="%Y"))) %>% 
  mutate(year_month = paste(year,month,sep="_")) |> 
  mutate(id.month = paste(id,year_month,sep=".")) %>% 
  filter(!is.na(x_)) %>% filter(!is.na(y_)) %>% group_by(id.month) %>% 
  filter(n() > 14) %>% ungroup() %>% dplyr::select(x_,y_,id.month) 

      mean.coord <- dat.mcp.monthly |> group_by(id.month) |> 
        mutate(mean.x = mean(x_, na.rm=T),
               mean.y = mean(y_, na.rm=T)) |> 
        dplyr::select(id.month, mean.x, mean.y) |> distinct()

      tryCatch({
    coordinates(dat.mcp.monthly) <- c("x_","y_")
    proj4string(dat.mcp.monthly) <- CRS("EPSG:4326")
    dat.mcp.monthly <- spTransform(dat.mcp.monthly,sp::CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"))
    }, error = function(e) {NA})

mcp.monthly <- 
  if(nrow(dat.mcp.monthly)==0) NULL else {
  mcp(dat.mcp.monthly, percent = 95, unout = c("m2")) %>% data.frame() %>%   
      left_join(mean.coord, by = c("id" = "id.month")) |> 
      mutate(month = as.numeric(stringr::str_extract(id, "(\\d+$)")),
         year_month = stringr::str_extract(id, "[^.]*$"),
         individual_id = str_extract(id, "[^.]+")) %>% 
  #dplyr::select(individual_id,month,year_month,area, mean.x, mean.y)|> 
  filter(!is.na(area)) }

if(is.null(mcp.monthly)) NULL else {
  # Spatial annotation 10km
  mean.coord$grid.id.10km <- dgGEO_to_SEQNUM(dggs.10, mean.coord$mean.x, mean.coord$mean.y)$seqnum
  mcp.monthly <- mcp.monthly |> left_join(mean.coord[,c("id.month","grid.id.10km")], by = c("id" = "id.month"))
  centers <- dgSEQNUM_to_GEO(dggs.10, mcp.monthly$grid.id.10km)
  mcp.monthly$lon.10km <- centers$lon_deg
  mcp.monthly$lat.10km <- centers$lat_deg
  rm(centers)
  
  # Spatial annotation 1km
  mean.coord$grid.id.1km <- dgGEO_to_SEQNUM(dggs.1, mean.coord$mean.x, mean.coord$mean.y)$seqnum
  mcp.monthly <- mcp.monthly |> left_join(mean.coord[,c("id.month","grid.id.1km")], by = c("id" = "id.month"))
  centers <- dgSEQNUM_to_GEO(dggs.1, mcp.monthly$grid.id.1km)
  mcp.monthly$lon.1km <- centers$lon_deg
  mcp.monthly$lat.1km <- centers$lat_deg
  rm(centers)
  
  mcp.monthly <- mcp.monthly |> 
    dplyr::select(individual_id,month,year_month,area, mean.x, mean.y, 
                  grid.id.10km, lon.10km,  lat.10km, 
                  grid.id.1km,  lon.1km,   lat.1km)
}

# #Visualization
# pts_sf <- st_as_sf(dat.mcp.monthly) |> st_transform(crs=4326) |>  filter(id.month %in% "3809647332.2009_9")
# pts_sf2 <- mean.coord |> filter(id.month %in% "3809647332.2009_9") |> sf::st_as_sf(coords = c("mean.x", "mean.y"))
# tmp1 <- mcp.monthly |> filter(year_month %in% "2009_9") |> dplyr::select("grid.id.10km")
# tmp2 <- mcp.monthly |> filter(year_month %in% "2009_9") |> dplyr::select("grid.id.1km")
# intersect_cells_10 <- dgcellstogrid(dggs.10, tmp1[[1]])
# intersect_cells_1 <- dgcellstogrid(dggs.1, tmp2[[1]])
# poly <- mcp(dat.mcp.monthly, percent = 95, unout = c("m2"))[5,]
# mapview(intersect_cells_10, map.types = c("OpenStreetMap.DE"), col.regions = "red")+
#   mapview(intersect_cells_1, alpha = 1, col.regions = "darkred")+
#   mapview(pts_sf2, col.regions = "yellow")+
#   mapview(poly, col.regions = "gray")+
#   mapview(pts_sf, cex= 0.2, col.regions = "black")

rm(dat.mcp.monthly);rm(mean.coord)

## ----Annual MCP-------------------------------------------------------------
dat.mcp.annual <- flatten(data_resampled[,"animlocs.weekly"]) %>% bind_rows() %>% 
  tibble() |>   mutate(year = as.numeric(strftime(t_,format="%Y"))) %>% 
  mutate(id.year = paste(id,year,sep=".")) %>% 
  filter(!is.na(x_)) %>% filter(!is.na(y_)) %>% group_by(id.year) %>% 
  filter(n() > 36) %>% ungroup() %>% dplyr::select(x_,y_,id.year) 

      mean.coord <- dat.mcp.annual |> group_by(id.year) |> 
        mutate(mean.x = mean(x_, na.rm=T),
               mean.y = mean(y_, na.rm=T)) |> 
        dplyr::select(id.year, mean.x, mean.y) |> distinct()
      
      tryCatch({
    coordinates(dat.mcp.annual) <- c("x_","y_")
    proj4string(dat.mcp.annual) <- CRS("EPSG:4326")
    dat.mcp.annual <- spTransform(dat.mcp.annual,sp::CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"))
    }, error = function(e) {NA})

mcp.annual <- 
  if(nrow(dat.mcp.annual)==0) NULL else {
  mcp(dat.mcp.annual, percent = 95, unout = c("m2")) %>% data.frame() %>% 
  left_join(mean.coord, by = c("id" = "id.year")) |> 
  mutate(id.year = id,
         year = stringr::str_extract(id, "[^.]*$"),
         individual_id = str_extract(id, "[^.]+")) %>% 
  #dplyr::select(individual_id, id.year, year,area, mean.x, mean.y)|> 
  filter(!is.na(area)) 
  }

if(is.null(mcp.annual)) NULL else {
  # Spatial annotation 10km
  mean.coord$grid.id.10km <- dgGEO_to_SEQNUM(dggs.10, mean.coord$mean.x, mean.coord$mean.y)$seqnum
  mcp.annual <- mcp.annual |> left_join(mean.coord[,c("id.year","grid.id.10km")], by = c("id" = "id.year"))
  centers <- dgSEQNUM_to_GEO(dggs.10, mcp.annual$grid.id.10km)
  mcp.annual$lon.10km <- centers$lon_deg
  mcp.annual$lat.10km <- centers$lat_deg
  rm(centers)
  
  # Spatial annotation 1km
  mean.coord$grid.id.1km <- dgGEO_to_SEQNUM(dggs.1, mean.coord$mean.x, mean.coord$mean.y)$seqnum
  mcp.annual <- mcp.annual |> left_join(mean.coord[,c("id.year","grid.id.1km")], by = c("id" = "id.year"))
  centers <- dgSEQNUM_to_GEO(dggs.1, mcp.annual$grid.id.1km)
  mcp.annual$lon.1km <- centers$lon_deg
  mcp.annual$lat.1km <- centers$lat_deg
  rm(centers)
  
  mcp.annual <- mcp.annual |> 
    dplyr::select(individual_id,id.year, year, area, mean.x, mean.y, 
                  grid.id.10km, lon.10km,  lat.10km, 
                  grid.id.1km,  lon.1km,   lat.1km)
}


# # #Visualization
# pts_sf <- st_as_sf(dat.mcp.annual) |> st_transform(crs=4326) |>  filter(id.year %in% "3809647332.2010")
# pts_sf2 <- mean.coord |> filter(id.year %in% "3809647332.2010") |> sf::st_as_sf(coords = c("mean.x", "mean.y"))
# tmp1 <- mcp.annual |> filter(year %in% "2010") |> dplyr::select("grid.id.10km")
# tmp2 <- mcp.annual |> filter(year %in% "2010") |> dplyr::select("grid.id.1km")
# intersect_cells_10 <- dgcellstogrid(dggs.10, tmp1[[1]])
# intersect_cells_1 <- dgcellstogrid(dggs.1, tmp2[[1]])
# poly <- mcp(dat.mcp.annual, percent = 95, unout = c("m2"))[1,]
# mapview(intersect_cells_10, map.types = c("OpenStreetMap.DE"), col.regions = "red")+
#   mapview(intersect_cells_1, alpha = 1, col.regions = "darkred")+
#   mapview(pts_sf2, col.regions = "yellow")+
#   mapview(poly, col.regions = "gray")+
#   mapview(pts_sf, cex= 0.2, col.regions = "black")

rm(dat.mcp.annual);rm(mean.coord)

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

## ----Summarize movement metrics per individual-------------------------------------------------------------

# ----Function d1h-------------------------------------------------------------
FDispl1h<-function(x)
  {
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n1h = NA,
                       d1h.mean = NA,d1h.median = NA,d1h.cv = NA,d1h.95 = NA,d1h.05 = NA)
  } else {
  
  individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
  
  # Get sample size per individual
  n1h<-as.numeric(with(x, tapply(as.character(x$t_),individual_id, length)))
  
   # 1hr Displacement
  d1h.mean <- as.numeric(with(x, tapply(x$d1h+0.001,individual_id, mean, na.rm=T)))
  d1h.median <- as.numeric(with(x, tapply(x$d1h+0.001,individual_id, median, na.rm=T)))
  d1h.cv <- as.numeric(with(x, tapply(x$d1h+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
  d1h.95 <- as.numeric(with(x, tapply(x$d1h+0.001,individual_id, quantile,.95, na.rm=T)))
  d1h.05 <- as.numeric(with(x, tapply(x$d1h+0.001,individual_id, quantile,.05, na.rm=T)))
  
  # build dataframe
  dats<-data.frame(individual_id,n1h,
                   d1h.mean,d1h.median,d1h.cv,d1h.95,d1h.05)
  
  return(dats)
}
}

Displ1h <- FDispl1h(animlocs.1hourly_sl[,c("individual_id","t_","d1h")])

## ----function Max24h Displacements------
FMaxDispl24h<-function(x)
  {
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n.dmax24h.days = NA,
                       dmax24h.mean = NA,dmax24h.median = NA,dmax24h.cv = NA,dmax24h.95 = NA,dmax24h.05 = NA)
  } else {
    
  individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
  
  # Get sample size per indivindividual_idual
  n.dmax24h.days<-as.numeric(with(x, tapply(x$ymd,individual_id, length)))
  
  # 24hr Displacement
  dmax24h.mean<-as.numeric(with(x, tapply(x$dmax24h+0.001,individual_id, mean, na.rm=T)))
  dmax24h.median<-as.numeric(with(x, tapply(x$dmax24h+0.001,individual_id, median, na.rm=T)))
  dmax24h.cv<-as.numeric(with(x, tapply(x$dmax24h+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
  dmax24h.95<-as.numeric(with(x, tapply(x$dmax24h+0.001,individual_id, quantile,.95, na.rm=T)))
  dmax24h.05<-as.numeric(with(x, tapply(x$dmax24h+0.001,individual_id, quantile,.05, na.rm=T)))
  
  # build dataframe
  dats<-data.frame(individual_id,n.dmax24h.days,
                   dmax24h.mean,dmax24h.median,dmax24h.cv,dmax24h.95,dmax24h.05)
  
  return(dats)
}
}

MaxDispl24h <- FMaxDispl24h(dmax24)

## ----function to summarize 24h Displacements---------
FDispl24h<-function(x)
  {
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n24h.days = NA,
                       d24h.mean = NA,d24h.median = NA,d24h.cv = NA,d24h.95 = NA,d24h.05 = NA)
  } else {
  
  individual_id<-with(x, tapply(as.character(x$individual_id),individual_id, unique))
  
  # Get sample size per indivindividual_idual
  n24h.days<-as.numeric(with(x, tapply(as.character(x$t_),individual_id, length)))
  
  # 24hr Displacement
  d24h.mean <- as.numeric(with(x, tapply(x$d24h+0.001,individual_id, mean, na.rm=T)))
  d24h.median <- as.numeric(with(x, tapply(x$d24h+0.001,individual_id, median, na.rm=T)))
  d24h.cv <- as.numeric(with(x, tapply(x$d24h+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
  d24h.95 <- as.numeric(with(x, tapply(x$d24h+0.001,individual_id, quantile,.95, na.rm=T)))
  d24h.05 <- as.numeric(with(x, tapply(x$d24h+0.001,individual_id, quantile,.05, na.rm=T)))
  
  # build dataframe
  dats<-data.frame(individual_id,n24h.days,
                   d24h.mean,d24h.median,d24h.cv,d24h.95,d24h.05)
  
  return(dats)
}
}

Displ24h <- FDispl24h(animlocs.daily_sl[,c("individual_id","t_","d24h")])

## ----function to summarize Max7d Displacements-------
FMaxDispl7d<-function(x)
  {
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n.dmax7d.weeks = NA,
                       dmax7d.mean = NA,dmax7d.median = NA,dmax7d.cv = NA,dmax7d.95 = NA,dmax7d.05 = NA)
  } else {
    
  individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
  
  # Get sample size per indivindividual_idual
  n.dmax7d.weeks<-as.numeric(with(x, tapply(x$year_week,individual_id, length)))
  
  # 24hr Displacement
  dmax7d.mean<-as.numeric(with(x, tapply(x$dmax7d+0.001,individual_id, mean, na.rm=T)))
  dmax7d.median<-as.numeric(with(x, tapply(x$dmax7d+0.001,individual_id, median, na.rm=T)))
  dmax7d.cv<-as.numeric(with(x, tapply(x$dmax7d+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
  dmax7d.95<-as.numeric(with(x, tapply(x$dmax7d+0.001,individual_id, quantile,.95, na.rm=T)))
  dmax7d.05<-as.numeric(with(x, tapply(x$dmax7d+0.001,individual_id, quantile,.05, na.rm=T)))
  
  # build dataframe
  dats<-data.frame(individual_id,n.dmax7d.weeks,
                   dmax7d.mean,dmax7d.median,dmax7d.cv,dmax7d.95,dmax7d.05)
  
  return(dats)
}
}

MaxDispl7d <- FMaxDispl7d(dmax7d)

## ----function to summarize Max12m Displacements------
FMaxDispl12m<-function(x)
  {
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n.max12m.years = NA,
                       dmax12m.mean = NA,dmax12m.median = NA,dmax12m.cv = NA,dmax12m.95 = NA,dmax12m.05 = NA)
  } else {
    
  individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
  
  # Get sample size per indivindividual_idual
  n.max12m.years<-as.numeric(with(x, tapply(x$year,individual_id, length)))
  
  # 24hr Displacement
  dmax12m.mean<-as.numeric(with(x, tapply(x$dmax12m+0.001,individual_id, mean, na.rm=T)))
  dmax12m.median<-as.numeric(with(x, tapply(x$dmax12m+0.001,individual_id, median, na.rm=T)))
  dmax12m.cv<-as.numeric(with(x, tapply(x$dmax12m+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
  dmax12m.95<-as.numeric(with(x, tapply(x$dmax12m+0.001,individual_id, quantile,.95, na.rm=T)))
  dmax12m.05<-as.numeric(with(x, tapply(x$dmax12m+0.001,individual_id, quantile,.05, na.rm=T)))
  
  # build dataframe
  dats<-data.frame(individual_id,n.max12m.years,
                   dmax12m.mean,dmax12m.median,dmax12m.cv,dmax12m.95,dmax12m.05)
  
  return(dats)
}
}

MaxDispl12m <- FMaxDispl12m(dmax12m)

## ----function to summarize 1d MCP--------------------
FSumMCP1d<-function(x)
  {
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n.mcp24h.days = NA,
                       mcp24h.mean = NA,mcp24h.median = NA,mcp24h.cv = NA,mcp24h.95 = NA,mcp24h.05 = NA)
  } else {
    
  individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
  
  # Get sample size per indivindividual_idual
  n.mcp24h.days<-as.numeric(with(x, tapply(x$individual_id,individual_id, length)))
  
  # 24MCP Displacement
  mcp24h.mean<-as.numeric(with(x, tapply(x$area+0.001,individual_id, mean, na.rm=T)))
  mcp24h.median<-as.numeric(with(x, tapply(x$area+0.001,individual_id, median, na.rm=T)))
  mcp24h.cv<-as.numeric(with(x, tapply(x$area+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
  mcp24h.95<-as.numeric(with(x, tapply(x$area+0.001,individual_id, quantile,.95, na.rm=T)))
  mcp24h.05<-as.numeric(with(x, tapply(x$area+0.001,individual_id, quantile,.05, na.rm=T)))
  
  # build dataframe
  dats<-data.frame(individual_id,n.mcp24h.days,
                   mcp24h.mean,mcp24h.median,mcp24h.cv,mcp24h.95,mcp24h.05)
  
  return(dats)
}
}
MCP24h <- FSumMCP1d(mcp.daily)

## ----function to summarize 7d MCP--------------------
FSumMCP7d<-function(x)
  {
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n.mcp7d.weeks = NA,
                       mcp7d.mean = NA,mcp7d.median = NA,mcp7d.cv = NA,mcp7d.95 = NA,mcp7d.05 = NA)
  } else {
    
  individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
  
  # Get sample size per indivindividual_idual
  n.mcp7d.weeks<-as.numeric(with(x, tapply(x$year_week,individual_id, length)))
  
  # 24mcp Displacement
  mcp7d.mean<-as.numeric(with(x, tapply(x$area+0.001,individual_id, mean, na.rm=T)))
  mcp7d.median<-as.numeric(with(x, tapply(x$area+0.001,individual_id, median, na.rm=T)))
  mcp7d.cv<-as.numeric(with(x, tapply(x$area+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
  mcp7d.95<-as.numeric(with(x, tapply(x$area+0.001,individual_id, quantile,.95, na.rm=T)))
  mcp7d.05<-as.numeric(with(x, tapply(x$area+0.001,individual_id, quantile,.05, na.rm=T)))
  
  # build dataframe
  dats<-data.frame(individual_id,n.mcp7d.weeks,
                   mcp7d.mean,mcp7d.median,mcp7d.cv,mcp7d.95,mcp7d.05)
  
  return(dats)
}
}

MCP7d <- FSumMCP7d(mcp.weekly)

## ----function to summarize 1m MCP--------------------
FSumMCP1m<-function(x)
  {
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n.mcp1m.months = NA,
                       mcp1m.mean = NA,mcp1m.median = NA,mcp1m.cv = NA,mcp1m.95 = NA,mcp1m.05 = NA)
  } else {
    
  individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))

  n.mcp1m.months<-as.numeric(with(x, tapply(x$year_month,individual_id, length)))
  
  mcp1m.mean<-as.numeric(with(x, tapply(x$area+0.001,individual_id, mean, na.rm=T)))
  mcp1m.median<-as.numeric(with(x, tapply(x$area+0.001,individual_id, median, na.rm=T)))
  mcp1m.cv<-as.numeric(with(x, tapply(x$area+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
  mcp1m.95<-as.numeric(with(x, tapply(x$area+0.001,individual_id, quantile,.95, na.rm=T)))
  mcp1m.05<-as.numeric(with(x, tapply(x$area+0.001,individual_id, quantile,.05, na.rm=T)))
  
  dats<-data.frame(individual_id,n.mcp1m.months,
                   mcp1m.mean,mcp1m.median,mcp1m.cv,mcp1m.95,mcp1m.05)
  
  return(dats)
}
}

MCP1m <- FSumMCP1m(mcp.monthly)

## ----function to summarize 12m MCP--------------------
FSumMCP12m<-function(x)
  {
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n.mcp12m.years = NA,
                       mcp12m.mean = NA,mcp12m.median = NA,mcp12m.cv = NA,mcp12m.95 = NA,mcp12m.05 = NA)
  } else {
    
  individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))

  n.mcp12m.years<-as.numeric(with(x, tapply(x$year,individual_id, length)))
  
  mcp12m.mean<-as.numeric(with(x, tapply(x$area+0.001,individual_id, mean, na.rm=T)))
  mcp12m.median<-as.numeric(with(x, tapply(x$area+0.001,individual_id, median, na.rm=T)))
  mcp12m.cv<-as.numeric(with(x, tapply(x$area+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
  mcp12m.95<-as.numeric(with(x, tapply(x$area+0.001,individual_id, quantile,.95, na.rm=T)))
  mcp12m.05<-as.numeric(with(x, tapply(x$area+0.001,individual_id, quantile,.05, na.rm=T)))
  
  dats<-data.frame(individual_id,n.mcp12m.years,
                   mcp12m.mean,mcp12m.median,mcp12m.cv,mcp12m.95,mcp12m.05)
  
  return(dats)
  }}

MCP12m <- FSumMCP12m(mcp.annual)

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
#' ## PROBLEM - this does not work when some metrics could not be computed
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

movedata2$grid.id.100km <- paste(grid.id.100km, collapse = ";")
movedata2$grid.id.10km <- paste(grid.id.10km, collapse = ";")

MoveTrait.v0.1 <- 
  if (is.null(MoveTrait.v0.1) | nrow(MoveTrait.v0.1) == 0) NULL else {
    movedata2 %>% 
    left_join(MoveTrait.v0.1, by = "individual_id") %>% 
    droplevels() }

## ----Save database with summaries--------------------
saveRDS(MoveTrait.v0.1, file=paste0(pthtraitsum,"3809257699_3809647332.rds"))

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

