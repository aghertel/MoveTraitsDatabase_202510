## ----Daily MCP-------------------------------------------------------------

#' Daily range used based on hourly relocations including only individuals with at least 12 locations on a given day.


calc_mcp24h <- function(trk, 
                         dggs_10, 
                         dggs_1, 
                         min_hours_n = 12) 
{
  dat.mcp.daily <- trk %>% 
  tibble() %>% mutate(ymd = as.character(format(as.Date(t_), "%Y-%m-%d")))  %>%
  filter(!is.na(x_)) %>% 
  filter(!is.na(y_)) %>% 
  mutate(id.day = paste(individual_id,ymd,sep=".")) %>% 
  group_by(id.day) %>% 
  filter(n() > min_hours_n) %>% 
  ungroup() %>% 
  dplyr::select(x_,y_,id.day) 

mean.coord <- dat.mcp.daily |> group_by(id.day) |> 
  mutate(mean.x = mean(x_, na.rm=T),
         mean.y = mean(y_, na.rm=T)) |> 
  dplyr::select(id.day, mean.x, mean.y) |> distinct()

tryCatch({
  coordinates(dat.mcp.daily) <- c("x_","y_")
  proj4string(dat.mcp.daily) <- CRS("EPSG:4326")
  dat.mcp.daily <- spTransform(dat.mcp.daily,sp::CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"))
}, error = function(e) {NA})

mcp.daily <- 
  if(nrow(dat.mcp.daily)==0) NULL else {
    mcp(dat.mcp.daily, percent = 95, unout = c( "m2")) %>% data.frame() %>% 
      left_join(mean.coord, by = c("id" = "id.day")) |> 
      mutate(ymd = str_split(id, '[.]', simplify = TRUE)[,2],
             individual_id = str_split(id, '[.]', simplify = TRUE)[,1]) %>% 
      filter(!is.na(area)) }

if(is.null(mcp.daily)) NULL else {
  # Spatial annotation 10km
  mean.coord$grid.id.10km <- dgGEO_to_SEQNUM(dggs.10, mean.coord$mean.x, mean.coord$mean.y)$seqnum
  mcp.daily <- mcp.daily |> left_join(mean.coord[,c("id.day","grid.id.10km")], by = c("id" = "id.day"))
  centers_10 <- dgSEQNUM_to_GEO(dggs.10, mcp.daily$grid.id.10km)
  mcp.daily$lon.10km <- centers_10$lon_deg
  mcp.daily$lat.10km <- centers_10$lat_deg
  
  # Spatial annotation 1km
  mean.coord$grid.id.1km <- dgGEO_to_SEQNUM(dggs.1, mean.coord$mean.x, mean.coord$mean.y)$seqnum
  mcp.daily <- mcp.daily |> left_join(mean.coord[,c("id.day","grid.id.1km")], by = c("id" = "id.day"))
  centers_1 <- dgSEQNUM_to_GEO(dggs.1, mcp.daily$grid.id.1km)
  mcp.daily$lon.1km <- centers_1$lon_deg
  mcp.daily$lat.1km <- centers_1$lat_deg
  
  mcp.daily <- mcp.daily |> 
    dplyr::select(individual_id,ymd,area, mean.x, mean.y,
                  grid.id.10km, lon.10km,  lat.10km, 
                  grid.id.1km,  lon.1km,   lat.1km) 
}

return(mcp.daily)
rm(mean.coord);rm(dat.mcp.daily);rm(centers_10);rm(centers_1)
}

## ----function to summarize 1d MCP--------------------
f_sum.ind.mcp24h<-function(x)
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


## ----summarize mcp24h at monthly individual level---------
f_sum.monthly.ind.mcp24h<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    return(data.frame(individual_id = NA,n.mcp24h.days = NA,
                      mcp24h.mean = NA,mcp24h.median = NA,mcp24h.cv = NA,mcp24h.95 = NA,mcp24h.05 = NA))
  } 
  
  # derive month and year from t_
  ym <- substr(x$ymd, 1, 7)  
  
  # group index: individual x month
  id_ym <- interaction(x$individual_id, ym, drop = TRUE)
  
  individual_id <- tapply(as.character(x$individual_id), id_ym, unique)
  ym_grp        <- tapply(ym, id_ym, unique)
  
  # sample size per individual-month
  n.mcp24h.days<-as.numeric(tapply(x$ymd, id_ym, length))
  
  # 24hr Displacement
  mcp24h.mean <- as.numeric(with(x, tapply(x$area+0.001,id_ym, mean, na.rm=T)))
  mcp24h.median <- as.numeric(with(x, tapply(x$area+0.001,id_ym, median, na.rm=T)))
  mcp24h.cv <- as.numeric(with(x, tapply(x$area+0.001,id_ym, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
  mcp24h.95 <- as.numeric(with(x, tapply(x$area+0.001,id_ym, quantile,.95, na.rm=T)))
  mcp24h.05 <- as.numeric(with(x, tapply(x$area+0.001,id_ym, quantile,.05, na.rm=T)))
  
  # build dataframe
  dats<-data.frame(individual_id,year,month,n.mcp24h.days,
                   mcp24h.mean,mcp24h.median,mcp24h.cv,mcp24h.95,mcp24h.05)
  
}


