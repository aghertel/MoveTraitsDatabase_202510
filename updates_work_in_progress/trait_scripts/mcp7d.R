## ----Weekly MCP-------------------------------------------------------------

calc_mcp7d <- function(trk, 
                        dggs_10, 
                        dggs_1, 
                        min_hours_n = 84) 
{
dat.mcp.weekly <- trk %>% 
  mutate(week = as.numeric(strftime(t_,format="%W")), 
         year = as.numeric(strftime(t_,format="%Y")),
         year_week = paste(year,week, sep="_")) %>% 
  mutate(id.week = paste(individual_id,year_week,sep=".")) %>% 
  filter(!is.na(x_)) %>% filter(!is.na(y_)) %>%
  group_by(id.week) %>% filter(n() > min_hours_n) %>% ungroup() %>% 
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
  centers_10 <- dgSEQNUM_to_GEO(dggs.10, mcp.weekly$grid.id.10km)
  mcp.weekly$lon.10km <- centers_10$lon_deg
  mcp.weekly$lat.10km <- centers_10$lat_deg

  
  # Spatial annotation 1km
  mean.coord$grid.id.1km <- dgGEO_to_SEQNUM(dggs.1, mean.coord$mean.x, mean.coord$mean.y)$seqnum
  mcp.weekly <- mcp.weekly |> left_join(mean.coord[,c("id.week","grid.id.1km")], by = c("id" = "id.week"))
  centers_1 <- dgSEQNUM_to_GEO(dggs.1, mcp.weekly$grid.id.1km)
  mcp.weekly$lon.1km <- centers_1$lon_deg
  mcp.weekly$lat.1km <- centers_1$lat_deg
  
  mcp.weekly <- mcp.weekly |> 
    dplyr::select(individual_id, week, year_week, area, mean.x, mean.y, 
                  grid.id.10km, lon.10km,  lat.10km, 
                  grid.id.1km,  lon.1km,   lat.1km)
}
return(mcp.weekly)
rm(dat.mcp.weekly);rm(mean.coord);rm(centers_10);rm(centers_1)
}

## ----function to summarize 7d MCP--------------------
f_sum.ind.mcp7d<-function(x)
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


## ----summarize mcp7d at monthly individual level---------
f_sum.monthly.ind.mcp7d<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,year=NA, month=NA,n.mcp7d.weeks = NA,
                       mcp7d.mean = NA,mcp7d.median = NA,mcp7d.cv = NA,mcp7d.95 = NA,mcp7d.05 = NA)
  } else {
    
    # derive month and year from t_
    year <- as.numeric(substr(x$year_week, 1, 4))  
    
    string_iso <- paste(y, sprintf("W%02d", as.numeric(x$week)), 1, sep="-")
    date <- ISOweek2date(string_iso)
    month <- format(date, "%m")
    ym <- paste(year,month,sep="-")
    
    # group index: individual x month
    id_ym <- interaction(x$individual_id, ym, drop = TRUE)
    
    individual_id <- tapply(as.character(x$individual_id), id_ym, unique)
    ym_grp        <- tapply(ym, id_ym, unique)
    
    # sample size per individual-month
    n.mcp7d.weeks <- as.numeric(tapply(x$year_week, id_ym, length))
    
    # 24hr Displacement
    mcp7d.mean<-as.numeric(with(x, tapply(x$area+0.001,id_ym, mean, na.rm=T)))
    mcp7d.median<-as.numeric(with(x, tapply(x$area+0.001,id_ym, median, na.rm=T)))
    mcp7d.cv<-as.numeric(with(x, tapply(x$area+0.001,id_ym, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
    mcp7d.95<-as.numeric(with(x, tapply(x$area+0.001,id_ym, quantile,.95, na.rm=T)))
    mcp7d.05<-as.numeric(with(x, tapply(x$area+0.001,id_ym, quantile,.05, na.rm=T)))
    
    year  <- as.numeric(sub(".*\\.(\\d{4})-\\d{2}$", "\\1", unique(id_ym)))
    month <- as.numeric(sub(".*\\.\\d{4}-(\\d{2})$", "\\1", unique(id_ym)))
    
    # build dataframe
    dats<-data.frame(individual_id,year,month,n.mcp7d.weeks,
                     mcp7d.mean,mcp7d.median,mcp7d.cv,mcp7d.95,mcp7d.05)
    
    return(dats)
  }
}
