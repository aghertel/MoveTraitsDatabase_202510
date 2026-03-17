## ----Monthly MCP-------------------------------------------------------------
calc_mcp1m <- function(trk, 
                       dggs_10, 
                       dggs_1, 
                       min_days_n = 14) 
{
  #trk <- animlocs.daily
  dat.mcp.monthly <- trk %>% 
  tibble() %>% 
  mutate(month = as.numeric(strftime(t_,format="%m")), 
         year = as.numeric(strftime(t_,format="%Y"))) %>% 
  mutate(year_month = paste(year,month,sep="_")) |> 
  mutate(id.month = paste(individual_id,year_month,sep=".")) %>% 
  filter(!is.na(x_)) %>% filter(!is.na(y_)) %>% group_by(id.month) %>% 
  filter(n() > min_days_n) %>% ungroup() %>% dplyr::select(x_,y_,id.month) 

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
      filter(!is.na(area)) }

if(is.null(mcp.monthly)) NULL else {
  # Spatial annotation 10km
  mean.coord$grid.id.10km <- dgGEO_to_SEQNUM(dggs.10, mean.coord$mean.x, mean.coord$mean.y)$seqnum
  mcp.monthly <- mcp.monthly |> left_join(mean.coord[,c("id.month","grid.id.10km")], by = c("id" = "id.month"))
  centers_10 <- dgSEQNUM_to_GEO(dggs.10, mcp.monthly$grid.id.10km)
  mcp.monthly$lon.10km <- centers_10$lon_deg
  mcp.monthly$lat.10km <- centers_10$lat_deg
  
  # Spatial annotation 1km
  mean.coord$grid.id.1km <- dgGEO_to_SEQNUM(dggs.1, mean.coord$mean.x, mean.coord$mean.y)$seqnum
  mcp.monthly <- mcp.monthly |> left_join(mean.coord[,c("id.month","grid.id.1km")], by = c("id" = "id.month"))
  centers_1 <- dgSEQNUM_to_GEO(dggs.1, mcp.monthly$grid.id.1km)
  mcp.monthly$lon.1km <- centers_1$lon_deg
  mcp.monthly$lat.1km <- centers_1$lat_deg
  
  mcp.monthly <- mcp.monthly |> 
    dplyr::select(individual_id,month,year_month,area, mean.x, mean.y, 
                  grid.id.10km, lon.10km,  lat.10km, 
                  grid.id.1km,  lon.1km,   lat.1km)
}
rm(dat.mcp.monthly);rm(mean.coord)
return(mcp.monthly);rm(centers_10);rm(centers_1)
}



## ----function to summarize 1m MCP--------------------
f_sum.ind.mcp1m<-function(x)
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
