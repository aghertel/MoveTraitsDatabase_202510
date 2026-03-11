## ----Annual MCP-------------------------------------------------------------
calc_mcp12m <- function(trk, 
                       dggs_10, 
                       dggs_1, 
                       min_weeks_n = 36) 
{
  dat.mcp.annual <- trk %>% 
  tibble() |>   mutate(year = as.numeric(strftime(t_,format="%Y"))) %>% 
  mutate(id.year = paste(individual_id,year,sep=".")) %>% 
  filter(!is.na(x_)) %>% filter(!is.na(y_)) %>% group_by(id.year) %>% 
  filter(n() > min_weeks_n) %>% ungroup() %>% dplyr::select(x_,y_,id.year) 

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
  centers_10 <- dgSEQNUM_to_GEO(dggs.10, mcp.annual$grid.id.10km)
  mcp.annual$lon.10km <- centers_10$lon_deg
  mcp.annual$lat.10km <- centers_10$lat_deg
  
  # Spatial annotation 1km
  mean.coord$grid.id.1km <- dgGEO_to_SEQNUM(dggs.1, mean.coord$mean.x, mean.coord$mean.y)$seqnum
  mcp.annual <- mcp.annual |> left_join(mean.coord[,c("id.year","grid.id.1km")], by = c("id" = "id.year"))
  centers_1 <- dgSEQNUM_to_GEO(dggs.1, mcp.annual$grid.id.1km)
  mcp.annual$lon.1km <- centers_1$lon_deg
  mcp.annual$lat.1km <- centers_1$lat_deg

  
  mcp.annual <- mcp.annual |> 
    dplyr::select(individual_id,id.year, year, area, mean.x, mean.y, 
                  grid.id.10km, lon.10km,  lat.10km, 
                  grid.id.1km,  lon.1km,   lat.1km)
}
return(mcp.annual)
rm(dat.mcp.annual);rm(mean.coord);rm(centers_1); rm(centers_10)
}


## ----function to summarize 12m MCP--------------------
f_sum.ind.mcp12m<-function(x)
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

