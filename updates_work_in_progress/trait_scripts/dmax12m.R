
## ----Maximum annual displacement distance-------------------------------------------------------------
#' Based on daily (weekly) relocations we calculated the maximum annual displacement 
#' distance from all pairwise comparisons. We included only individuals with at least 36 weeks (9 months) of data

calc_dmax12m <- function(trk, 
                         dggs_10, 
                         dggs_1, 
                         min_weeks_n = 36) 
{
  
locs7d <- trk %>%
  mutate(year = as.numeric(strftime(t_,format="%Y")),
         id.year = paste(individual_id,year,sep="_")) |>   
  group_by(id.year)  |> filter(n() >= min_weeks_n) |> ungroup()

if (nrow(locs7d) == 0) return(NULL)

# Calculate mean coordinates per day
mean.coord <- locs7d |> group_by(id.year) |> 
  mutate(mean.x = mean(x_, na.rm=T),
         mean.y = mean(y_, na.rm=T)) |> 
  dplyr::select(id.year, mean.x, mean.y) |> distinct()

# Convert to sf and calculate max pairwise distance per day
locs7d.sf <- sf::st_as_sf(locs7d,
                          coords = c("x_", "y_"),
                          crs = 4326)

moveObjSplitTime <- split(locs7d.sf, locs7d.sf$id.year)
maxNetDispL <- lapply(moveObjSplitTime, function(x){max(sf::st_distance(x))})
maxNetDisp <- do.call("rbind",maxNetDispL)

# Process results
dmax12m <- 
  if(is.null(maxNetDisp)) NULL else {
    data.frame(keyName=row.names(maxNetDisp), dmax12m=maxNetDisp[,1], row.names=NULL) |> 
      filter(!is.na(dmax12m)) |> 
      mutate(year = str_split(keyName, "_", simplify = TRUE)[,2],
             individual_id = str_split(keyName, "_", simplify = TRUE)[,1]) |> 
      mutate(id.year = paste(individual_id,year,sep="_")) |>  
      left_join(mean.coord, by = "id.year") |> 
      dplyr::select(individual_id,year,dmax12m,mean.x, mean.y) }

# Final NULL checks
if (is.null(dmax12m) || nrow(dmax12m) == 0) return(NULL)

if(is.null(dmax12m)) NULL else {
  # Spatial annotation 10km
  cell_info_10 <- dgGEO_to_SEQNUM(dggs.10, mean.coord$mean.x, mean.coord$mean.y)
  dmax12m$grid.id.10km <- cell_info_10$seqnum
  centers_10 <- dgSEQNUM_to_GEO(dggs.10, dmax12m$grid.id.10km)
  dmax12m$lon.10km <- centers_10$lon_deg
  dmax12m$lat.10km <- centers_10$lat_deg
  
  # Spatial annotation 1km
  cell_info_1 <- dgGEO_to_SEQNUM(dggs.1, mean.coord$mean.x, mean.coord$mean.y)
  dmax12m$grid.id.1km <- cell_info_1$seqnum
  centers_1 <- dgSEQNUM_to_GEO(dggs.1, dmax12m$grid.id.1km)
  dmax12m$lon.1km <- centers_1$lon_deg
  dmax12m$lat.1km <- centers_1$lat_deg
  
}
return(dmax12m)

rm(cell_info);rm(centers)
rm(cell_info_10);rm(centers_10)
rm(mean.coord);rm(locs7d.sf)
rm(moveObjSplitTime);rm(maxNetDispL);rm(locs7d)

}



## ----function to summarize Max12m Displacements------
f_sum.ind.dmax12m<-function(x)
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
