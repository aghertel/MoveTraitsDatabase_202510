## ----Maximum 7day displacement distance-------------------------------------------------------------
calc_dmax7d <- function(trk, 
                         dggs_10, 
                         dggs_1, 
                         min_weekly_n = 5, 
                         min_week_n = 10) 
{

  locs24h <- trk |> 
  mutate(week = as.numeric(strftime(t_,format="%W")), 
         year = as.numeric(strftime(t_,format="%Y")),
         year_week = paste(year,week, sep="_")) |> 
  mutate(id.year_week = paste(individual_id,year_week,sep="_")) |>  
  group_by(individual_id,year_week) |> filter(n() >= min_weekly_n) |> ungroup()

if (nrow(locs24h) == 0) return(NULL)

# Calculate mean coordinates per week
mean.coord <- locs24h |> group_by(id.year_week) |> 
  mutate(mean.x = mean(x_, na.rm=T),
         mean.y = mean(y_, na.rm=T)) |> 
  dplyr::select(id.year_week, mean.x, mean.y) |> distinct()

# Convert to sf and calculate max pairwise distance per week
locs24h.sf <- sf::st_as_sf(locs24h,
                           coords = c("x_", "y_"),
                           crs = 4326)

moveObjSplitTime <- split(locs24h.sf, locs24h.sf$id.year_week)
maxNetDispL <- lapply(moveObjSplitTime, function(x){max(sf::st_distance(x))})
maxNetDisp <- do.call("rbind",maxNetDispL)

# Process results
dmax7d <- 
  if(is.null(maxNetDisp)) NULL else {
    data.frame(keyName=row.names(maxNetDisp), dmax7d=maxNetDisp[,1], row.names=NULL) |> 
      mutate(week = str_split(keyName, "_", simplify = TRUE)[,3],
             year_week = paste(str_split(keyName, "_", simplify = TRUE)[,2],
                               str_split(keyName, "_", simplify = TRUE)[,3], sep ="_"),
             individual_id = str_split(keyName, "_", simplify = TRUE)[,1]) |> 
      filter(!is.na(dmax7d))|> 
      group_by(individual_id) |>  filter(n() >= min_week_n) |> ungroup() |>  
      mutate(id.year_week = paste(individual_id,year_week,sep="_")) |> 
      left_join(mean.coord, by = "id.year_week") |> 
      dplyr::select(individual_id,week, year_week,dmax7d,mean.x, mean.y) }

# Final NULL checks
if (is.null(dmax7d) || nrow(dmax7d) == 0) return(NULL)

if(is.null(dmax7d)) NULL else {
  # Spatial annotation 10km
  cell_info_10 <- dgGEO_to_SEQNUM(dggs.10, mean.coord$mean.x, mean.coord$mean.y)
  dmax7d$grid.id.10km <- cell_info_10$seqnum
  centers_10 <- dgSEQNUM_to_GEO(dggs.10, dmax7d$grid.id.10km)
  dmax7d$lon.10km <- centers_10$lon_deg
  dmax7d$lat.10km <- centers_10$lat_deg
  
  
  # Spatial annotation 1km
  cell_info_1 <- dgGEO_to_SEQNUM(dggs.1, mean.coord$mean.x, mean.coord$mean.y)
  dmax7d$grid.id.1km <- cell_info_1$seqnum
  centers_1 <- dgSEQNUM_to_GEO(dggs.1, dmax7d$grid.id.1km)
  dmax7d$lon.1km <- centers_1$lon_deg
  dmax7d$lat.1km <- centers_1$lat_deg

  rm(moveObjSplitTime);rm(maxNetDispL)
  rm(mean.coord);rm(locs24h.sf)
  rm(cell_info_10);rm(centers_10);rm(cell_info_1);rm(centers_1)
  
}
return(dmax7d)
}

## ----function to summarize Max7d Displacements-------
f_sum.ind.dmax7d<-function(x)
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


## ----function to summarize Max7d Displacements-------
f_sum.monthly.ind.dmax7d<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,year = NA,month=NA,n.dmax7d.weeks = NA,
                       dmax7d.mean = NA,dmax7d.median = NA,dmax7d.cv = NA,dmax7d.95 = NA,dmax7d.05 = NA)
  } else {
    
    # derive month and year from t_
    y <- as.numeric(substr(x$year_week, 1, 4))  
    string_iso <- paste(y, sprintf("W%02d", as.numeric(x$week)), 1, sep="-")
    date <- ISOweek2date(string_iso)
    m <- format(date, "%m")
    ym <- paste(y,m,sep="-")
    
    # group index: individual x month
    id_ym <- interaction(x$individual_id, ym, drop = TRUE)
    
    individual_id <- tapply(as.character(x$individual_id), id_ym, unique)
    ym_grp        <- tapply(ym, id_ym, unique)
    
    # sample size per individual-month
    n.dmax7d.weeks <- as.numeric(tapply(x$year_week, id_ym, length))

    # 24hr Displacement
    dmax7d.mean<-as.numeric(with(x, tapply(x$dmax7d+0.001,id_ym, mean, na.rm=T)))
    dmax7d.median<-as.numeric(with(x, tapply(x$dmax7d+0.001,id_ym, median, na.rm=T)))
    dmax7d.cv<-as.numeric(with(x, tapply(x$dmax7d+0.001,id_ym, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
    dmax7d.95<-as.numeric(with(x, tapply(x$dmax7d+0.001,id_ym, quantile,.95, na.rm=T)))
    dmax7d.05<-as.numeric(with(x, tapply(x$dmax7d+0.001,id_ym, quantile,.05, na.rm=T)))
    
    year  <- as.numeric(sub(".*\\.(\\d{4})-\\d{2}$", "\\1", unique(id_ym)))
    month <- as.numeric(sub(".*\\.\\d{4}-(\\d{2})$", "\\1", unique(id_ym)))
    
    # build dataframe
    dats<-data.frame(individual_id,year,month,n.dmax7d.weeks,dmax7d.mean,dmax7d.median,dmax7d.cv,dmax7d.95,dmax7d.05)
    
    return(dats)
  }
  rm(year);rm(month);rm(string_iso);rm(date);rm(m);rm(ym);rm(y);rm(ym_grp);rm(dats);rm(individual_id)
}
