## ----Maximum 24hr displacement-------------------------------------------------------------

calc_dmax24h <- function(trk, 
                         dggs_10, 
                         dggs_1, 
                         min_daily_n = 12, 
                         min_ymd_n = 7) 
  {
  
locs1h <- trk %>% 
  mutate(ymd = as.character(format(as.Date(t_), "%Y-%m-%d"))) |> 
  mutate(id.ymd = paste(individual_id,ymd,sep="_")) |> 
  group_by(id.ymd) %>% 
  filter(n() >= min_daily_n) %>%
  ungroup() 

if (nrow(locs1h) == 0) return(NULL)

# Calculate mean coordinates per day
mean.coord <- locs1h |> group_by(id.ymd) |> 
  mutate(mean.x = mean(x_, na.rm=T),
         mean.y = mean(y_, na.rm=T)) |> 
  dplyr::select(id.ymd, mean.x, mean.y) |> distinct()

# Convert to sf and calculate max pairwise distance per day
locs1h.sf <- sf::st_as_sf(locs1h,
                          coords = c("x_", "y_"),
                          crs = 4326)

moveObjSplitTime <- split(locs1h.sf, locs1h$id.ymd)
maxNetDispL <- lapply(moveObjSplitTime, function(x){max(sf::st_distance(x))})
maxNetDisp <- do.call("rbind",maxNetDispL)

# Process results
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

  # Final NULL checks
  if (is.null(dmax24) || nrow(dmax24) == 0) return(NULL)
 
if(is.null(dmax24)) NULL else {
  # Spatial annotation 10km
  cell_info_10 <- dgGEO_to_SEQNUM(dggs.10, mean.coord$mean.x, mean.coord$mean.y)
  dmax24$grid.id.10km <- cell_info_10$seqnum
  centers_10 <- dgSEQNUM_to_GEO(dggs.10, dmax24$grid.id.10km)
  dmax24$lon.10km <- centers_10$lon_deg
  dmax24$lat.10km <- centers_10$lat_deg
  
  # Spatial annotation 1km
  cell_info_1 <- dgGEO_to_SEQNUM(dggs.1, mean.coord$mean.x, mean.coord$mean.y)
  dmax24$grid.id.1km <- cell_info_1$seqnum
  centers_1 <- dgSEQNUM_to_GEO(dggs.1, dmax24$grid.id.1km)
  dmax24$lon.1km <- centers_1$lon_deg
  dmax24$lat.1km <- centers_1$lat_deg
}

return(dmax24)

rm(moveObjSplitTime);rm(maxNetDispL);rm(maxNetDisp)
rm(mean.coord);rm(locs1h.sf);rm(cell_info_1);rm(cell_info_10);
rm(centers_10);rm(centers_1)

}

# ----summarize dmax24h at individual level
f_sum.ind.dmax24h<-function(x)
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


# ----summarize dmax24h at individual level
f_sum.monthly.ind.dmax24h<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,month = NA,
                       year = NA,n.dmax24h.days = NA,
                       dmax24h.mean = NA,dmax24h.median = NA,dmax24h.cv = NA,dmax24h.95 = NA,dmax24h.05 = NA)
  } else {
    
    # derive month and year from t_
    ym <- substr(x$ymd, 1, 7)  
    
    # group index: individual x month
    id_ym <- interaction(x$individual_id, ym, drop = TRUE)
    
    individual_id <- tapply(as.character(x$individual_id), id_ym, unique)
    ym_grp        <- tapply(ym, id_ym, unique)
    
    # sample size per individual-month
    n.dmax24h.days <- as.numeric(tapply(x$ymd, id_ym, length))
    
    # 24hr Displacement
    dmax24h.mean<-as.numeric(with(x, tapply(x$dmax24h+0.001,id_ym, mean, na.rm=T)))
    dmax24h.median<-as.numeric(with(x, tapply(x$dmax24h+0.001,id_ym, median, na.rm=T)))
    dmax24h.cv<-as.numeric(with(x, tapply(x$dmax24h+0.001,id_ym, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
    dmax24h.95<-as.numeric(with(x, tapply(x$dmax24h+0.001,id_ym, quantile,.95, na.rm=T)))
    dmax24h.05<-as.numeric(with(x, tapply(x$dmax24h+0.001,id_ym, quantile,.05, na.rm=T)))
    
    ym <- as.character(unlist(ym_grp))
    ym_split <- do.call(rbind, strsplit(ym, "-", fixed = TRUE))
    year  <- as.integer(ym_split[, 1])
    month <- as.integer(ym_split[, 2])
    
    # build dataframe
    dats<-data.frame(individual_id,month,
                     year, n.dmax24h.days,
                     dmax24h.mean,dmax24h.median,dmax24h.cv,dmax24h.95,dmax24h.05)
    
    return(dats)
  }
}


