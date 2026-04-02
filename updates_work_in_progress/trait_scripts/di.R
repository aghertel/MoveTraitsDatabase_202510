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


calc_di <- function(trk, 
                        dggs_10, 
                        dggs_1, 
                        min_hours_n = 12) 
{

di <- 
  if(is.null(trk) ) NULL else {
    
    diurn.ind <- trk |> 
      mutate(date=as.Date(timestamp)) |> 
      dplyr::select(lon,lat,date,individual_id,timestamp,ymd,d1h) 
    
    sun_all <- getSunlightTimes(data = diurn.ind, 
                                tz = "UTC",
                                keep = c("sunrise", "sunset"))
    
    diurn.ind$sunrise <- sun_all$sunrise
    diurn.ind$sunset<- sun_all$sunset
    
    diurn.ind$daytime <- 
      ifelse(diurn.ind$timestamp < diurn.ind$sunrise |
               diurn.ind$timestamp > diurn.ind$sunset, "night","day")
    
    mean.coord <- diurn.ind |> 
      mutate(id.day = paste(individual_id,ymd,sep=".")) |> group_by(id.day) |> 
      mutate(mean.x = mean(lon, na.rm=T),
             mean.y = mean(lat, na.rm=T)) |> 
      dplyr::select(id.day, mean.x, mean.y) |> distinct()
    
    di <- diurn.ind %>%
      data.frame %>% 
      group_by(individual_id, ymd)%>%
      summarise(dist.sum.day = sum(d1h[daytime=="day"]),
                dist.sum.night = sum(d1h[daytime=="night"]),
                dist.length.day = sum(daytime=="day"),
                dist.length.night = sum(daytime=="night"),
                total.daylength = dist.length.day + dist.length.night) |> 
      data.frame()
    
    di$diurnality <- f.diurn(di$dist.sum.day, di$dist.length.day, di$dist.sum.night, 
                             di$dist.length.night)
    
    di <- di |> 
      mutate(id.day = paste(individual_id,ymd,sep=".")) |> 
      left_join(mean.coord,by = "id.day") |> 
      filter(total.daylength > 19)  |> 
      filter(!is.na(diurnality)) |> 
      dplyr::select( individual_id,ymd,diurnality,mean.x,mean.y)
    
    di <- if(nrow(di)==0) NULL else di
  }

if(is.null(di)) NULL else {
  # Spatial annotation 10km
  cell_info_10 <- dgGEO_to_SEQNUM(dggs.10, di$mean.x, di$mean.y)
  di$grid.id.10km <- cell_info_10$seqnum
  centers_10 <- dgSEQNUM_to_GEO(dggs.10, di$grid.id.10km)
  di$lon.10km <- centers_10$lon_deg
  di$lat.10km <- centers_10$lat_deg
  
  # Spatial annotation 1km
  cell_info_1 <- dgGEO_to_SEQNUM(dggs.1, di$mean.x, di$mean.y)
  di$grid.id.1km <- cell_info_1$seqnum
  centers_1 <- dgSEQNUM_to_GEO(dggs.1, di$grid.id.1km)
  di$lon.1km <- centers_1$lon_deg
  di$lat.1km <- centers_1$lat_deg
 
}
rm(cell_info_10);rm(centers_10);rm(cell_info_1);rm(centers_1);rm(mean.coord)
return(di)
}

## ----function to summarize Diurnality Index--------------------
f_sum.ind.di<-function(x)
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


## ----function for monthly summaries of Diurnality Index--------------------
f_sum.monthly.ind.di<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,year = NA,month=NA,n.di.days = NA,
                       di.mean = NA,di.median = NA,di.cv = NA,di.95 = NA,di.05 = NA)
  } else {

    # derive month and year from t_
    ym <- substr(x$ymd, 1, 7)
    
    # group index: individual x month
    id_ym <- interaction(x$individual_id, ym, drop = TRUE)
    
    individual_id <- tapply(as.character(x$individual_id), id_ym, unique)
    ym_grp        <- tapply(ym, id_ym, unique)
    
    # How many days per indivindividual_idual
    n.di.days<-as.numeric(with(x, tapply(x$ymd,id_ym, length)))
    
    # Diurnality
    di.mean<-as.numeric(with(x, tapply(x$diurnality,id_ym, mean, na.rm=T)))
    di.median<-as.numeric(with(x, tapply(x$diurnality,id_ym, median, na.rm=T)))
    di.cv<-as.numeric(with(x, tapply(x$diurnality,id_ym, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
    di.95<-as.numeric(with(x, tapply(x$diurnality,id_ym, quantile,.95, na.rm=T)))
    di.05<-as.numeric(with(x, tapply(x$diurnality,id_ym, quantile,.05, na.rm=T)))
    
    year  <- as.numeric(sub(".*\\.(\\d{4})-\\d{2}$", "\\1", unique(id_ym)))
    month <- as.numeric(sub(".*\\.\\d{4}-(\\d{2})$", "\\1", unique(id_ym)))
    
    # build dataframe
    dats<-data.frame(individual_id,n.di.days,year,month,
                     di.mean,di.median,di.cv,di.95,di.05)
    
    return(dats)
    rm(year);rm(month);rm(ym);rm(id_ym);rm(ym_grp);rm(dats);rm(individual_id)
    
  }
}

