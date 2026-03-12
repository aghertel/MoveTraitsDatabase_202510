
## ----24hr displacement distance-------------------------------------------------------------

calc_d24h <- function(trk, 
                         dggs_10, 
                         dggs_1, 
                         dist_factor = 100000,
                         mindt = 20, 
                         maxdt = 28,
                         min_n = 30) 
{ 
  out <- trk %>% 
  mutate(d24h = step_lengths(.)*dist_factor) |> 
  mutate(time_diff = as.numeric(difftime(lead(t_),t_, units = "hours"))) |> 
  filter(time_diff >= mindt & time_diff <= maxdt) |> 
  dplyr::select(individual_id, t_, d24h, x_, y_) |> 
  mutate(ymd = as.character(format(as.Date(t_), "%Y-%m-%d")),
         month = lubridate::month(t_),
         year = lubridate::year(t_))|> 
  filter(!is.na(d24h))

  # drop if too short
  if (nrow(out) < min_n) return(NULL)
  
#if(is.null(out)) NULL else {
  # Spatial annotation 10km
  cell_info <- dgGEO_to_SEQNUM(dggs.10, out$x_, out$y_)
  out$grid.id.10km <- cell_info$seqnum
  centers <- dgSEQNUM_to_GEO(dggs.10, out$grid.id.10km)
  out$lon.10km <- centers$lon_deg
  out$lat.10km <- centers$lat_deg
  
  # Spatial annotation 1km
  cell_info <- dgGEO_to_SEQNUM(dggs.1, out$x_, out$y_)
  out$grid.id.1km <- cell_info$seqnum
  centers <- dgSEQNUM_to_GEO(dggs.1, out$grid.id.1km)
  out$lon.1km <- centers$lon_deg
  out$lat.1km <- centers$lat_deg
  
  out
}

## ----function to summarize 24h Displacements---------
f_sum.ind.d24h<-function(x)
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


## ----summarize d24h at monthly individual level---------
f_sum.monthly.ind.d24h<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    return(data.frame(individual_id = NA,n24h.days = NA,month= NA,
                       year= NA,d24h.mean = NA,d24h.median = NA,d24h.cv = NA,d24h.95 = NA,d24h.05 = NA))
  } 
  
  # derive month and year from t_
  ym <- format(x$t_, "%Y-%m")  
  
  # group index: individual x month
  id_ym <- interaction(x$individual_id, ym, drop = TRUE)
  
  individual_id <- tapply(as.character(x$individual_id), id_ym, unique)
  ym_grp        <- tapply(ym, id_ym, unique)
  
  # sample size per individual-month
  n24h.days<-as.numeric(tapply(x$t_, id_ym, length))
    
    # 24hr Displacement
    d24h.mean <- as.numeric(with(x, tapply(x$d24h+0.001,id_ym, mean, na.rm=T)))
    d24h.median <- as.numeric(with(x, tapply(x$d24h+0.001,id_ym, median, na.rm=T)))
    d24h.cv <- as.numeric(with(x, tapply(x$d24h+0.001,id_ym, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
    d24h.95 <- as.numeric(with(x, tapply(x$d24h+0.001,id_ym, quantile,.95, na.rm=T)))
    d24h.05 <- as.numeric(with(x, tapply(x$d24h+0.001,id_ym, quantile,.05, na.rm=T)))
    
    year  <- as.numeric(sub(".*\\.(\\d{4})-\\d{2}$", "\\1", unique(id_ym)))
    month <- as.numeric(sub(".*\\.\\d{4}-(\\d{2})$", "\\1", unique(id_ym)))
    
    # build dataframe
    dats<-data.frame(individual_id,year,month,n24h.days,
                     d24h.mean,d24h.median,d24h.cv,d24h.95,d24h.05)

}



