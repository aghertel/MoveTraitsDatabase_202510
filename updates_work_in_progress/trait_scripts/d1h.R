
## ----1h displacement-------------------------------------------------------------
calc_d1h <- function(trk,
                     dggs_10,
                     dggs_1,
                     min_n = 167,
                     dist_factor = 100000,
                     mindt = 45,
                     maxdt = 75) {
  out <- trk %>% 
    mutate(d1h = step_lengths(.) * dist_factor)  |>  
    mutate(time_diff = as.numeric(difftime(lead(t_), t_, units = "mins")))  |>  
    filter(time_diff >= mindt & time_diff <= maxdt) |>  
    filter(!is.na(d1h)) |>     
    mutate(t_ = floor_date(t_, "hour"),
           ymd = as.character(format(as.Date(t_), "%Y-%m-%d")),
           individual_id = as.character(individual_id)) |>
    dplyr::select(individual_id, t_, ymd, d1h, x_, y_) |>  
    rename("timestamp" = "t_",
           "lon" = "x_",
           "lat" = "y_") 
    
    
  
  # drop if too short
  if (nrow(out) < min_n) return(NULL)
  
  # spatial annotation 10 km
  cell_info_10 <- dgGEO_to_SEQNUM(dggs_10, out$lon, out$lat)
  out$grid.id.10km <- cell_info_10$seqnum
  centers_10 <- dgSEQNUM_to_GEO(dggs_10, out$grid.id.10km)
  out$lon.10km <- centers_10$lon_deg
  out$lat.10km <- centers_10$lat_deg
  
  # spatial annotation 1 km
  cell_info_1 <- dgGEO_to_SEQNUM(dggs_1, out$lon, out$lat)
  out$grid.id.1km <- cell_info_1$seqnum
  centers_1 <- dgSEQNUM_to_GEO(dggs_1, out$grid.id.1km)
  out$lon.1km <- centers_1$lon_deg
  out$lat.1km <- centers_1$lat_deg
  
  out
}


# ----summarize d1h at individual level
f_sum.ind.d1h <- function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n1h = NA,
                       d1h.mean = NA,d1h.median = NA,d1h.cv = NA,d1h.95 = NA,d1h.05 = NA)
  } else {
    individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
    
    # Get sample size per individual
    n1h<-as.numeric(with(x, tapply(as.character(x$timestamp),individual_id, length)))
    
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


# ----summarize d1h at monthly individual level
f_sum.monthly.ind.d1h <- function(x) {
  # If input is NULL, return a 1-row NA template (including month)
  if (is.null(x)) {
    return(data.frame(
      individual_id = NA,
      month         = NA,
      year          = NA,
      n1h           = NA,
      d1h.mean      = NA,
      d1h.median    = NA,
      d1h.cv        = NA,
      d1h.95        = NA,
      d1h.05        = NA
    ))
  }
  
  # derive month and year from t_
  ym <- format(x$timestamp, "%Y-%m")  
  
  # group index: individual x month
  id_ym <- interaction(x$individual_id, ym, drop = TRUE)
  
  individual_id <- tapply(as.character(x$individual_id), id_ym, unique)
  ym_grp        <- tapply(ym, id_ym, unique)
  
  # sample size per individual-month
  n1h <- as.numeric(tapply(x$timestamp, id_ym, length))
  
  # 1h displacement summaries per individual-month
  d1h.mean   <- as.numeric(tapply(x$d1h + 0.001, id_ym, mean,    na.rm = TRUE))
  d1h.median <- as.numeric(tapply(x$d1h + 0.001, id_ym, median,  na.rm = TRUE))
  d1h.cv     <- as.numeric(tapply(x$d1h + 0.001, id_ym, function(z) sd(z, na.rm = TRUE) / mean(z, na.rm = TRUE)))
  d1h.95     <- as.numeric(tapply(x$d1h + 0.001, id_ym, quantile, 0.95, na.rm = TRUE))
  d1h.05     <- as.numeric(tapply(x$d1h + 0.001, id_ym, quantile, 0.05, na.rm = TRUE))
  
  year  <- as.numeric(sub(".*\\.(\\d{4})-\\d{2}$", "\\1", unique(id_ym)))
  month <- as.numeric(sub(".*\\.\\d{4}-(\\d{2})$", "\\1", unique(id_ym)))
  
  data.frame(
    individual_id = unlist(individual_id),
    month         = month,
    year          = year,
    n1h           = n1h,
    d1h.mean      = d1h.mean,
    d1h.median    = d1h.median,
    d1h.cv        = d1h.cv,
    d1h.95        = d1h.95,
    d1h.05        = d1h.05,
    row.names     = NULL
  )
}


