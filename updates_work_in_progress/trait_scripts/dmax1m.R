## ----Maximum monthly displacement distance-------------------------------------------------------------
calc_dmax1m <- function(trk, 
                         dggs_10, 
                         dggs_1, 
                         min_monthly_n = 18 # 20 days? hours?
                         ) 
{
  locs24h <- trk |> 
  mutate(month = as.numeric(strftime(t_,format="%m")), 
         year = as.numeric(strftime(t_,format="%Y")),
         year_month = paste(year,month,sep="_")) |> 
  mutate(id.month = paste(individual_id,year_month,sep="_")) %>% 
  group_by(id.month) |> filter(n() >= min_monthly_n) |> ungroup()

if (nrow(locs24h) == 0) return(NULL)

# Convert to sf and calculate max pairwise distance per week
locs24h.sf <- sf::st_as_sf(locs24h,
                           coords = c("x_", "y_"),
                           crs = 4326)

dmax1m <- split(locs24h.sf, locs24h.sf$id.month) |>
  imap_dfr(function(x, nm) {
    d <- st_distance(x)
    diag(d) <- NA
    idx <- which(d == max(d, na.rm = TRUE), arr.ind = TRUE)[1, ]
    
    coords <- st_coordinates(x)
    
    tibble(
      id.month = nm,
      dmax1m = as.numeric(d[idx[1], idx[2]]),
      lon_start = coords[idx[1], "X"],
      lat_start = coords[idx[1], "Y"],
      lon_end = coords[idx[2], "X"],
      lat_end = coords[idx[2], "Y"]
    )
  })

dmax1m <- 
  dmax1m |> 
  mutate(month = as.numeric(str_split(id.month, "_", simplify = TRUE)[,3]),
         year = as.numeric(str_split(id.month, "_", simplify = TRUE)[,2]),
         individual_id = str_split(id.month, "_", simplify = TRUE)[,1]) |> 
      filter(!is.na(dmax1m)) |> 
  dplyr::select(individual_id,month, year,dmax1m, lon_start, lat_start, lon_end, lat_end)

# Final NULL checks
if (is.null(dmax1m) || nrow(dmax1m) == 0) return(NULL)

if(is.null(dmax1m)) NULL else {
  # Spatial annotation 10km
  cell_info_10.a <- dgGEO_to_SEQNUM(dggs.10, dmax1m$lon_start, dmax1m$lat_start)
  cell_info_10.b <- dgGEO_to_SEQNUM(dggs.10, dmax1m$lon_end, dmax1m$lat_end)
  dmax1m$grid.id.10km <- cell_info_10.a$seqnum
  dmax1m$grid.id.10km <- paste(dmax1m$grid.id.10km,cell_info_10.b$seqnum,sep=";")
  
  # Spatial annotation 1km
  cell_info_1.a <- dgGEO_to_SEQNUM(dggs.1, dmax1m$lon_start, dmax1m$lat_start)
  cell_info_1.b <- dgGEO_to_SEQNUM(dggs.1, dmax1m$lon_end, dmax1m$lat_end)
  dmax1m$grid.id.1km <- cell_info_1.a$seqnum
  dmax1m$grid.id.1km <- paste(dmax1m$grid.id.1km,cell_info_1.b$seqnum,sep=";")
  
rm(locs24h.sf)

}
return(dmax1m)
}

## ----function to summarize Max1m Displacements-------
f_sum.ind.dmax1m<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n.dmax1m.months = NA,
                       dmax1m.mean = NA,dmax1m.median = NA,dmax1m.cv = NA,dmax1m.95 = NA,dmax1m.05 = NA)
  } else {
    
    individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
    
    # Get sample size per individual
    n.dmax1m.months<-as.numeric(with(x, tapply(x$month,individual_id, length)))
    
    # 24hr Displacement
    dmax1m.mean<-as.numeric(with(x, tapply(x$dmax1m+0.001,individual_id, mean, na.rm=T)))
    dmax1m.median<-as.numeric(with(x, tapply(x$dmax1m+0.001,individual_id, median, na.rm=T)))
    dmax1m.cv<-as.numeric(with(x, tapply(x$dmax1m+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
    dmax1m.95<-as.numeric(with(x, tapply(x$dmax1m+0.001,individual_id, quantile,.95, na.rm=T)))
    dmax1m.05<-as.numeric(with(x, tapply(x$dmax1m+0.001,individual_id, quantile,.05, na.rm=T)))
    
    # build dataframe
    dats<-data.frame(individual_id,n.dmax1m.months,
                     dmax1m.mean,dmax1m.median,dmax1m.cv,dmax1m.95,dmax1m.05)
    
    return(dats)
  }
}

## ----summarize max1m at monthly individual level---------
f_sum.monthly.ind.dmax1m<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    return(data.frame(individual_id = NA,month= NA,
                      year= NA, iou1m = NA ))
  } 
  
  individual_id <- x$individual_id
  year <- x$year
  month <- x$month
  dmax1m <- x$dmax1m
  
  # build dataframe
  dats<-data.frame(individual_id,year,month,dmax1m)
  
}
