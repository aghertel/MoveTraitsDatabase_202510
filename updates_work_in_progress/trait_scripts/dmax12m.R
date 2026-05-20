
## ----Maximum annual displacement distance-------------------------------------------------------------
# Based on daily relocations we calculated the maximum annual displacement 

# the animal hast to have data recorded in at least 9 months of the year
# minimum daily locations - 100

calc_dmax12m <- function(trk, 
                         dggs_10, 
                         dggs_1, 
                         min_locs_n = 150,
                         min_months = 9) 
{
  locs24h <- trk %>%
  mutate(year = as.numeric(strftime(t_,format="%Y")),
         month = as.numeric(strftime(t_,format="%m")),
         id.year = paste(individual_id,year,sep="_")) |> 
  group_by(individual_id, year) |> mutate(n_months = n_distinct(month)) |> 
  ungroup() |>   filter(n_months >= min_months) |>
  group_by(id.year)  |> filter(n() >= min_locs_n) |> ungroup() |> 
  dplyr::select(-c(n_months))

if (nrow(locs24h) == 0) return(NULL)

# Convert to sf and calculate max pairwise distance per day
locs24h.sf <- sf::st_as_sf(locs24h,
                          coords = c("x_", "y_"),
                          crs = 4326)

dmax12m <- split(locs24h.sf, locs24h.sf$id.year) |>
  imap_dfr(function(x, nm) {
    d <- st_distance(x)
    diag(d) <- NA
    idx <- which(d == max(d, na.rm = TRUE), arr.ind = TRUE)[1, ]
    
    coords <- st_coordinates(x)
    
    tibble(
      id.year = nm,
      dmax12m = as.numeric(d[idx[1], idx[2]]),
      lon_start = coords[idx[1], "X"],
      lat_start = coords[idx[1], "Y"],
      lon_end = coords[idx[2], "X"],
      lat_end = coords[idx[2], "Y"]
    )
  })

dmax12m <- 
  dmax12m |> 
  mutate(year = str_split(id.year, "_", simplify = TRUE)[,2],
         individual_id = str_split(id.year, "_", simplify = TRUE)[,1]) |> 
  filter(!is.na(dmax12m)) |> 
  dplyr::select(individual_id,year,dmax12m, lon_start, lat_start, lon_end, lat_end)

# Final NULL checks
if (is.null(dmax12m) || nrow(dmax12m) == 0) return(NULL)

if(is.null(dmax12m)) NULL else {
  # Spatial annotation 10km
  cell_info_10.a <- dgGEO_to_SEQNUM(dggs.10, dmax12m$lon_start, dmax12m$lat_start)
  cell_info_10.b <- dgGEO_to_SEQNUM(dggs.10, dmax12m$lon_end, dmax12m$lat_end)
  dmax12m$grid.id.10km <- cell_info_10.a$seqnum
  dmax12m$grid.id.10km <- paste(dmax12m$grid.id.10km,cell_info_10.b$seqnum,sep=";")
  
  # Spatial annotation 1km
  cell_info_1.a <- dgGEO_to_SEQNUM(dggs.1, dmax12m$lon_start, dmax12m$lat_start)
  cell_info_1.b <- dgGEO_to_SEQNUM(dggs.1, dmax12m$lon_end, dmax12m$lat_end)
  dmax12m$grid.id.1km <- cell_info_1.a$seqnum
  dmax12m$grid.id.1km <- paste(dmax12m$grid.id.1km,cell_info_1.b$seqnum,sep=";")
  
  rm(locs24h.sf)
  
}

return(dmax12m)

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
