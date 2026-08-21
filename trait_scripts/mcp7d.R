crs_4326  <- sp::CRS("EPSG:4326")
crs_moll <- sp::CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")

## ----Function to get corner coordinates (vertices)-------------------------------------------------------------
get_polygon_vertices <- function(mcp_spdf) {
  out <- lapply(seq_along(mcp_spdf@polygons), function(i) {
    coords <- mcp_spdf@polygons[[i]]@Polygons[[1]]@coords
    if (all(coords[1, ] == coords[nrow(coords), ])) {
      coords <- coords[-nrow(coords), , drop = FALSE]
    }
    
    tmp <- SpatialPoints(coords, proj4string = crs_moll)
    tmp_ll <- spTransform(tmp, crs_4326)
    coords_ll <- coordinates(tmp_ll)
    
    data.frame(
      id = mcp_spdf@data$id[i],
      x_vertices = paste(coords_ll[, 1], collapse = ";"),
      y_vertices = paste(coords_ll[, 2], collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

## ----Weekly MCP-------------------------------------------------------------

calc_mcp7d <- function(trk, 
                        dggs_10, 
                        dggs_1, 
                        min_hours_n = 84) 
{
dat.mcp.weekly <- trk %>% 
  mutate(week = as.numeric(strftime(t_,format="%W")), 
         year = as.numeric(strftime(t_,format="%Y")),
         year_week = paste(year,week, sep="_")) %>% 
  mutate(id.week = paste(individual_id,year_week,sep=".")) %>% 
  filter(!is.na(x_)) %>% filter(!is.na(y_)) %>%
  group_by(id.week) %>% filter(n() > min_hours_n) %>% ungroup() %>% 
  dplyr::select(x_,y_,id.week) 

tryCatch({
  coordinates(dat.mcp.weekly) <- c("x_","y_")
  proj4string(dat.mcp.weekly) <- crs_4326
  dat.mcp.weekly <- spTransform(dat.mcp.weekly,
                                crs_moll)
}, error = function(e) {NA})

mcp.weekly <- if (nrow(dat.mcp.weekly) == 0) {
  NULL
} else {
  mcp_spdf <- mcp(dat.mcp.weekly, percent = 95, unout = c("m2"))
  
  mcp_df <- data.frame(mcp_spdf) |>
    mutate(
      week = as.numeric(stringr::str_extract(id, "(\\d+$)")),
      year_week = stringr::str_extract(id, "[^.]*$"),
      individual_id = str_extract(id, "[^.]+")) |>
    filter(!is.na(area))
  
  verts_df <- get_polygon_vertices(mcp_spdf)
  
  mcp_df |>
    left_join(verts_df, by = "id")
}

if (is.null(mcp.weekly)) {
  NULL
} else {
  vertices_long <- mcp.weekly |>
    dplyr::select(id, individual_id, week, year_week, area, x_vertices, y_vertices) |>
    mutate(row_id = row_number()) |>
    separate_rows(x_vertices, y_vertices, sep = ";") |>
    mutate(
      x_vertices = as.numeric(x_vertices),
      y_vertices = as.numeric(y_vertices)
    )
  
  grid_10 <- dgGEO_to_SEQNUM(dggs.10, vertices_long$x_vertices, vertices_long$y_vertices)$seqnum
  grid_1  <- dgGEO_to_SEQNUM(dggs.1,  vertices_long$x_vertices, vertices_long$y_vertices)$seqnum
  
  vertices_long <- vertices_long |>
    mutate(
      grid.id.10km = grid_10,
      grid.id.1km = grid_1
    )
  
  grid_summary <- vertices_long |>
    group_by(id, individual_id, week, year_week, area) |>
    summarise(
      x_vertices = paste(x_vertices, collapse = ";"),
      y_vertices = paste(y_vertices, collapse = ";"),
      grid.id.10km = paste(sort(unique(grid.id.10km)), collapse = ";"),
      grid.id.1km = paste(sort(unique(grid.id.1km)), collapse = ";"),
      .groups = "drop"
    )
  
  mcp.weekly <- grid_summary
}
return(mcp.weekly)
}

## ----function to summarize 7d MCP--------------------
f_sum.ind.mcp7d<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n.mcp7d.weeks = NA,
                       mcp7d.mean = NA,mcp7d.median = NA,mcp7d.cv = NA,mcp7d.95 = NA,mcp7d.05 = NA)
  } else {
    
    individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
    
    # Get sample size per indivindividual_idual
    n.mcp7d.weeks<-as.numeric(with(x, tapply(x$year_week,individual_id, length)))
    
    # 24mcp Displacement
    mcp7d.mean<-as.numeric(with(x, tapply(x$area+0.001,individual_id, mean, na.rm=T)))
    mcp7d.median<-as.numeric(with(x, tapply(x$area+0.001,individual_id, median, na.rm=T)))
    mcp7d.cv<-as.numeric(with(x, tapply(x$area+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
    mcp7d.95<-as.numeric(with(x, tapply(x$area+0.001,individual_id, quantile,.95, na.rm=T)))
    mcp7d.05<-as.numeric(with(x, tapply(x$area+0.001,individual_id, quantile,.05, na.rm=T)))
    
    # build dataframe
    dats<-data.frame(individual_id,n.mcp7d.weeks,
                     mcp7d.mean,mcp7d.median,mcp7d.cv,mcp7d.95,mcp7d.05)
    
    return(dats)
  }
}


## ----summarize mcp7d at monthly individual level---------
f_sum.monthly.ind.mcp7d<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,year=NA, month=NA,n.mcp7d.weeks = NA,
                       mcp7d.mean = NA,mcp7d.median = NA,mcp7d.cv = NA,mcp7d.95 = NA,mcp7d.05 = NA)
  } else {
    
    # derive month and year from t_
    year <- as.numeric(substr(x$year_week, 1, 4))  
    
    string_iso <- paste(year, sprintf("W%02d", as.numeric(x$week)), 1, sep="-")
    date <- ISOweek2date(string_iso)
    month <- format(date, "%m")
    ym <- paste(year,month,sep="-")
    
    # group index: individual x month
    id_ym <- interaction(x$individual_id, ym, drop = TRUE)
    
    individual_id <- tapply(as.character(x$individual_id), id_ym, unique)
    ym_grp        <- tapply(ym, id_ym, unique)
    
    # sample size per individual-month
    n.mcp7d.weeks <- as.numeric(tapply(x$year_week, id_ym, length))
    
    # 24hr Displacement
    mcp7d.mean<-as.numeric(with(x, tapply(x$area+0.001,id_ym, mean, na.rm=T)))
    mcp7d.median<-as.numeric(with(x, tapply(x$area+0.001,id_ym, median, na.rm=T)))
    mcp7d.cv<-as.numeric(with(x, tapply(x$area+0.001,id_ym, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
    mcp7d.95<-as.numeric(with(x, tapply(x$area+0.001,id_ym, quantile,.95, na.rm=T)))
    mcp7d.05<-as.numeric(with(x, tapply(x$area+0.001,id_ym, quantile,.05, na.rm=T)))
    
    year  <- as.numeric(sub(".*\\.(\\d{4})-\\d{2}$", "\\1", unique(id_ym)))
    month <- as.numeric(sub(".*\\.\\d{4}-(\\d{2})$", "\\1", unique(id_ym)))
    
    # build dataframe
    dats<-data.frame(individual_id,year,month,n.mcp7d.weeks,
                     mcp7d.mean,mcp7d.median,mcp7d.cv,mcp7d.95,mcp7d.05)
    
    return(dats)
  }
}
