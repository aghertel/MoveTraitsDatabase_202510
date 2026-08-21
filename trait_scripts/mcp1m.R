
crs_4326  <- sp::CRS("EPSG:4326")
crs_moll <- sp::CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")

## ----Function to get corner coordinates (vertices)-------------------------------------------------------------

get_polygon_vertices <- function(mcp_spdf) {
  if (is.null(mcp_spdf) || length(mcp_spdf@polygons) == 0) return(NULL)
  
  out <- lapply(seq_along(mcp_spdf@polygons), function(i) {
    poly <- mcp_spdf@polygons[[i]]
    if (length(poly@Polygons) == 0) return(NULL)
    
    coords <- poly@Polygons[[1]]@coords
    if (is.null(coords) || nrow(coords) < 3) return(NULL)
    
    if (all(coords[1, ] == coords[nrow(coords), ])) {
      coords <- coords[-nrow(coords), , drop = FALSE]
    }
    if (nrow(coords) < 3) return(NULL)
    
    tmp <- sp::SpatialPoints(coords, proj4string = crs_moll)
    if (length(tmp) == 0) return(NULL)
    
    tmp_ll <- tryCatch(
      sp::spTransform(tmp, crs_4326),
      error = function(e) NULL
    )
    if (is.null(tmp_ll)) return(NULL)
    
    coords_ll <- sp::coordinates(tmp_ll)
    if (is.null(coords_ll) || nrow(coords_ll) < 3) return(NULL)
    
    data.frame(
      id = mcp_spdf@data$id[i],
      x_vertices = paste(coords_ll[, 1], collapse = ";"),
      y_vertices = paste(coords_ll[, 2], collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  
  out <- Filter(Negate(is.null), out)
  if (length(out) == 0) return(NULL)
  
  do.call(rbind, out)
}

## ----Monthly MCP-------------------------------------------------------------
calc_mcp1m <- function(trk, 
                       dggs_10, 
                       dggs_1, 
                       min_days_n = 18) 
{
  dat.mcp.monthly <- trk %>% 
  tibble() %>% 
  mutate(month = as.numeric(strftime(t_,format="%m")), 
         year = as.numeric(strftime(t_,format="%Y"))) %>% 
  mutate(year_month = paste(year,month,sep="_")) |> 
  mutate(id.month = paste(individual_id,year_month,sep="_")) %>% 
  filter(!is.na(x_)) %>% filter(!is.na(y_)) %>% group_by(id.month) %>% 
  filter(n() > min_days_n) %>% ungroup() %>% dplyr::select(x_,y_,id.month) 
  
tryCatch({
  coordinates(dat.mcp.monthly) <- c("x_","y_")
  proj4string(dat.mcp.monthly) <- crs_4326
  dat.mcp.monthly <- spTransform(dat.mcp.monthly,crs_moll)
}, error = function(e) {NA})

mcp.monthly <- if (nrow(dat.mcp.monthly) == 0) {
    NULL
  } else {
    mcp_spdf <- mcp(dat.mcp.monthly, percent = 95, unout = c("m2"))
    
    mcp_df <- data.frame(mcp_spdf) |>
      mutate(
        month = as.numeric(stringr::str_extract(id, "\\d{1,2}$")),
        year_month = stringr::str_extract(id, "\\d{4}_\\d{1,2}$"),
        individual_id = stringr::str_remove(id, "_\\d{4}_\\d{1,2}$")) |>
      filter(!is.na(area))
    
    verts_df <- get_polygon_vertices(mcp_spdf)
    
    mcp_df |>
      left_join(verts_df, by = "id")
  }
  

if (is.null(mcp.monthly)) {
  NULL
} else {
  vertices_long <- mcp.monthly |>
    dplyr::select(id, individual_id, month, year_month, area, x_vertices, y_vertices) |>
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
    group_by(id, individual_id, month, year_month, area) |>
    summarise(
      x_vertices = paste(x_vertices, collapse = ";"),
      y_vertices = paste(y_vertices, collapse = ";"),
      grid.id.10km = paste(sort(unique(grid.id.10km)), collapse = ";"),
      grid.id.1km = paste(sort(unique(grid.id.1km)), collapse = ";"),
      .groups = "drop"
    )
  
  mcp.monthly <- grid_summary
}
}

## ----function to summarize 1m MCP--------------------
f_sum.ind.mcp1m<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA_character_,n.mcp1m.months = NA,
                       mcp1m.mean = NA,mcp1m.median = NA,mcp1m.cv = NA,mcp1m.95 = NA,mcp1m.05 = NA)
  } else {
    
    individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
    
    n.mcp1m.months<-as.numeric(with(x, tapply(x$year_month,individual_id, length)))
    
    mcp1m.mean<-as.numeric(with(x, tapply(x$area+0.001,individual_id, mean, na.rm=T)))
    mcp1m.median<-as.numeric(with(x, tapply(x$area+0.001,individual_id, median, na.rm=T)))
    mcp1m.cv<-as.numeric(with(x, tapply(x$area+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
    mcp1m.95<-as.numeric(with(x, tapply(x$area+0.001,individual_id, quantile,.95, na.rm=T)))
    mcp1m.05<-as.numeric(with(x, tapply(x$area+0.001,individual_id, quantile,.05, na.rm=T)))
    
    dats<-data.frame(individual_id,n.mcp1m.months,
                     mcp1m.mean,mcp1m.median,mcp1m.cv,mcp1m.95,mcp1m.05)
    
    return(dats)
  }
}


## ----summarize mcp1m at monthly individual level---------
f_sum.monthly.ind.mcp1m<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    return(data.frame(individual_id = NA_character_,month= NA,
                      year= NA, mcp1m = NA ))
  } 
  
  individual_id <- x$individual_id
  year <- as.numeric(substr(x$year_month, 1, 4))
  month <- x$month
  mcp1m <- x$area
  
  # build dataframe
  dats<-data.frame(individual_id,year,month,mcp1m)
  
}




