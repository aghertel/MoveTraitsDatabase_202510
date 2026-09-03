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

## ----Daily MCP-------------------------------------------------------------
calc_mcp24h <- function(trk, 
                         dggs_100, 
                         dggs_10, 
                         dggs_1, 
                         min_hours_n = 12) 
{
  dat.mcp.daily <- trk %>% 
  tibble() %>% mutate(ymd = as.character(format(as.Date(t_), "%Y-%m-%d")))  %>%
  filter(!is.na(x_)) %>% 
  filter(!is.na(y_)) %>% 
  mutate(id.day = paste(individual_id,ymd,sep="_")) %>% 
  group_by(id.day) %>% 
  filter(n() > min_hours_n) %>% 
  ungroup() %>% 
  dplyr::select(x_,y_,id.day) 

tryCatch({
  coordinates(dat.mcp.daily) <- c("x_","y_")
  proj4string(dat.mcp.daily) <- crs_4326
  dat.mcp.daily <- spTransform(dat.mcp.daily,crs_moll)
}, error = function(e) {NA})

mcp.daily <- if (nrow(dat.mcp.daily) == 0) {
  NULL
} else {
  mcp_spdf <- mcp(dat.mcp.daily, percent = 95, unout = c("m2"))
  
  mcp_df <- data.frame(mcp_spdf) |>
    mutate(
      ymd = str_extract(id, "\\d{4}-\\d{2}-\\d{2}$"),
      individual_id = stringr::str_remove(id, "_\\d{4}-\\d{2}-\\d{2}$")) |>
    filter(!is.na(area))
  
  verts_df <- get_polygon_vertices(mcp_spdf)
  
  mcp_df |>
    left_join(verts_df, by = "id")
}

if (is.null(mcp.daily)) {
  NULL
} else {
  vertices_long <- mcp.daily |>
    dplyr::select(id, individual_id, ymd, area, x_vertices, y_vertices) |>
    mutate(row_id = row_number()) |>
    separate_rows(x_vertices, y_vertices, sep = ";") |>
    mutate(
      x_vertices = as.numeric(x_vertices),
      y_vertices = as.numeric(y_vertices)
    )
  
  grid_100 <- dgGEO_to_SEQNUM(dggs_100, vertices_long$x_vertices, vertices_long$y_vertices)$seqnum
  grid_10 <- dgGEO_to_SEQNUM(dggs_10, vertices_long$x_vertices, vertices_long$y_vertices)$seqnum
  grid_1  <- dgGEO_to_SEQNUM(dggs_1,  vertices_long$x_vertices, vertices_long$y_vertices)$seqnum
  
  vertices_long <- vertices_long |>
    mutate(
      grid.id.100km = grid_100,
      grid.id.10km = grid_10,
      grid.id.1km = grid_1
    )
  
  grid_summary <- vertices_long |>
    group_by(id, individual_id, ymd, area) |>
    summarise(
      x_vertices = paste(x_vertices, collapse = ";"),
      y_vertices = paste(y_vertices, collapse = ";"),
      grid.id.100km = paste(sort(unique(grid.id.100km)), collapse = ";"),
      grid.id.10km = paste(sort(unique(grid.id.10km)), collapse = ";"),
      grid.id.1km = paste(sort(unique(grid.id.1km)), collapse = ";"),
      .groups = "drop"
    )
  
  mcp.daily <- grid_summary
}
return(mcp.daily)
}

## ----function to summarize 1d MCP--------------------
f_sum.ind.mcp24h<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA_character_,n.mcp24h.days = NA,
                       mcp24h.mean = NA,mcp24h.median = NA,mcp24h.cv = NA,mcp24h.95 = NA,mcp24h.05 = NA)
  } else {
    
    individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
    
    # Get sample size per indivindividual_idual
    n.mcp24h.days<-as.numeric(with(x, tapply(x$individual_id,individual_id, length)))
    
    # 24MCP Displacement
    mcp24h.mean<-as.numeric(with(x, tapply(x$area+0.001,individual_id, mean, na.rm=T)))
    mcp24h.median<-as.numeric(with(x, tapply(x$area+0.001,individual_id, median, na.rm=T)))
    mcp24h.cv<-as.numeric(with(x, tapply(x$area+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
    mcp24h.95<-as.numeric(with(x, tapply(x$area+0.001,individual_id, quantile,.95, na.rm=T)))
    mcp24h.05<-as.numeric(with(x, tapply(x$area+0.001,individual_id, quantile,.05, na.rm=T)))
    
    # build dataframe
    dats<-data.frame(individual_id,n.mcp24h.days,
                     mcp24h.mean,mcp24h.median,mcp24h.cv,mcp24h.95,mcp24h.05)
    
    return(dats)
  }
}


## ----summarize mcp24h at monthly individual level---------
f_sum.monthly.ind.mcp24h<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    return(data.frame(individual_id = NA_character_,year=NA,month=NA,n.mcp24h.days = NA,
                      mcp24h.mean = NA,mcp24h.median = NA,mcp24h.cv = NA,mcp24h.95 = NA,mcp24h.05 = NA))
  } 
  
  # derive month and year from t_
  ym <- substr(x$ymd, 1, 7)  
  
  # group index: individual x month
  id_ym <- interaction(x$individual_id, ym, drop = TRUE)
  
  individual_id <- tapply(as.character(x$individual_id), id_ym, unique)
  ym_grp        <- tapply(ym, id_ym, unique)
  
  # sample size per individual-month
  n.mcp24h.days<-as.numeric(tapply(x$ymd, id_ym, length))
  
  # 24hr Displacement
  mcp24h.mean <- as.numeric(with(x, tapply(x$area+0.001,id_ym, mean, na.rm=T)))
  mcp24h.median <- as.numeric(with(x, tapply(x$area+0.001,id_ym, median, na.rm=T)))
  mcp24h.cv <- as.numeric(with(x, tapply(x$area+0.001,id_ym, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
  mcp24h.95 <- as.numeric(with(x, tapply(x$area+0.001,id_ym, quantile,.95, na.rm=T)))
  mcp24h.05 <- as.numeric(with(x, tapply(x$area+0.001,id_ym, quantile,.05, na.rm=T)))
  
  year  <- as.numeric(sub(".*\\.(\\d{4})-\\d{2}$", "\\1", unique(id_ym)))
  month <- as.numeric(sub(".*\\.\\d{4}-(\\d{2})$", "\\1", unique(id_ym)))
  
  # build dataframe
  dats<-data.frame(individual_id,year,month,n.mcp24h.days,
                   mcp24h.mean,mcp24h.median,mcp24h.cv,mcp24h.95,mcp24h.05)
  
}


