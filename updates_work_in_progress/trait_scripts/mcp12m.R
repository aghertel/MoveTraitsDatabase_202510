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
## ----Annual MCP-------------------------------------------------------------
calc_mcp12m <- function(trk, 
                       dggs_10, 
                       dggs_1, 
                       min_weeks_n = 36) 
{
  dat.mcp.annual <- trk %>% 
  tibble() |>   mutate(year = as.numeric(strftime(t_,format="%Y"))) %>% 
  mutate(id.year = paste(individual_id,year,sep=".")) %>% 
  filter(!is.na(x_)) %>% filter(!is.na(y_)) %>% group_by(id.year) %>% 
  filter(n() > min_weeks_n) %>% ungroup() %>% dplyr::select(x_,y_,id.year) 

mean.coord <- dat.mcp.annual |> group_by(id.year) |> 
  mutate(mean.x = mean(x_, na.rm=T),
         mean.y = mean(y_, na.rm=T)) |> 
  dplyr::select(id.year, mean.x, mean.y) |> distinct()

tryCatch({
  coordinates(dat.mcp.annual) <- c("x_","y_")
  proj4string(dat.mcp.annual) <- crs_4326  
  dat.mcp.annual <- spTransform(dat.mcp.annual,
                                crs_moll)
}, error = function(e) {NA})
# 
# mcp.annual <- 
#   if(nrow(dat.mcp.annual)==0) NULL else {
#     mcp(dat.mcp.annual, percent = 95, unout = c("m2")) %>% data.frame() %>% 
#       left_join(mean.coord, by = c("id" = "id.year")) |> 
#       mutate(id.year = id,
#              year = stringr::str_extract(id, "[^.]*$"),
#              individual_id = str_extract(id, "[^.]+")) %>% 
#       #dplyr::select(individual_id, id.year, year,area, mean.x, mean.y)|> 
#       filter(!is.na(area)) 
#   }

mcp.annual <- if (nrow(dat.mcp.annual) == 0) {
  NULL
} else {
  mcp_spdf <- mcp(dat.mcp.annual, percent = 95, unout = c("m2"))
  
  mcp_df <- data.frame(mcp_spdf) |>
    mutate(
      id.year = id,
      year = stringr::str_extract(id, "[^.]*$"),
      individual_id = str_extract(id, "[^.]+")) |>
    filter(!is.na(area))
  
  verts_df <- get_polygon_vertices(mcp_spdf)
  
  mcp_df |>
    left_join(verts_df, by = "id")
}

if (is.null(mcp.annual)) {
  NULL
} else {
  vertices_long <- mcp.annual |>
    dplyr::select(id, individual_id, year, area, x_vertices, y_vertices) |>
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
    group_by(id, individual_id, year, area) |>
    summarise(
      x_vertices = paste(x_vertices, collapse = ";"),
      y_vertices = paste(y_vertices, collapse = ";"),
      grid.id.10km = paste(sort(unique(grid.id.10km)), collapse = ";"),
      grid.id.1km = paste(sort(unique(grid.id.1km)), collapse = ";"),
      .groups = "drop"
    )
  
  mcp.annual <- grid_summary
}
# 
# if(is.null(mcp.annual)) NULL else {
#   # Spatial annotation 10km
#   mean.coord$grid.id.10km <- dgGEO_to_SEQNUM(dggs.10, mean.coord$mean.x, mean.coord$mean.y)$seqnum
#   mcp.annual <- mcp.annual |> left_join(mean.coord[,c("id.year","grid.id.10km")], by = c("id" = "id.year"))
#   centers_10 <- dgSEQNUM_to_GEO(dggs.10, mcp.annual$grid.id.10km)
#   mcp.annual$lon.10km <- centers_10$lon_deg
#   mcp.annual$lat.10km <- centers_10$lat_deg
#   
#   # Spatial annotation 1km
#   mean.coord$grid.id.1km <- dgGEO_to_SEQNUM(dggs.1, mean.coord$mean.x, mean.coord$mean.y)$seqnum
#   mcp.annual <- mcp.annual |> left_join(mean.coord[,c("id.year","grid.id.1km")], by = c("id" = "id.year"))
#   centers_1 <- dgSEQNUM_to_GEO(dggs.1, mcp.annual$grid.id.1km)
#   mcp.annual$lon.1km <- centers_1$lon_deg
#   mcp.annual$lat.1km <- centers_1$lat_deg
# 
#   
#   mcp.annual <- mcp.annual |> 
#     dplyr::select(individual_id,id.year, year, area, mean.x, mean.y, 
#                   grid.id.10km, lon.10km,  lat.10km, 
#                   grid.id.1km,  lon.1km,   lat.1km)
#}

return(mcp.annual)
}


## ----function to summarize 12m MCP--------------------
f_sum.ind.mcp12m<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n.mcp12m.years = NA,
                       mcp12m.mean = NA,mcp12m.median = NA,mcp12m.cv = NA,mcp12m.95 = NA,mcp12m.05 = NA)
  } else {
    
    individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
    
    n.mcp12m.years<-as.numeric(with(x, tapply(x$year,individual_id, length)))
    
    mcp12m.mean<-as.numeric(with(x, tapply(x$area+0.001,individual_id, mean, na.rm=T)))
    mcp12m.median<-as.numeric(with(x, tapply(x$area+0.001,individual_id, median, na.rm=T)))
    mcp12m.cv<-as.numeric(with(x, tapply(x$area+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
    mcp12m.95<-as.numeric(with(x, tapply(x$area+0.001,individual_id, quantile,.95, na.rm=T)))
    mcp12m.05<-as.numeric(with(x, tapply(x$area+0.001,individual_id, quantile,.05, na.rm=T)))
    
    dats<-data.frame(individual_id,n.mcp12m.years,
                     mcp12m.mean,mcp12m.median,mcp12m.cv,mcp12m.95,mcp12m.05)
    
    return(dats)
  }}
