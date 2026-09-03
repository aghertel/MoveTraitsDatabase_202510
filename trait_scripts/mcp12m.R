crs_4326  <- sp::CRS("EPSG:4326")
crs_moll <- sp::CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")

## ----Function to get corner coordinates (vertices)-------------------------------------------------------------
# get_polygon_vertices <- function(mcp_spdf) {
#   out <- lapply(seq_along(mcp_spdf@polygons), function(i) {
#     coords <- mcp_spdf@polygons[[i]]@Polygons[[1]]@coords
#     if (all(coords[1, ] == coords[nrow(coords), ])) {
#       coords <- coords[-nrow(coords), , drop = FALSE]
#     }
#     
#     tmp <- SpatialPoints(coords, proj4string = crs_moll)
#     tmp_ll <- spTransform(tmp, crs_4326)
#     coords_ll <- coordinates(tmp_ll)
#     
#     data.frame(
#       id = mcp_spdf@data$id[i],
#       x_vertices = paste(coords_ll[, 1], collapse = ";"),
#       y_vertices = paste(coords_ll[, 2], collapse = ";"),
#       stringsAsFactors = FALSE
#     )
#   })
#   do.call(rbind, out)
# }
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

## ----Annual MCP-------------------------------------------------------------
calc_mcp12m <- function(trk, 
                       dggs_100,
                       dggs_10, 
                       dggs_1, 
                       min_locs_n = 150,
                       min_months = 9) 
{
  out <- tryCatch({
    
  dat.mcp.annual <- trk %>% 
  tibble() |>   
  filter(!is.na(x_)) %>% filter(!is.na(y_)) |>  
  mutate(year = as.numeric(strftime(t_,format="%Y")),
         month = as.numeric(strftime(t_,format="%m")),
         id.year = paste(individual_id,year,sep=".")) |> 
    group_by(id.year) |> mutate(n_months = n_distinct(month)) |> 
    ungroup() |>   filter(n_months >= min_months) |>
    group_by(id.year)  |> filter(n() >= min_locs_n) |> ungroup() |> 
    dplyr::select(-c(n_months)) |> 
    dplyr::select(x_,y_,id.year) 

  if (nrow(dat.mcp.annual) == 0) return(NULL)
  
  coordinates(dat.mcp.annual) <- c("x_","y_")
  proj4string(dat.mcp.annual) <- crs_4326  
  dat.mcp.annual <- tryCatch(
    sp::spTransform(dat.mcp.annual, crs_moll),
    error = function(e) NULL) 
  if (is.null(dat.mcp.annual)) return(NULL)
  
  mcp_spdf <- tryCatch(
    adehabitatHR::mcp(dat.mcp.annual, percent = 95, unout = c("m2")),
    error = function(e) NULL) 
  if (is.null(mcp_spdf) || length(mcp_spdf@polygons) == 0) return(NULL)
  
  mcp_df <- data.frame(mcp_spdf) |>
    mutate(
      id.year = id,
      year = stringr::str_extract(id, "[^.]*$"),
      individual_id = str_extract(id, "[^.]+")) |>
    filter(!is.na(area))
  
  verts_df <- get_polygon_vertices(mcp_spdf)
  if (is.null(verts_df) || nrow(verts_df) == 0) return(NULL)
  
  mcp_out <- dplyr::left_join(mcp_df, verts_df, by = "id")
  if (nrow(mcp_out) == 0) return(NULL)
  
  vertices_long <- mcp_out |>
    dplyr::select(id, individual_id, year, area, x_vertices, y_vertices) |>
    dplyr::mutate(row_id = dplyr::row_number()) |>
    tidyr::separate_rows(x_vertices, y_vertices, sep = ";") |>
    dplyr::mutate(
      x_vertices = as.numeric(x_vertices),
      y_vertices = as.numeric(y_vertices)
    ) |>
    dplyr::filter(!is.na(x_vertices), !is.na(y_vertices))
  
  if (nrow(vertices_long) == 0) return(NULL)
  
  grid_100 <- dgGEO_to_SEQNUM(dggs_100, vertices_long$x_vertices, vertices_long$y_vertices)$seqnum
  grid_10 <- dgGEO_to_SEQNUM(dggs_10, vertices_long$x_vertices, vertices_long$y_vertices)$seqnum
  grid_1  <- dgGEO_to_SEQNUM(dggs_1,  vertices_long$x_vertices, vertices_long$y_vertices)$seqnum
  
  vertices_long <- vertices_long |>
    dplyr::mutate(
      grid.id.100km = grid_100,
      grid.id.10km = grid_10,
      grid.id.1km = grid_1
    )
  
  grid_summary <- vertices_long |>
    dplyr::group_by(id, individual_id, year, area) |>
    dplyr::summarise(
      x_vertices = paste(x_vertices, collapse = ";"),
      y_vertices = paste(y_vertices, collapse = ";"),
      grid.id.100km = paste(sort(unique(grid.id.100km[!is.na(grid.id.100km)])), collapse = ";"),
      grid.id.10km = paste(sort(unique(grid.id.10km[!is.na(grid.id.10km)])), collapse = ";"),
      grid.id.1km = paste(sort(unique(grid.id.1km[!is.na(grid.id.1km)])), collapse = ";"),
      .groups = "drop"
    )
  
  if (nrow(grid_summary) == 0) return(NULL)
  grid_summary
  
  }, error = function(e) {
    NULL
  })
  
  out
}

## ----function to summarize 12m MCP--------------------
f_sum.ind.mcp12m<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA_character_,n.mcp12m.years = NA,
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
