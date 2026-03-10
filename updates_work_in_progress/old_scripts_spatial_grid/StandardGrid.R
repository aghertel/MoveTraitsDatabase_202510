library(raster)
library(terra)
library(sf)
library(stars)
library(ggplot2)
library(dggridR)


##-- Approach 1.  standard aquare grid

# 'makeRasterTemplate()' creates the template raster with a specific grid topology defined by the 
# resolution and the
# projection. This empty raster Will be used to match all the subsequent data in relation to spatial topology and to rasterize
# the species range data

# 1. Create Mollweide raster directly with precise extent
moll_crs <- "ESRI:54009"  # Mollweide projection code

res <- c(100000000, 100000000)
projString <-  moll_crs
makeRasterTemplate <- function(res, projString){
  template <- st_as_stars(rast(nrows = 1, ncols = 1, xmin = -180, xmax = 180, ymin = -90, ymax = 90,
                               crs = "EPSG:4326", vals=as.integer(1)))
  tempvec <- st_contour(template, breaks=0:1)
  templproj <- st_transform(tempvec, projString)
  templMoll <- rast(ext=ext(templproj), resolution=res, crs=projString)
  return(templMoll)
}

# 1km
OneKM <- makeRasterTemplate(res=c(1000, 1000))
values(OneKM) <- 1:ncell(OneKM)

# 5km
FiveKM <- makeRasterTemplate(res=c(5000, 5000))
values(FiveKM) <- 1:ncell(FiveKM)

# 30km 
ThirtyKM <- makeRasterTemplate(res=c(30000, 30000))
values(ThirtyKM) <- 1:ncell(ThirtyKM)

#writeRaster(OneKM, "StandardGrid/global_1km_mollweide.tif", overwrite=TRUE)
#writeRaster(FiveKM, "StandardGrid/global_5km_mollweide.tif", overwrite=TRUE)
#writeRaster(ThirtyKM, "StandardGrid/global_30km_mollweide.tif", overwrite=TRUE)

# create a lookup dataframe containing the center coordinate of each cell in WGS84
# this dataframe can be used to look up coordinates at 1, 5, 30 km resolution

makeCoordinateLookup <- function(data){
  
  # 2. Get cell centers in Mollweide coordinates
  coords_moll <- crds(data, df = TRUE)
  # 3. Convert to WGS84 using sf
  coords_sf <- st_as_sf(
    coords_moll, 
    coords = c("x", "y"), 
    crs = moll_crs
  ) |>
    st_transform(4326)
  
  # 4. Extract coordinates and create lookup table
  coords_wgs84 <- st_coordinates(coords_sf)
  
  lookup_df <- data.frame(
    cell_id = 1:ncell(data),
    longitude = coords_wgs84[, "X"],
    latitude = coords_wgs84[, "Y"]
  )
  return(lookup_df)
  }

# why are there so many NA values for latitude and longitude?

#lookup_df.1km <- makeCoordinateLookup(data=OneKM)
#lookup_df.5km <- makeCoordinateLookup(data=FiveKM)
lookup_df.30km <- makeCoordinateLookup(data=ThirtyKM)
summary(lookup_df.30km)

#write_csv(lookup_df.1km, "StandardGrid/global_1km_coordinates.csv")
#write_csv(lookup_df.5km, "StandardGrid/global_5km_coordinates.csv")
#write_csv(lookup_df.30km, "StandardGrid/global_30km_coordinates.csv")

ger_sf <- map_data("world") %>%
  filter(region == "Germany") |> 
  st_as_sf(coords = c("long", "lat"), crs = 4326, agr = "constant") %>%
  group_by(group) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") %>%
  st_union() %>%
  st_sf() |> 
  st_transform(moll_crs)

cropped_raster <- st_as_sf(as.polygons(crop(ThirtyKM, ger_sf)))

centroids.30.sf <- lookup_df.30km |> 
  filter(cell_id %in% cropped_raster[["layer"]]) |> 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |> 
  st_transform(moll_crs)


ggplot() +
  geom_sf     (data=ger_sf, color=alpha("black"), fill = NA) +
  geom_sf(data = cropped_raster, fill = NA, color = "black", size = 0.2) +
  geom_sf  (data = centroids.30.sf, col = "red", size = 0.1) +
  coord_sf(crs = moll_crs)+
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

##-- Approach 2.  hex grid

# Alternative - use a hexagonal grid
# dggridR builds discrete global grids which partition the surface of the Earth into 
# hexagonal, triangular, or diamond cells, all of which have the same size. 

#Construct a global grid with cells approximately 
# 1 x 1 km = 1 km^2 - res=16
# 5 x 5 km = 25 km^2 -res=13
# 30 x 30 km = 900 km^2 -res=10

# cell area size 900 square km
dggs.30          <- dgconstruct(projection = "ISEA",
                             area = 900, 
                             resround='nearest')
res.30          <- dg_closest_res_to_area(dggs.30,900)
dggs.30         <- dgsetres(dggs.30,res.30)

# extract cell centers
maxcell.30 <- dgmaxcell(dggs.30) 
cell_ids.30 <- 1:maxcell.30
centroids.30 <- dgSEQNUM_to_GEO(dggs.30, cell_ids.30)

lookup_df.hex.30km <- data.frame(
  running_id = seq_along(cell_ids.30),
  lon = centroids.30$lon_deg,
  lat = centroids.30$lat_deg
)

#Get the grid cell boundaries
grid.30          <- dgearthgrid(dggs.30, savegrid = NA, return_sf = TRUE)

# The grid cells include a running number under "seqnum"
head(grid.30$seqnum)

# Plot everything on a flat map - as an exmaple only for Germany because grid is too fine to 
ger_sf <- map_data("world") %>%
  filter(region == "Germany") |> 
  st_as_sf(coords = c("long", "lat"), crs = 4326, agr = "constant") %>%
  group_by(group) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") %>%
  st_union() %>%
  st_sf() 

cropped_poly <- st_intersection(grid.30, ger_sf)

ggplot() +
  geom_sf     (data=ger_sf, color=alpha("black"), fill = NA) +
  geom_sf     (data=cropped_poly, color=alpha("black"), fill = NA) +
  theme_minimal() +
  geom_point  (data = lookup_df.hex.30km[lookup_df.hex.30km$running_id %in% cropped_poly$seqnum,], 
               aes(x=lon, y=lat), col = "red", size = 0.1) 


# save grid and lookup dataframe
saveRDS(grid.30, "StandardGrid/global_900sqkm_hex.rds")
saveRDS(lookup_df.hex.30km, "StandardGrid/global_900sqkm_coordinates.rds")


# # repeat for finer grids
# # cell area size 1 square km
dggs.1          <- dgconstruct(projection = "ISEA",
                               area = 1,
                               metric = FALSE,
                               resround='nearest')
res.1          <- dg_closest_res_to_area(dggs.1,1)
dggs.1         <- dgsetres(dggs.1,res.1)

# extract cell centers
maxcell.1 <- dgmaxcell(dggs.1) 
cell_ids.1 <- 1:maxcell.1
centroids.1 <- dgSEQNUM_to_GEO(dggs.1, cell_ids.1)

lookup_df.hex.1km <- data.frame(
  running_id = seq_along(cell_ids.1),
  lon = centroids.1$lon_deg,
  lat = centroids.1$lat_deg
)

#Get the grid cell boundaries
grid.1          <- dgearthgrid(dggs.1, savegrid = NA)

saveRDS(grid.1, "StandardGrid/global_1sqkm_hex.rds")
saveRDS(lookup_df.hex.1km, "StandardGrid/global_1sqkm_coordinates.rds")


## cell area size 25 square km
dggs.5          <- dgconstruct(projection = "ISEA",
                               area = 25,
                               metric = FALSE,
                               resround='nearest')
res.5          <- dg_closest_res_to_area(dggs.5,25)
dggs.5         <- dgsetres(dggs.5,res.5)

# extract cell centers
maxcell.5 <- dgmaxcell(dggs.5) 
cell_ids.5 <- 1:maxcell.5
centroids.5 <- dgSEQNUM_to_GEO(dggs.5, cell_ids.5)

lookup_df.hex.5km <- data.frame(
  running_id = seq_along(cell_ids.5),
  lon = centroids.5$lon_deg,
  lat = centroids.5$lat_deg
)

#Get the grid cell boundaries
grid.5          <- dgearthgrid(dggs.5, savegrid = NA, return_sf = TRUE)

saveRDS(grid.5, "StandardGrid/global_25sqkm_hex.rds")
saveRDS(lookup_df.hex.5km, "StandardGrid/global_25sqkm_coordinates.rds")

