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
OneKM <- makeRasterTemplate(res=c(1000, 1000),projString=moll_crs)
values(OneKM) <- 1:ncell(OneKM)

# 5km
FiveKM <- makeRasterTemplate(res=c(5000, 5000),projString=moll_crs)
values(FiveKM) <- 1:ncell(FiveKM)

# 30km 
ThirtyKM <- makeRasterTemplate(res=c(30000, 30000),projString=moll_crs)
values(ThirtyKM) <- 1:ncell(ThirtyKM)

# 1000Km -- for overview map code!!
FivehunKM <- makeRasterTemplate(res=c(500000, 500000),projString=moll_crs)
values(FivehunKM) <- 1:ncell(FivehunKM)

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
    crs = moll_crs) |>
    st_transform("EPSG:4326")
  
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
#lookup_df.30km <- makeCoordinateLookup(data=ThirtyKM)
lookup_df.500km <- makeCoordinateLookup(data=FivehunKM)

summary(lookup_df.30km)

#write_csv(lookup_df.1km, "StandardGrid/global_1km_coordinates.csv")
#write_csv(lookup_df.5km, "StandardGrid/global_5km_coordinates.csv")
#write_csv(lookup_df.30km, "StandardGrid/global_30km_coordinates.csv")

FivehunKMcheck <- FivehunKM
FivehunKMcheck[] <- sample(1:ncell(FivehunKMcheck),ncell(FivehunKMcheck),replace=T)
FivehunKMcheck[lookup_df.500km$cell_id[is.na(lookup_df.500km$longitude)]] <- NA
#plot(FivehunKMcheck,colNA="red" )

saveRDS(FivehunKMcheck, "StandardGrid/global_500km_mollweide.rds")

 # example map
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


