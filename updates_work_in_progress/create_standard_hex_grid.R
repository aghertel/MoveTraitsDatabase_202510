library(raster)
library(terra)
library(sf)
library(stars)
library(ggplot2)
library(dggridR)

# Use a hexagonal, equal area  grid to spatially annotate traits 
# dggridR builds discrete global grids which partition the surface of the Earth into 
# hexagonal cells, all of which have the same size. 

#Construct a global grid with cells approximately 
# 1 x 1 km = 1 km^2 - res=16
# 10 x 10 km = 25 km^2 -res=13
# 100 x 100 km = 10000 km^2 -res=8

# cell area size 10000 square km
dggs.100          <- dgconstruct(projection = "ISEA",
                             area = 10000, 
                             resround='nearest')
res.100          <- dg_closest_res_to_area(dggs.100,10000)
dggs.100         <- dgsetres(dggs.100,res.100)

# extract cell centers
maxcell.100 <- dgmaxcell(dggs.100) 
cell_ids.100 <- 1:maxcell.100
centroids.100 <- dgSEQNUM_to_GEO(dggs.100, cell_ids.100)

lookup_df.hex.100km <- data.frame(
  running_id = seq_along(cell_ids.100),
  lon = centroids.100$lon_deg,
  lat = centroids.100$lat_deg
)

#Get the grid cell boundaries
grid.100          <- dgearthgrid(dggs.100, savegrid = NA, return_sf = TRUE)

# The grid cells include a running number under "seqnum"
head(grid.100$seqnum)

# save grid and lookup dataframe
saveRDS(grid.100, "StandardGrid/global_10000sqkm_hex.rds")
saveRDS(lookup_df.hex.100km, "StandardGrid/global_10000sqkm_coordinates.rds")


#-------------------------------------------------------------------------------

## cell area size 100 square km
dggs.10          <- dgconstruct(projection = "ISEA",
                                area = 100,
                                metric = FALSE,
                                resround='nearest')
res.10          <- dg_closest_res_to_area(dggs.10,100)
dggs.10         <- dgsetres(dggs.10,res.10)

# extract cell centers
maxcell.10 <- dgmaxcell(dggs.10) 
cell_ids.10 <- 1:maxcell.10
centroids.10 <- dgSEQNUM_to_GEO(dggs.10, cell_ids.10)

lookup_df.hex.10km <- data.frame(
  running_id = seq_along(cell_ids.10),
  lon = centroids.10$lon_deg,
  lat = centroids.10$lat_deg
)

#Get the grid cell boundaries
grid.10          <- dgearthgrid(dggs.10, savegrid = NA, return_sf = TRUE)

saveRDS(grid.10, "StandardGrid/global_100sqkm_hex.rds")
saveRDS(lookup_df.hex.10km, "StandardGrid/global_100sqkm_coordinates.rds")

#-------------------------------------------------------------------------------

## cell area size 1 square km
dggs.1          <- dgconstruct(projection = "ISEA",
                               area = 1,
                               metric = FALSE,
                               resround='nearest')
res.1          <- dg_closest_res_to_area(dggs.1,1)
dggs.1         <- dgsetres(dggs.1,res.1)

# extract cell centers
maxcell.1 <- dgmaxcell(dggs.1) 

# cell_ids.1 <- 1:107616803
# cell_ids.2 <- 107616804:215233606
# cell_ids.3 <- 215233607:322850409
# cell_ids.4 <- 322850410:430467212
# 
# centroids.1 <- dgSEQNUM_to_GEO(dggs.1, cell_ids.4)
# 
# lookup_df.hex.1km.b <- data.frame(
#   running_id = 322850410:430467212,
#   lon = centroids.1$lon_deg,
#   lat = centroids.1$lat_deg
# )
# 
# saveRDS(lookup_df.hex.1km.b, "StandardGrid/global_1sqkm_coordinates_d.rds")

# a<-readRDS("StandardGrid/global_1sqkm_coordinates_a.rds")
# b<-readRDS("StandardGrid/global_1sqkm_coordinates_b.rds")
# c<-readRDS("StandardGrid/global_1sqkm_coordinates_c.rds")
# d<-readRDS("StandardGrid/global_1sqkm_coordinates_d.rds")
# 
# lookup_df.hex.1km <- rbind(a,b,c,d)

saveRDS(lookup_df.hex.1km, "StandardGrid/global_1sqkm_coordinates.rds")

#Get the grid cell boundaries
#grid.1          <- dgearthgrid(dggs.1, savegrid = NA)

#saveRDS(grid.1, "StandardGrid/global_1sqkm_hex.rds")

