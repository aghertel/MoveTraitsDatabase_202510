library(tidyverse)
library(dggridR)
library(collapse)
library(rnaturalearth)
library(sf)
library(ggspatial)

## 1. Map Species and Individual locations based on 1hr data
wind.1 <- readRDS("OSF/OUTPUT/MoveTrait.v0.1_withinindividual_20251011.rds")

wind <- wind.1 |>
  dplyr::select(species, individual_id,mean.longitude, mean.latitude, class, displ.1h)

# 5687 individuals
# 243 studies
# 152 species
d1 <- wind |> 
  unnest(displ.1h) |>
  dplyr::select("individual_id","species","class","x_","y_") |>
  rename(x = x_)|>
  rename(y = y_)

length(levels(factor(d1$species)))


####
# # Loxodonta africana
# d1 <- wind.1 |> 
#   dplyr::select("individual_id","species") |>
#   unnest(mcp.1m) |>
#   dplyr::select("individual_id","species","mcp.1m","month") |>
#   rename(x = x_)|>
#   rename(y = y_)



# cell area size 10000 square km
dggs          <- dgconstruct(projection = "ISEA",
                                area = 10000, 
                                resround='nearest')
# area: 7,774.20548	
# spacing of center nodes: min=82.31100,	max=104.47000, mean=95.26360	
res.100          <- dg_closest_res_to_area(dggs,10000)
dggs         <- dgsetres(dggs,res.100)

#Get the corresponding grid cells
d1$cell <- dgGEO_to_SEQNUM(dggs, d1$x, d1$y)$seqnum

## 3. Species map
#Get the number of unique SPECIES
spcounts   <- d1 |> 
  distinct(species, cell, .keep_all = TRUE) |> 
  collapse::fcount(cell)

#Get the grid cell boundaries for cells which had tracking data
grid          <- dgcellstogrid(dggs,spcounts$cell)
grid          <- merge(grid,spcounts,by.x="seqnum",by.y="cell")

wrapped_grid = sf::st_wrap_dateline(grid, options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"), quiet = TRUE)

# Define custom colors using the viridis palette
custom_colors <- viridis::viridis(n = 5, option = "D")

# Manually bin the data into categories
wrapped_grid$bin <- cut(
  wrapped_grid$N,
  breaks = c(-Inf, 1, 3, 6, 12, Inf),
  labels = c("1", "2 - 3", "4 - 6", "7 - 12", "> 12"))

world_map <- ne_countries(scale = "medium", returnclass = "sf")
bbox <- sf::st_bbox(world_map)

map_species <-
ggplot() +
  geom_sf(data = world_map, fill = "gray90", color = "darkgray") +  # Base map
  geom_sf(data=wrapped_grid, aes(fill=bin, alpha = 0.5), color=alpha("white", 0.4)) +
  coord_sf(xlim = c(bbox["xmin"]+15, bbox["xmax"]-15), ylim = c(bbox["ymin"]+7, bbox["ymax"]-5)) +  # Constrain extent
  scale_fill_manual(
    name = "# of species",## of species
    values = custom_colors, 
    labels = c("1", "2 - 3", "4 - 6", "7 - 12", ">12"))+
  labs(title=expression("MoveTraits species traits"[italic(" (effective October 2025)")] )) +
  xlim(expand=c(0,0)) + ylim(expand = c(0,0)) +
  theme_minimal() +
  scale_alpha(guide = 'none')+
  theme(panel.background = element_rect(fill = "#aad3df"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.justification = "left",
        text = element_text(family = "Avenir", size = 25),
        legend.key.size = unit(0.4, "cm")) 

ggsave(
  filename = "traits_map_species_100.png",  
  plot = map_species,             
  width = 24, #1920               
  height = 15,#1200                
  units = "cm", #px              
  dpi = 600                  
)

## 4. One species map
#Get the number of unique INDIVIDUALS in each equally-sized cell
reindeercounts   <- d1 |> 
  filter(species == "Rangifer tarandus") |> 
  distinct(species, cell, .keep_all = TRUE) 

#Get the grid cell boundaries for cells which had tracking data
grid2          <- dgcellstogrid(dggs,reindeercounts$cell)
grid2          <- merge(grid2,reindeercounts,by.x="seqnum",by.y="cell")

wrapped_grid2 = sf::st_wrap_dateline(grid2, options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"), quiet = TRUE)

# 2306, 2387
map_rangifertarandus_100 <- 
ggplot() +
  geom_sf(data = world_map, fill = "gray90", color = "gray90") +  # Base map
  geom_sf     (data=grid2, aes(fill=species, alpha = 0.5), color=alpha("white", 0.4)) +
  geom_sf     (data=grid2[grid2$seqnum %in% c(2306, 2387),], 
               color="red", fill = NA, linewidth = 0.8) +
  coord_sf(crs = "ESRI:102016", 
           xlim = c(-3500000, 3000000), 
           ylim = c(-3500000, 3000000))+
  ggtitle("Species summary",
          subtitle = expression(italic("Rangifer tarandus")))+
  theme_bw() +
  scale_alpha(guide = 'none')+
  theme(panel.background = element_rect(fill = "#aad3df"),#
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.justification = "left",
        text = element_text(family = "Avenir", size = 14),
        legend.key.size = unit(0.4, "cm"),
        legend.title=element_blank(),
        axis.text=element_blank()) +
  scale_fill_manual(labels = expression ("10000"~km^2),values = c("#D55E00"))
#  scale_fill_manual(labels = "10000 km",values = c("#D55E00"))

ggsave(
  filename = "map_rangifertarandus_100.png",   # File name and format
  plot = map_rangifertarandus_100,             # ggplot object to save
  width = 1500, #9
  height = 1500,
  units = "px", #cm
  dpi = 300
)

## 4. One individual map
#Get the number of unique INDIVIDUALS in each equally-sized cell
one_reindeer_points   <- d1 |> 
  filter(species == "Rangifer tarandus") |> 
  filter(individual_id == "BWCA18601_Rangifer tarandus")  
one_reindeer_points_sf <- st_as_sf(one_reindeer_points,                         
               coords = c("x", "y"),
               crs = 4326)
one_reindeer_points_sf <- st_transform(one_reindeer_points_sf, crs = "ESRI:102016")

one_reindeer   <- d1 |> 
  filter(species == "Rangifer tarandus") |> 
  filter(individual_id == "BWCA18601_Rangifer tarandus") |> 
  distinct(individual_id, cell, .keep_all = TRUE) 

#Get the grid cell boundaries for cells which had tracking data
grid2          <- dgcellstogrid(dggs,one_reindeer$cell)
grid2          <- merge(grid2,one_reindeer,by.x="seqnum",by.y="cell")

# cell area size 100 square km
dggs.10          <- dgconstruct(projection = "ISEA",
                                area = 100, 
                                resround='nearest')
res.10          <- dg_closest_res_to_area(dggs.10,100)
dggs.10         <- dgsetres(dggs.10,res.10)

#Get the corresponding grid cells
one_reindeer_points$cell.10 <- dgGEO_to_SEQNUM(dggs.10, one_reindeer_points$x, one_reindeer_points$y)$seqnum
one_reindeer_hex   <- one_reindeer_points |> 
  distinct(individual_id, cell.10, .keep_all = TRUE) 
grid10          <- dgcellstogrid(dggs.10,one_reindeer_hex$cell.10)
grid10          <- merge(grid10,one_reindeer_hex,by.x="seqnum",by.y="cell.10")


x_center <- -2762964
y_center <- 1354626
window_size <- 100000  # adjust as needed for your zoom level


map_one.rangifertarandus_100_10 <- 
  ggplot() +
  geom_sf     (data=grid2, aes(fill=species, alpha = 0.3), color="red") + #alpha("white", 0.4)
    geom_sf     (data=one_reindeer_points_sf, size = 0.3, color = "darkgray") +
    geom_sf     (data=grid10, fill=NA, aes(color="#0072B2", 1.5)) +
  geom_sf     (data=grid10[grid10$seqnum == 188407,], fill=NA, color="blue",linewidth=1) +
  coord_sf(crs = "ESRI:102016",
           xlim = c(x_center - window_size, x_center + window_size),
           ylim = c(y_center - window_size, y_center + window_size))+
    ggtitle("Individual summary",
            subtitle = expression("One " * italic("Rangifer tarandus") ))+
    theme_bw() +
  scale_alpha(guide = 'none')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.justification = "left",
        text = element_text(family = "Avenir", size = 14),
        legend.key.size = unit(0.4, "cm"),
        legend.title=element_blank(),
        axis.title = element_blank(),
        axis.text=element_blank()) +
  annotation_scale(
    location = "bl",
    width_hint = 0.5)+
    scale_fill_manual(labels = expression ("10000"~km^2),values = c("#D55E00")) +
    scale_color_manual(labels = expression ("100"~km^2),values = c("#0072B2"))

ggsave(
  filename = "map_one _rangifertarandus_100_10.png",   # File name and format
  plot = map_one.rangifertarandus_100_10,             # ggplot object to save
  width = 1500, #9
  height = 1500,
  units = "px", #cm
  dpi = 300
)

## 4. One reinder, one day map
d1.reindeer.ind <- wind |> 
  unnest(displ.1h) |>
  dplyr::select("individual_id","species","class","x_","y_","t_") |>
  filter(species == "Rangifer tarandus") |> 
  filter(individual_id == "BWCA18601_Rangifer tarandus") |> 
  rename(x = x_)|>
  rename(y = y_)

d1.reindeer.ind <- d1.reindeer.ind[c(1,25),]

d1.reindeer.ind_sf <- st_as_sf(d1.reindeer.ind,                         
                                   coords = c("x", "y"),
                                   crs = 4326)
d1.reindeer.ind_sf <- st_transform(d1.reindeer.ind_sf, crs = "ESRI:102016")

#Get the corresponding grid cells at 100km^2
d1.reindeer.ind$cell.10 <- dgGEO_to_SEQNUM(dggs.10, d1.reindeer.ind$x, d1.reindeer.ind$y)$seqnum
one_reindeer_day_hex   <- d1.reindeer.ind |> 
  distinct(individual_id, cell.10, .keep_all = TRUE) 
grid10          <- dgcellstogrid(dggs.10,one_reindeer_day_hex$cell.10)
grid10          <- merge(grid10,one_reindeer_day_hex,by.x="seqnum",by.y="cell.10")

# cell area size 1 square km
dggs.1          <- dgconstruct(projection = "ISEA",
                                area = 1, 
                                resround='nearest')
res.1          <- dg_closest_res_to_area(dggs.1,1)
dggs.1         <- dgsetres(dggs.1,res.1)

#Get the corresponding grid cells
d1.reindeer.ind$cell.1 <- dgGEO_to_SEQNUM(dggs.1, d1.reindeer.ind$x, d1.reindeer.ind$y)$seqnum
one_reindeer_day_hex.1   <- d1.reindeer.ind |> 
  distinct(individual_id, cell.1, .keep_all = TRUE) 
grid1          <- dgcellstogrid(dggs.1,one_reindeer_day_hex.1$cell.1)
grid1          <- merge(grid1,one_reindeer_day_hex.1,by.x="seqnum",by.y="cell.1")

x_center <- -2783855
y_center <- 1369000
window_size <- 8000  # adjust as needed for your zoom level

map_one.rangifertarandus_one.day_10_1 <- 
  ggplot() +
  geom_sf     (data=grid10, aes(fill=species, alpha = 0.3), color="blue") + #alpha("white", 0.4)
  geom_sf     (data=grid1, fill=NA, aes(color="#CC79A7", 1.5),linewidth=1.5) +
  geom_sf     (data=d1.reindeer.ind_sf, size = 1.5, aes(shape = individual_id)) +
  coord_sf(crs = "ESRI:102016",
           xlim = c(x_center - window_size, x_center + window_size),
           ylim = c(y_center - window_size, y_center + window_size))+
    ggtitle("Within-individual",
            subtitle = expression("One " * italic("Rangifer tarandus") * ", one day"))+
  theme_bw() +
  scale_alpha(guide = 'none')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.justification = "left",
        text = element_text(family = "Avenir", size = 14),
        legend.key.size = unit(0.4, "cm"),
        legend.title=element_blank(),
        axis.title = element_blank(),
        legend.key=element_blank(),
        axis.text=element_blank()) +
  annotation_scale(
    location = "bl",
    width_hint = 0.6)+
  scale_color_manual(labels = expression ("1"~km^2),values = c("#CC79A7")) +
  scale_shape_manual(labels = c("Observation scale"),values = c(16)) +
scale_fill_manual(labels = expression ("100"~km^2),values = c("#0072B2")) 
  
    
ggsave(
  filename = "map_one.rangifer_one.day_10_1.png",   # File name and format
  plot = map_one.rangifertarandus_one.day_10_1,             # ggplot object to save
  width = 1500, #9
  height = 1500,
  units = "px", #cm
  dpi = 300
)
  
