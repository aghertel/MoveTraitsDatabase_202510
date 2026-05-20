library(dplyr)
library(tidyr)

## ----Species level Database-------------------------------------------------------------
MoveTrait <- readRDS("./updates_work_in_progress/DATA/8.MoveTraits_db/MoveTrait.v0.1_withinindividual_20260519.rds")

dmax24h.sp.mean <- 
  MoveTrait |> 
  unnest(dmax24h) |> 
  dplyr::select("individual_id","species","common_name","dmax24h") |> 
  group_by(species) |> 
  mutate(dmax24h.mean = mean(dmax24h),
         dmax24h.sd = sd(dmax24h),
         n = tally(individual_id)) |> 
  distinct()








## species summaries
# Function to calculate the coefficient of variation
cv <- function(x, na.rm = TRUE) {
  # Calculate standard deviation and mean
  sd_value <- sd(x, na.rm = na.rm)
  mean_value <- mean(x, na.rm = na.rm)
  cv_value <- sd_value / abs(mean_value)
  return(cv_value)
}

MoveTrait.v0.1.sp <- MoveTrait.v0.1 |> 
  mutate(common_name = recode(common_name, "reindeer" = "reindeer/caribou")) |> 
  mutate(common_name = recode(common_name, "elk" = "red deer/elk")) |> 
  mutate(common_name = recode(common_name, "red deer" = "red deer/elk"))

MoveTrait.v0.1.sp <-
  MoveTrait.v0.1.sp |> 
  dplyr::select(3:6,11:90,95) 

MoveTrait.v0.1.sp2 <- 
  MoveTrait.v0.1.sp |> 
  group_by(species) |>
  mutate(species = unique(species),
         common_name = unique(common_name),
         class = unique(class),
         movement.mode = unique(movement.mode),
         grid.id.100km = paste(unique(trimws(unlist(strsplit(na.omit(grid.id.100km), ";")))),collapse = ";"),
         grid.id.10km = paste(unique(trimws(unlist(strsplit(na.omit(grid.id.100km), ";")))),collapse = ";"),
         across(c("n1h","n24h.days","n.dmax24h.days","n.dmax1m.weeks",
                  "n.max12m.years","n.mcp24h.days","n.mcp1m.months",
                  "n.mcp12m.years","n.iou24h.days","n.iou1m.month", "n.iou12m.year", 
                  "n.di.days"), 
                sum, na.rm = TRUE),
         across(c("d1h.mean",       "d1h.median",     "d1h.cv",         "d1h.95",        
                  "d1h.05",         "d24h.mean",      "d24h.median",    "d24h.cv",       
                  "d24h.95",        "d24h.05",        "dmax24h.mean",   "dmax24h.median",
                  "dmax24h.cv",     "dmax24h.95",     "dmax24h.05",     "dmax1m.mean",   
                  "dmax1m.median",  "dmax1m.cv",      "dmax1m.95",      "dmax1m.05",     
                  "dmax12m.mean",   "dmax12m.median", "dmax12m.cv",     "dmax12m.95",    
                  "dmax12m.05",     "mcp24h.mean",    "mcp24h.median",  "mcp24h.cv",     
                  "mcp24h.95",      "mcp24h.05",      "mcp1m.mean",    
                  "mcp1m.median",   "mcp1m.cv",       "mcp1m.95",       "mcp1m.05",      
                  "mcp12m.mean",    "mcp12m.median",  "mcp12m.cv",      "mcp12m.95",     
                  "mcp12m.05",      "iou24h.mean",    "iou24h.median",  "iou24h.cv",     
                  "iou24h.95",      "iou24h.05",      "iou1m.mean",     "iou1m.median",  
                  "iou1m.cv",       "iou1m.95",       "iou1m.05",       "iou12m.mean",   
                  "iou12m.median",  "iou12m.cv",      "iou12m.95",      "iou12m.05",     
                  "di.mean",        "di.median",      "di.cv",          "di.95",         
                  "di.05"), 
                mean, na.rm = TRUE),
         contact_person_name = paste(unique(contact_person_name), collapse = ", ")) |> 
  distinct()


saveRDS(MoveTrait.v0.1.sp2, file=paste0(pthdb,"MoveTrait.v0.1_species.sum_20251011.rds"))
