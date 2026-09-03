## ----Annual IOU-------------------------------------------------------------

calc_iou12m <- function(area, trk) 
{

tmp.mcp12m <- 
  if(is.null(area)) NULL else {
    area|>
      mutate(id_year = paste(individual_id,year,sep="."))
  }

iou12m <- 
  if(is.null(trk) | is.null(tmp.mcp12m)) NULL else {
    trk  |>   
      mutate(id_year = paste(individual_id,year,sep="."))  |>  group_by(id_year) |>  
      mutate(cumsum.d24h = sum(d24h,na.rm=T)) |>  
      dplyr::select(id_year,year,individual_id,cumsum.d24h) |>  
      distinct() |> 
      left_join(tmp.mcp12m[,c("id_year","area", "x_vertices", "y_vertices",
                              "grid.id.100km", "grid.id.10km", "grid.id.1km")],
                by = "id_year") |>  
      mutate(iou12m = cumsum.d24h/sqrt(area)) |> 
      filter(!is.na(iou12m)) |>  ungroup() |> 
    dplyr::select(individual_id,year,iou12m,
                  x_vertices, y_vertices, grid.id.100km, grid.id.10km, grid.id.1km)
  }

if(is.null(iou12m) || nrow(iou12m) == 0) return(NULL)

return(iou12m)
}
## ----function to summarize IoU12m--------------------

f_sum.ind.iou12m<-function(x)
{
  # Check if the input is NULL
  if (is.null(x) || nrow(x) == 0) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA_character_,n.iou12m.year = NA,
                       iou12m.mean = NA,iou12m.median = NA,iou12m.cv = NA,iou12m.95 = NA,iou12m.05 = NA)
  } else {
    
    individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
    
    # Get sample size per indivindividual_idual
    n.iou12m.year<-as.numeric(with(x, tapply(x$year,individual_id, length)))
    
    # 24hr Displacement
    iou12m.mean<-as.numeric(with(x, tapply(x$iou12m+0.001,individual_id, mean, na.rm=T)))
    iou12m.median<-as.numeric(with(x, tapply(x$iou12m+0.001,individual_id, median, na.rm=T)))
    iou12m.cv<-as.numeric(with(x, tapply(x$iou12m+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
    iou12m.95<-as.numeric(with(x, tapply(x$iou12m+0.001,individual_id, quantile,.95, na.rm=T)))
    iou12m.05<-as.numeric(with(x, tapply(x$iou12m+0.001,individual_id, quantile,.05, na.rm=T)))
    
    # build dataframe
    dats<-data.frame(individual_id,n.iou12m.year,
                     iou12m.mean,iou12m.median,iou12m.cv,iou12m.95,iou12m.05)
    
    return(dats)
  }}

