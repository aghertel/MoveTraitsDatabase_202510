## ----Annual IOU-------------------------------------------------------------

calc_iou12m <- function(area,
                        trk,
                        dggs_10, 
                        dggs_1, 
                        min_weeks_n = 36) 
{
 area <-mcp12m
 trk <-d24h

tmp.mcp12m <- 
  if(is.null(area)) NULL else {
    area|>
      mutate(id_year = paste(individual_id,year,sep="."))
  }

iou12m <- 
  if(is.null(trk) | is.null(tmp.mcp12m)) NULL else {
    trk  |>   
      mutate(id_year = paste(individual_id,year,sep="."))  |>  group_by(id_year) |>  
      mutate(cumsum.d24h = sum(d24h,na.rm=T),
             mean.x = mean(x_),
             mean.y = mean(y_)) |>  
      dplyr::select(id_year,year,individual_id,cumsum.d24h,mean.x,mean.y) |>  distinct() |> 
      left_join(tmp.mcp12m[,c("id.year","area")],by = c("id_year"="id.year")) |>  
      mutate(iou12m = cumsum.d24h/sqrt(area)) |> 
      filter(!is.na(iou12m)) |>  ungroup() |> 
    dplyr::select(individual_id,year,iou12m,mean.x,mean.y)
  }

if(is.null(iou12m)) NULL else {
  # Spatial annotation 10km
  cell_info_10 <- dgGEO_to_SEQNUM(dggs.10, iou12m$mean.x, iou12m$mean.y)
  iou12m$grid.id.10km <- cell_info_10$seqnum
  centers_10 <- dgSEQNUM_to_GEO(dggs.10, iou12m$grid.id.10km)
  iou12m$lon.10km <- centers_10$lon_deg
  iou12m$lat.10km <- centers_10$lat_deg

  # Spatial annotation 1km
  cell_info_1 <- dgGEO_to_SEQNUM(dggs.1, iou12m$mean.x, iou12m$mean.y)
  iou12m$grid.id.1km <- cell_info_1$seqnum
  centers_1 <- dgSEQNUM_to_GEO(dggs.1, iou12m$grid.id.1km)
  iou12m$lon.1km <- centers_1$lon_deg
  iou12m$lat.1km <- centers_1$lat_deg

}
rm(cell_info_1);rm(centers_1)  
rm(cell_info_10);rm(centers_10)

return(iou12m)
}
## ----function to summarize IoU12m--------------------

f_sum.ind.iou12m<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n.iou12m.year = NA,
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

