## ----Monthly IOU-------------------------------------------------------------
calc_iou1m <- function(area,
                        trk,
                        dggs_10, 
                        dggs_1, 
                        min_weeks_n = 36) 
{
  tmp.mcp1m <- 

    if(is.null(area)) NULL else {
    area|>
      mutate(id_month = paste(individual_id,year_month,sep="."))
  }

iou1m <- 
  if(is.null(trk) | is.null(tmp.mcp1m)) NULL else {
    trk  %>%  
      mutate(year_month = paste(year,month,sep="_"),
             id_month = paste(individual_id,year_month,sep="."))  |>  group_by(id_month)  |>  
      mutate(cumsum.d24h = sum(d24h,na.rm=T),
             mean.x = mean(x_),
             mean.y = mean(y_))  |>  
      dplyr::select(individual_id,id_month,month,year,year_month,cumsum.d24h,mean.x,mean.y)  |>
      distinct() |> 
      left_join(tmp.mcp1m[,c("id_month","area")],by = "id_month")  |>  
      mutate(iou1m = cumsum.d24h/sqrt(area))  |>  
      filter(!is.na(iou1m)) |> ungroup() |> 
      dplyr::select(individual_id,month,year,year_month,iou1m,mean.x,mean.y) |> 
      mutate(individual_id = as.character(individual_id))
  }

if(is.null(iou1m)) NULL else {
  # Spatial annotation 10km
  cell_info_10 <- dgGEO_to_SEQNUM(dggs.10, iou1m$mean.x, iou1m$mean.y)
  iou1m$grid.id.10km <- cell_info_10$seqnum
  centers_10 <- dgSEQNUM_to_GEO(dggs.10, iou1m$grid.id.10km)
  iou1m$lon.10km <- centers_10$lon_deg
  iou1m$lat.10km <- centers_10$lat_deg
  
  # Spatial annotation 1km
  cell_info_1 <- dgGEO_to_SEQNUM(dggs.1, iou1m$mean.x, iou1m$mean.y)
  iou1m$grid.id.1km <- cell_info_1$seqnum
  centers_1 <- dgSEQNUM_to_GEO(dggs.1, iou1m$grid.id.1km)
  iou1m$lon.1km <- centers_1$lon_deg
  iou1m$lat.1km <- centers_1$lat_deg
}
rm(cell_info_10);rm(centers_10)
rm(cell_info_1);rm(centers_1)

return(iou1m)
}

## ----function to summarize IoU1m--------------------
f_sum.ind.iou1m<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n.iou1m.month = NA,
                       iou1m.mean = NA,iou1m.median = NA,iou1m.cv = NA,iou1m.95 = NA,iou1m.05 = NA)
  } else {
    
    individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
    
    # Get sample size per indivindividual_idual
    n.iou1m.month<-as.numeric(with(x, tapply(x$year_month,individual_id, length)))
    
    # 24hr Displacement
    iou1m.mean<-as.numeric(with(x, tapply(x$iou1m+0.001,individual_id, mean, na.rm=T)))
    iou1m.median<-as.numeric(with(x, tapply(x$iou1m+0.001,individual_id, median, na.rm=T)))
    iou1m.cv<-as.numeric(with(x, tapply(x$iou1m+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
    iou1m.95<-as.numeric(with(x, tapply(x$iou1m+0.001,individual_id, quantile,.95, na.rm=T)))
    iou1m.05<-as.numeric(with(x, tapply(x$iou1m+0.001,individual_id, quantile,.05, na.rm=T)))
    
    # build dataframe
    dats<-data.frame(individual_id,n.iou1m.month,
                     iou1m.mean,iou1m.median,iou1m.cv,iou1m.95,iou1m.05)
    
    return(dats)
  }}

