## ----Daily IOU-------------------------------------------------------------
calc_iou24h <- function(area,
                        trk,
                        dggs_10, 
                        dggs_1, 
                        min_weeks_n = 36) 
{

tmp.mcp24h <- 
  if(is.null(area)) NULL else {
    area|>
      mutate(id_ymd = paste(individual_id,ymd,sep="."))
  }

iou24h <- 
  if(is.null(trk) | is.null(area)) NULL else {
    trk |>  
      group_by(ymd) %>% 
      mutate(cumsumD1h = sum(d1h,na.rm=T),
             mean.x = mean(x_),
             mean.y = mean(y_)) %>% 
      dplyr::select(individual_id,ymd,cumsumD1h,mean.x,mean.y) %>% distinct() %>%
      left_join(tmp.mcp24h[,c("ymd","area")],by = "ymd") %>% 
      mutate(iou24h = cumsumD1h/sqrt(area)) %>% 
      filter(!is.na(iou24h)) |> 
      dplyr::select(individual_id,ymd,iou24h,mean.x, mean.y)
  }

if(is.null(iou24h)) NULL else {
  # Spatial annotation 10km
  cell_info_10 <- dgGEO_to_SEQNUM(dggs.10, iou24h$mean.x, iou24h$mean.y)
  iou24h$grid.id.10km <- cell_info_10$seqnum
  centers_10 <- dgSEQNUM_to_GEO(dggs.10, iou24h$grid.id.10km)
  iou24h$lon.10km <- centers_10$lon_deg
  iou24h$lat.10km <- centers_10$lat_deg

  
  # Spatial annotation 1km
  cell_info_1 <- dgGEO_to_SEQNUM(dggs.1, iou24h$mean.x, iou24h$mean.y)
  iou24h$grid.id.1km <- cell_info_1$seqnum
  centers_1 <- dgSEQNUM_to_GEO(dggs.1, iou24h$grid.id.1km)
  iou24h$lon.1km <- centers_1$lon_deg
  iou24h$lat.1km <- centers_1$lat_deg

}
rm(cell_info_10);rm(centers_1);  rm(cell_info_1);rm(centers_10);rm(tmp.mcp24h)
return(iou24h)
}

## ----function to summarize IoU24h--------------------
f_sum.ind.iou24h<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    # Create a placeholder dataframe with NA values
    dats <- data.frame(individual_id = NA,n.iou24h.days = NA,
                       iou24h.mean = NA,iou24h.median = NA,iou24h.cv = NA,iou24h.95 = NA,iou24h.05 = NA)
  } else {
    
    individual_id <- with(x, tapply(as.character(x$individual_id),individual_id, unique))
    
    # Get sample size per indivindividual_idual
    n.iou24h.days<-as.numeric(with(x, tapply(x$ymd,individual_id, length)))
    
    # 24hr Displacement
    iou24h.mean<-as.numeric(with(x, tapply(x$iou24h+0.001,individual_id, mean, na.rm=T)))
    iou24h.median<-as.numeric(with(x, tapply(x$iou24h+0.001,individual_id, median, na.rm=T)))
    iou24h.cv<-as.numeric(with(x, tapply(x$iou24h+0.001,individual_id, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
    iou24h.95<-as.numeric(with(x, tapply(x$iou24h+0.001,individual_id, quantile,.95, na.rm=T)))
    iou24h.05<-as.numeric(with(x, tapply(x$iou24h+0.001,individual_id, quantile,.05, na.rm=T)))
    
    # build dataframe
    dats<-data.frame(individual_id,n.iou24h.days,
                     iou24h.mean,iou24h.median,iou24h.cv,iou24h.95,iou24h.05)
    
    return(dats)
  }}


## ----summarize iou24h at monthly individual level---------
f_sum.monthly.ind.iou24h<-function(x)
{
  # Check if the input is NULL
  if (is.null(x)) {
    return(data.frame(individual_id = NA,year=NA,month=NA,n.iou24h.days = NA,
                      iou24h.mean = NA,iou24h.median = NA,iou24h.cv = NA,iou24h.95 = NA,iou24h.05 = NA))
  } 
  
  # derive month and year from t_
  ym <- substr(x$ymd, 1, 7)  
  
  # group index: individual x month
  id_ym <- interaction(x$individual_id, ym, drop = TRUE)
  
  individual_id <- tapply(as.character(x$individual_id), id_ym, unique)
  ym_grp        <- tapply(ym, id_ym, unique)
  
  # sample size per individual-month
  n.iou24h.days<-as.numeric(tapply(x$ymd, id_ym, length))
  
  iou24h.mean <- as.numeric(with(x, tapply(x$iou24h+0.001,id_ym, mean, na.rm=T)))
  iou24h.median <- as.numeric(with(x, tapply(x$iou24h+0.001,id_ym, median, na.rm=T)))
  iou24h.cv <- as.numeric(with(x, tapply(x$iou24h+0.001,id_ym, function(x) sd(x, na.rm=T) / mean(x, na.rm=T))))
  iou24h.95 <- as.numeric(with(x, tapply(x$iou24h+0.001,id_ym, quantile,.95, na.rm=T)))
  iou24h.05 <- as.numeric(with(x, tapply(x$iou24h+0.001,id_ym, quantile,.05, na.rm=T)))
  
  year  <- as.numeric(sub(".*\\.(\\d{4})-\\d{2}$", "\\1", unique(id_ym)))
  month <- as.numeric(sub(".*\\.\\d{4}-(\\d{2})$", "\\1", unique(id_ym)))
  
  # build dataframe
  dats<-data.frame(individual_id,year,month,n.iou24h.days,
                   iou24h.mean,iou24h.median,iou24h.cv,iou24h.95,iou24h.05)
  
}


