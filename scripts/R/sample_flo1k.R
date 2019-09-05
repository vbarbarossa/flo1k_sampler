
# takes xy_coord file in input, reallocate them and appends flo1k values
# the xy_coord file must be in csv format and with xy coords labeled as "long" and "lat" and
# upstream area as "area" (if reallocation is toggled on)

# results are stored in flo1k_sampler/out/
# for implementation of a parallelized version see flo1k_vec

xy_coord <- 'xLivneh_manually_adjusted'
REALLOCATE = TRUE

NC <- 10
year_seq <- 1960:2015
Qnames = c('qav','qma','qmi')
dir_flo1k <- '/vol/milkundata/FLO1K/TIFF/'

library(raster); library(sf); library(foreach)

# REALLOCATE POINTS #################################################################################################

#read upstreamArea layer
ras <- raster('/vol/milkunarc/vbarbarossa/xGlobio-aqua/flo1k_hydrography/upstreamArea.tif')
# set unit cells
unit_cell <- 0.008333333

#load csv with stations coordiantes and area
points <- read.csv(paste0(xy_coord,'.csv'))

if(REALLOCATE){
  
  source('scripts/R/reallocate.R')
  
  # reallocate(1)
  
  results_matrix <- as.data.frame(do.call('rbind',parallel::mclapply(1:nrow(points),reallocate,mc.cores = NC)))
  
  colnames(results_matrix) <- c('new_lon','new_lat','new_area')
  points = cbind(points,results_matrix)
  
  points$area_diff_perc = abs(points$area - points$new_area)*100/points$area
  
  ################################################################################################################
  # QUALITY ASSESSMENT
  total = which(!is.na(points$area_diff_perc)) #8221
  high_quality = which(points$area_diff_perc <= 5 & !is.na(points$area_diff_perc)) #5369
  medium_quality = which(points$area_diff_perc > 5 & points$area_diff_perc <= 10 & !is.na(points$area_diff_perc)) #1206
  low_quality = which(points$area_diff_perc > 10 & points$area_diff_perc <= 50 & !is.na(points$area_diff_perc)) #1646
  
  #number of stations out of the relocation domain
  not_relocated = which(is.na(points$area_diff_perc)) #746
  
  points$quality = NA
  points$quality[high_quality] <- 'high'
  points$quality[medium_quality] <- 'medium'
  points$quality[low_quality] <- 'low'
  points$quality[not_relocated] <- 'not_relocated'
  ###################################################################################################################################
  ###################################################################################################################################
  
}


# SAMPLE FLO1K ######################################################################################################

tab <- points
tab[is.na(tab$area),c('new_lon','new_lat')] <- tab[is.na(tab$area),c('long','lat')]


to.export <- list()
for(var.n in 1:length(Qnames)){
  
  Qname = Qnames[var.n]
  dir.rasters <- paste0(dir_flo1k,Qname,'/')
  ext1year <- function(year,xytab){
    
    e <- round(extract(raster(paste0(dir.rasters,year,'.tif')),xytab),3)
    e[which(e < 0)] <- NA
    
    return(e)
  }
  
  
  ext <- parallel::mcmapply(ext1year,year = year_seq,xytab = rep(list(tab[,c('new_lon','new_lat')]),length(year_seq)),SIMPLIFY = FALSE,mc.cores=NC)
  df <- data.frame(matrix(unlist(ext),nrow=length(ext[[1]])))
  colnames(df) <- year_seq
  
  df$av <- apply(df,1,function(x) round(mean(x,na.rm = TRUE),3))
  df$me <- apply(df,1,function(x) round(median(x,na.rm = TRUE),3))
  df$sd <- apply(df,1,function(x) round(sd(x,na.rm = TRUE),3))
  
  colnames(df) <- paste0(Qname,'_',colnames(df))
  
  
  to.export[[var.n]] <- df
}

write.csv(cbind(tab,to.export[[1]],to.export[[2]],to.export[[3]]),
          paste0(xy_coord,'_w_flo1k.csv'))

