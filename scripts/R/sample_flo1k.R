# Valerio Barbarosa 29/7/2019
# takes xy_coord file in input, reallocate them and appends flo1k values
# the xy_coord file must be in csv format and with xy coords labeled as "long" and "lat" and
# upstream area as "area" (if reallocation is toggled on)

# results are stored in flo1k_sampler/out/

# USER SETTINGS -------------------------------------------------------------------------------------------------
# make sure the following libraries are available
library(raster); library(sf); library(foreach)

# input file (without the extension)
xy_coord <- 'infile'

# reallocate points based on upstream area? T/F
REALLOCATE = TRUE

# number of cores for parallelized execution (useful for > ~1,000 points [not benchmarked])
NC <- 10

# time span (in years) and variables needed
year_seq <- 1960:2015
Qnames = c('qav','qma','qmi')

# directory with FLO1K layers in GeoTIFF format (one layer per year-variable)
# these can be exported from the netCDF files available with the publication Barbarossa et al. 2018 https://www.nature.com/articles/sdata201852
# or can be kindly requested via email to v.barbarossa@fnwi.ru.nl
dir_flo1k <- '/path_to_FLO1K_TIFF_directory/'

# this can be derived from the global hydrography (see Barbarossa et al. 2018)
# or can be kindly requested via email to v.barbarossa@fnwi.ru.nl
ras <- raster('/path_to_directory/upstreamArea.tif')

# set unit cells [not tested for different resolutions]
unit_cell <- 0.008333333
#---------------------------------------------------------------------------------------------------------------


# REALLOCATE POINTS #################################################################################################


#load csv with stations coordiantes and area
points <- read.csv(paste0(xy_coord,'.csv'))

if(REALLOCATE){
  
  source('scripts/R/reallocate.R')
  
  results_matrix <- as.data.frame(do.call('rbind',parallel::mclapply(1:nrow(points),reallocate,mc.cores = NC)))
  
  colnames(results_matrix) <- c('new_lon','new_lat','new_area')
  points = cbind(points,results_matrix)
  
  points$area_diff_perc = abs(points$area - points$new_area)*100/points$area
  
  ################################################################################################################
  # QUALITY ASSESSMENT
  total = which(!is.na(points$area_diff_perc))
  high_quality = which(points$area_diff_perc <= 5 & !is.na(points$area_diff_perc))
  medium_quality = which(points$area_diff_perc > 5 & points$area_diff_perc <= 10 & !is.na(points$area_diff_perc))
  low_quality = which(points$area_diff_perc > 10 & points$area_diff_perc <= 50 & !is.na(points$area_diff_perc))
  
  #number of stations out of the reallocation domain
  not_relocated = which(is.na(points$area_diff_perc))
  
  points$quality = NA
  points$quality[high_quality] <- 'high'
  points$quality[medium_quality] <- 'medium'
  points$quality[low_quality] <- 'low'
  points$quality[not_relocated] <- 'not_relocated'
  
}


# SAMPLE FLO1K ######################################################################################################

tab <- points
tab[is.na(tab$area),c('new_lon','new_lat')] <- tab[is.na(tab$area),c('long','lat')]


to.export <- list()
for(var.n in 1:length(Qnames)){
  
  # read variables name
  Qname = Qnames[var.n]
  
  # raster layers directory
  dir.rasters <- paste0(dir_flo1k,Qname,'/')
  
  #function that extract one year-variable value for all the points
  ext1year <- function(year,xytab){
    
    e <- round(extract(raster(paste0(dir.rasters,year,'.tif')),xytab),3)
    e[which(e < 0)] <- NA
    
    return(e)
  }
  
  # run the extraction function in parallel and collapse results in a data frame
  ext <- parallel::mcmapply(ext1year,year = year_seq,xytab = rep(list(tab[,c('new_lon','new_lat')]),length(year_seq)),SIMPLIFY = FALSE,mc.cores=NC)
  df <- data.frame(matrix(unlist(ext),nrow=length(ext[[1]])))
  colnames(df) <- year_seq
  
  # calculate additional statistics (average, median, standard deviation)
  df$av <- apply(df,1,function(x) round(mean(x,na.rm = TRUE),3))
  df$me <- apply(df,1,function(x) round(median(x,na.rm = TRUE),3))
  df$sd <- apply(df,1,function(x) round(sd(x,na.rm = TRUE),3))
  
  colnames(df) <- paste0(Qname,'_',colnames(df))
  
  # return one data frame per variable
  to.export[[var.n]] <- df
}

# create out directory if it does not exists
if(!dir.exists('out')) dir.create('out')

# write the output in CSV format
write.csv(cbind(tab,do.call('cbind',to.export)),
          paste0('out/',xy_coord,'_w_flo1k.csv'))

