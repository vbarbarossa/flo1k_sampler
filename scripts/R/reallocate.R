# Valerio Barbarosa 29/7/2019
# reallocation routine is described in Barbarossa et al. 2018 https://www.nature.com/articles/sdata201852

#define function to transform lon, lat to col,row
lonlat2colrow <- function(lon,lat,unit_cell){
  #translate
  Tlon = lon + 180
  Tlat = 90 - lat
  
  #convert to ncol, nrow
  Ncol = Tlon/unit_cell
  Nrow = Tlat/unit_cell
  
  #center
  Ccol = trunc(Ncol + 1)
  Crow = trunc(Nrow + 1)
  
  return(c(Ccol,Crow))
  
}

#define function to back transform col,row to lon,lat
colrow2lonlat <- function(col,row,unit_cell){
  backT.lon = ((col - 0.5)*unit_cell) - 180
  backT.lat = 90 - ((row - 0.5)*unit_cell)
  
  return(c(backT.lon,backT.lat))
}


###################################################################################################################################
###################################################################################################################################

# core function that scrolls through rows of the table and reallocate the points
# takes in input the lonlat table row and number of cells search radius
reallocate <- function(i,N_ring=5){
  
  # define side dimension of square with "radius" N_ring
  dim_side = N_ring*2+1
  
  # convert latlong to cell locations (row,col)
  convert = lonlat2colrow(points[i,"long"],points[i,"lat"],unit_cell)
  ind_col = convert[1]
  ind_row = convert[2]
  
  # read reported upstream area
  area = points[i,"area"]
  # read matrix of area values surrounding the point location
  mat = matrix(ras[(ind_row-N_ring):(ind_row+N_ring),(ind_col-N_ring):(ind_col+N_ring)],
               nrow=dim_side,ncol=dim_side,byrow = TRUE)
  
  # read areas and calculate euclidean distance at each point location of the surrounding square -----
  ind = N_ring+1
  
  temp_area = numeric(dim_side**2)
  temp_dist = numeric(dim_side**2)
  
  v = -N_ring:N_ring
  
  for(r in 1:dim_side){
    for(c in 1:dim_side) {
      temp_area[c+dim_side*(r-1)] = mat[ind+v[r],ind+v[c]]
      temp_dist[c+dim_side*(r-1)] = sqrt(v[r]**2 + v[c]**2)
      
    }
  }
  neigh = data.frame(Area = temp_area, Dist = temp_dist)
  # --------------------------------------------------------------------------------------------------
  
  # calculate differences and rank
  neigh$diff = abs(area - neigh[,1])
  neigh$diffp = abs(area - neigh[,1])*100/area
  neigh$rank = neigh$diffp + 2*10*neigh$Dist
  neigh$rank[neigh$diffp > 50] <- NA
  mat.assign = matrix(1:dim_side**2,dim_side,dim_side,byrow = TRUE)
  
  # search for best match
  if((sum(neigh$diff,na.rm = TRUE) != 0) & (min(neigh$diffp,na.rm = TRUE) <= 50)){
    
    best.match = which(neigh$rank == min(neigh$rank,na.rm = TRUE))
    best.match = best.match[which(neigh$Dist[best.match] == min(neigh$Dist[best.match]))]
    
    # when solution not unique, choose closest cells first
    if(length(best.match) > 1)  best.match = sample(best.match,1)
    
    new_row = which(mat.assign == best.match,arr.ind = TRUE)[1] + (ind_row - ind)
    new_col = which(mat.assign == best.match,arr.ind = TRUE)[2] + (ind_col - ind)
    
    back_convert = colrow2lonlat(new_col,new_row,unit_cell)
    
    final = c(back_convert,neigh[best.match,1]) #lon,lat,area
    
  }
  
  # report missing value if no match was found
  if((sum(neigh$diff,na.rm = TRUE) == 0) | (min(neigh$diffp,na.rm = TRUE) > 50)) {
    
    final = rep(NA,3)
  }

  return(final)
  
}
