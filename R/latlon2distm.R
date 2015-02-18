# Computing a distance matrix in meters among points in decimal degrees
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015


latlon2distm <- function(XY) {
  out <- matrix(,nrow(XY), nrow(XY)) 
  for(i in 1:nrow(XY)) {
    for(j in 1:nrow(XY)) {
      lat1 <- XY[i, 1] 
      lon1 <- XY[i, 2] 
      lat2 <- XY[j, 1] 
      lon2 <- XY[j, 2] 
      R <- 6371                                
      dLat <- (lat2 - lat1) * pi / 180
      dLon <- (lon2 - lon1) * pi / 180
      a <- sin((dLat/2)) ^ 2 + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * (sin(dLon / 2)) ^2
      c <- 2 * atan2(sqrt(a), sqrt(1-a))
      d <- R * c  
      out[i, j] <- d
    }
  }
  rownames(out) <- rownames(XY)
  colnames(out) <- rownames(XY)
  as.dist(out)
}

