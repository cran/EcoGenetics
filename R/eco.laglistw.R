# Obtention of a list with weight matrices for lag distance classes
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015

setGeneric("eco.laglistw", function(xy, int, smax, w =c("W", "B")) {
  
  w <- match.arg(w)
  
  distancia <- as.matrix(dist(xy, upper = T))
  
  laglw <- list()
  j <- 1
  for (i in seq(int, smax, int)) {
    temp <- which((distancia <= i) & (distancia > i - int))
    dummy <- distancia
    dummy <- dummy - distancia
    dummy[temp] <- 1
    laglw[[j]] <- dummy
    j <- j+1 
  }
  
  
  if(w == "W") {
    laglist <- list()
    laglist <- lapply(laglw, function(y) y/apply(y, 1, sum))
    for(i in 1:length(laglist)) {
      laglist[[i]][is.na(laglist[[i]])] <- 0
    }
  } else {
    laglist <-laglw
  }
  
  laglist
})
