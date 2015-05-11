# Spatial weights

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015

setGeneric("eco.weight", function(XY,
                                  method = c("circle", "knearest", "inverse",  
                                             "circle.inverse", "exponential", 
                                             "circle.exponential"),
                                  d1 = 0, d2 = NULL,  k = NULL,  p = 1, alpha = 1, 
                                  dist.method = "euclidean",
                                  row.sd = FALSE, self = FALSE, latlon = FALSE) {
  
  method <- match.arg(method)
  dist.method <- match.arg(dist.method)
  
  #distance configuration
  if(latlon == FALSE) {
    distancia <- as.matrix(dist(XY), upper = T, method = dist.method)
  } else {
    distancia <- dist(SoDA::geoXY(XY[,2], XY[,1], unit=1))
    distancia <- as.matrix(distancia, upper = T, method = dist.method)
  }
  #####
  
  if(method == "inverse") {
    y <- as.matrix(dist(XY, upper = T, method = dist.method))
    y <- 1 / (y ^ p)
    diag(y) <- 0
    
    param <- "p"
    param.values <- p
    
  } else if(method == "circle") {
    
    if(is.null(d2)) {
      stop("a d2 argument must be given")
    }
    temp <- which((distancia <= d2) & (distancia > d1))
    y <- distancia
    y <- y - distancia
    y[temp] <- 1
    if(self) {
      diag(y) <- 1
    }
    
    param <- c("d1", "d2")
    param.values <- c(d1, d2)
    
  } else if(method == "circle.inverse") {
    
    temp <- which((distancia <= d2) & (distancia >= d1))
    dummy <- distancia - distancia
    dummy[temp] <- 1
    y <- 1 / (distancia ^ p)
    y <- dummy * y
    diag(y) <- 0
    
    param <- c("d1", "d2", "p")
    param.values <- c(d1, d2, p)
    
  } else if(method == "exponential") {
    y <- 1 / exp(alpha * distancia)
    diag(y) <- 0
    
    param <- "alpha"
    param.values <- alpha
    
  } else if(method == "circle.exponential") {
    
    if(is.null(d2)) {
      stop("a d2 argument must be given")
    }
    
    temp <- which((distancia <= d2) & (distancia > d1))
    dummy <- distancia - distancia
    dummy[temp] <- 1
    y <- 1 / exp(alpha * distancia)
    y <- dummy * y
    diag(y) <- 0
    
    param <- c("d1", "d2", "alpha")
    param.values <- c(d1, d2, alpha)
    
  } else if(method == "knearest") {
    
    if(is.null(k)) {
      stop("a k argument must be given")
    }
    
    y <- t(apply(distancia, 1, function(x) {
      x[x != 0] <- rank(x[x != 0], ties.method = "random")
      return(x)}))
    
    y <- t(apply(y, 1, function(x) as.numeric(x <= k & x != 0)))
    if(self) {
      diag(y) <- 1
    }
    
    param <- "k"
    param.values <- k
    
  } 
  
  #row standardization
  if(row.sd) {
    y <- y/apply(y, 1, sum)
    y[is.na(y)] <- 0
  } 
  
  #output construction
  out <- new("eco.weight")
  out@METHOD <- method
  out@ROW.SD <- row.sd
  if(method == "circle" | method =="knearest") {
    out@SELF <- self
  } else {
    out@SELF <- NA
  }
  out@W <- y
  out@XY <- XY
  
  y2 <- y
  diag(y2) <- 0
  out@NONZERO <- round(100* sum(y2 != 0) / (nrow(y2)^2 - nrow(y2)), 1)
  out@NONZEROIND <- round(100 * sum(apply(y2, 1, sum) != 0) / nrow(y2), 1)
  out@AVG <- round(sum(apply(y2, 1, sum))  / nrow(y2), 1)
  out@PAR <- param
  out@PAR.VAL <- param.values
  
  out
  
})
