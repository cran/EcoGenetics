# Leandro Roser leandroroser@ege.fcen.uba.ar
# June 17, 2015 


# Obtention of a list of spatial weights for classes


setGeneric("eco.lagweight", 
          function(XY, 
                   int = NULL, 
          				 smin = 0,
          				 smax = NULL,
          				 kmax = NULL,
          				 nclass = NULL,
          				 seqvec = NULL,
          				 size = NULL,
          				 bin = c("sturges", "FD"),
          				 cummulative = FALSE,
					 				 row.sd = FALSE,
					 				 self = FALSE,
					 				 latlon = FALSE) {
  
  bin <- match.arg(bin)				 	

  
	#variables definitions
  if(ncol(XY) != 2) {
    message("more than 3 columns in coordinates. The first two will we
								taken as X-Y data for estimating distance intervals")
    XY <- XY[,1:2]
  }
  
	  if(latlon == FALSE) {
	    distancia <- dist(XY)
	  } else {
	    distancia <- dist(SoDA::geoXY(XY[,2], XY[,1], unit=1))
	  }
  	distancia <- as.matrix(distancia)
  	logdistancia <- log(distancia)

#method control
match.control <- sum(!is.null(smax), !is.null(kmax), !is.null(seqvec))

if(match.control == 0) {
	smax <- max(distancia)
} else if (match.control != 1) {
	stop("Only one of smax, kmax, nclass or seqvec should be given")
}


	meandist <- vector()
  logdist <- vector()
	cardinal <- vector()
	laglw <- list()
	j <- 1
	
	
	#-----------computation of lag matrices
	
	###computation based in distance
	
	#based on distance between individuals, different size
	if(is.null(kmax) & is.null(size)) {
	  
	  
	  input <- int.break(XY = XY, 
	                     int = int, 
	                     smin = smin,
	                     smax = smax,
	                     nclass = nclass,
	                     seqvec = seqvec,
	                     latlon = latlon,
	                     bin = bin)
	  breakpoints <- input$breakpoints
	  method <- input$method
	  
	  #  intervals (a, b]
	  
	  for(i in 2:length(breakpoints)) {
	    #cummulative distance
	    if(cummulative) {
	      temp <- which(distancia <= breakpoints[i] & (distancia > smin))
	    } else {
	      temp <- which((distancia <= breakpoints[i]) & (distancia > breakpoints[i-1]))
	    }
	    
	    distemp <- distancia[temp]
	    logdistemp <- logdistancia[temp]
	    meandist[j] <- mean(distemp)
	    logdist[j] <- mean(logdistemp)
	    dummy <- distancia
	    dummy <- dummy - distancia
	    dummy[temp] <- 1
	    laglw[[j]] <- dummy
	    cardinal[j] <- sum(dummy) / 2
	    j <- j+1
	  }
	  
	  if(self) {
	    dummy.self <- distancia - distancia
	    diag(dummy.self) <- 1
	    cardinal <- c(length(diag(dummy.self)), cardinal)
	    dummy.self <- list(dummy.self)
	    laglw <- append(dummy.self, laglw)
	    meandist <- c(0, meandist)
	    
	  }
	} 
	
	#based on distance between individuals, equal size
	if(is.null(kmax) & !is.null(size)) {
	  size2 <- 2 * size
	  method <- "equal.size"
	  vec.mat <- as.vector(distancia)
	  largo <- length(vec.mat)
	  names(vec.mat) <- 1:largo
	  vec.mat.nodiag <- vec.mat[vec.mat != 0]
	  vec.sort <- sort(vec.mat.nodiag)
	  cortes <- seq(size2, length(vec.mat.nodiag), size2)
	  cortes <- c(0, cortes)
	  breakpoints <- 0
	  
	  for(i in 2:length(cortes)) {
	    #cummulative distance
	    if(cummulative) {
	      temp <- as.numeric(names(vec.sort[1:cortes[i]]))
	    } else {
	      temp <- as.numeric(names(vec.sort[(cortes[i - 1] +1):cortes[i]]))
	    }
	    
	    distemp <- distancia[temp]
	    logdistemp <- logdistancia[temp]
	    control.corte <- max(distemp)
	    
	    #control
	    if(control.corte > smax) {
	      break
	    }
	    
	    if(i >2) {
	      temp0 <- as.numeric(names(vec.sort[(cortes[i - 2] +1):cortes[i-1]]))
	      distemp0 <- distancia[temp0]
	      control.corte0 <- max(distemp0)
	      if(control.corte0 == control.corte) {
	        stop(paste("size not appropiated, max distance non different between consecutive breaks.
	                   Increase size"))
	      }
	      }
	    ##
	    
	    breakpoints[j] <- control.corte
	    meandist[j] <- mean(distemp)
	    logdist[j] <- mean(logdistemp)
	    dummy <- distancia
	    dummy <- dummy - distancia
	    dummy[temp] <- 1
	    laglw[[j]] <- dummy
	    cardinal[j] <- sum(dummy) / 2
	    j <- j+1
	    }
	  breakpoints <- c(min(distancia), breakpoints)
	  
	  if(self) {
	    dummy.self <- distancia - distancia
	    diag(dummy.self) <- 1
	    cardinal <- c(length(diag(dummy.self)), cardinal)
	    dummy.self <- list(dummy.self)
	    laglw <- append(dummy.self, laglw)
	    meandist <- c(0, meandist)
	    
	  }
	}
	
	
	#case k neighbors
	
	if(!is.null(kmax)) {
	  #computation based in k max
	  method <- "kmax"
	  smin <- NULL
	  cummulative <- TRUE
	  if(class(XY) == "dist") {
	    stop("XY is a distance matrix. kmax require a matrix XY with coordinates")
	  }
	  for(i in 1:kmax) {
	    laglw[[i]] <- (eco.weight(XY, method = "knearest", k=i))@W
	    npair <- sum(laglw[[i]])
	    meandist[i] <- sum(distancia * laglw[[i]]) / npair
	    logdist[i] <- sum(logdistancia * laglw[[i]]) / npair
	    cardinal[i] <- sum(laglw[[i]]) / 2
	  }
	}
	
	#control
	if(any(cardinal == 0)) {
	  stop("empty classes. Change parameters setting")
	}
	  
	#conversion to row standardized weights
	if(row.sd) {
	  laglist <- list()
	  laglist <- lapply(laglw, function(y) y/apply(y, 1, sum))
	  for(i in 1:length(laglist)) {
	    laglist[[i]][is.na(laglist[[i]])] <- 0
	  }
	} else {
	  laglist <-laglw
	}
	
	if(method != "kmax") {
	  if(!self){
	    breaks <- breakpoints 
	  } else {
	    breaks <- c(0, breakpoints)
	  }
	} else {
	  breaks <- 1:kmax
	}

	#parameters for seqvec
	
	if(!is.null(seqvec)) {
	  smax <- max(seqvec)
	  smin <- min(seqvec)
	  nclass <- length(seqvec) - 1
	}
	  
	param <- c("int", "smin", "smax", "kmax", "nclass", "size")
	param.val <- c(is.null(int), is.null(smin), is.null(smax), 
	               is.null(kmax), is.null(nclass), is.null(size))
	cuales <- which(!param.val)
	PAR <- param[cuales]
	PAR.VAL <- c(int, smin, smax, kmax, nclass, size)
	
	
	res <- new("eco.lagweight")
	res@W <- laglist
	res@XY <- data.frame(XY)
	res@PAR <- PAR
	res@PAR.VAL <- PAR.VAL
	res@ROW.SD <- row.sd
	res@SELF <- self
	res@CUMMUL <- cummulative
	res@MEAN <- meandist
	res@LOGMEAN <- logdist
	res@CARDINAL <- cardinal
	res@BREAKS <- breaks
	res@METHOD <- method
	
	
	res
	
          })

