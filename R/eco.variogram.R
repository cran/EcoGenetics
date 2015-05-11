# Empirical variogram

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015 

setGeneric("eco.variogram",  
           function(Z, XY, 
                    int = NULL,
                    smin = 0,
                    smax = NULL,
                    nclass = NULL,
                    seqvec = NULL,
                    size = NULL,
                    bin = c("sturges", "FD"),
                    row.sd = FALSE,
                    latlon = FALSE) {
             
             bin <- match.arg(bin)
             
             #CHECKING XY DATA
             
             if(ncol(XY) > 2) {
               message("XY slot with > 2 columns. The first two are taken as X-Y coordinates")
               XY <- XY[,1:2]
             } 
             
             if(latlon == TRUE) {
               XY <- SoDA::geoXY(XY[,2], XY[,1], unit=1)
             } 
             
             
             ####
             
             mat <- as.matrix(dist(Z))
             
             listaw <- eco.lagweight(XY, 
                                     int = int, 
                                     smin = smin,
                                     smax = smax, 
                                     nclass = nclass,
                                     size = size,
                                     seqvec = seqvec,
                                     row.sd = row.sd,
                                     bin = bin)
             
             wg <- listaw@W
             
             breakpoints <- listaw@BREAKS
             d.max <- round(breakpoints[-1], 3)
             d.min <- round(breakpoints[-length(breakpoints)], 3)
             classint <- listaw@MEAN
             classint <- round(classint, 3)
             cardinal <- listaw@CARDINAL
             
             dist.dat<-paste("d=", d.min, "-", d.max)
             
             d.mean <- listaw@MEAN
             
             mat2 <- mat ^ 2
             wsub <- (2 * sapply(wg, sum))
             est <- sapply(wg, function(x) sum(x * mat2)) / wsub 
             
             
             tab <- data.frame(matrix(nrow = length(dist.dat), ncol = 2))
             rownames(tab) <- dist.dat
             tab[, 1] <- d.mean
             tab[, 2] <- est
             colnames(tab) <- c("d.mean","obs")
             
             salida <- new("eco.correlog")
             
             salida@OUT <- list(tab)
             salida@IN <- list(XY = XY, Z = Z)
             salida@BREAKS <- breakpoints
             salida@CARDINAL <- cardinal
             salida@METHOD <- "empirical variogram"
             salida@DISTMETHOD <- listaw@METHOD
             
             salida
             
           })

