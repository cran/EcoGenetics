# Moran's I, Geary's C and bivariate Moran's Ixy correlograms

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015

setGeneric("eco.correlog", 
           function(Z, XY, Y = NULL, 
                    int = NULL,
                    smin = 0,
                    smax = NULL, nclass = NULL,
                    size = NULL,
                    seqvec = NULL,
                    method = c("I", "C", "CC"),
                    nsim = 99,
                    test = c("permutation", "bootstrap"),
                    alternative = c("auto", "two.sided", 
                                    "greater", "less"),
                    adjust = "holm",
                    sequential = TRUE, 
                    include.zero = TRUE,
                    cummulative = FALSE,
                    bin = c("sturges", "FD"),
                    row.sd = FALSE,
                    latlon = FALSE) {
             
             
             # We start with some checks.
             
             
             cat("\n")
             
             method <- match.arg(method)
             alternative <- match.arg(alternative)
             test <- match.arg(test)  
             bin <- match.arg(bin)
             
             
             
             Z.class <- class(Z)
             if(Z.class != "numeric" & Z.class != "integer" &  Z.class != "matrix" & Z.class != "data.frame") {
               stop("Z must me a numeric vector, a matrix or a data.frame")
             }
             if(!is.null(Y)) {
               Y.class <- class(Y)
               if(Y.class != "numeric" & Y.class != "integer" & Y.class != "vector") {
                 stop("Y must me a numeric vector")
                 if(is.null(method)) {
                   method <- "CC"
                 }
               }
             }
             
             
             Z <- as.data.frame(Z)
             nvar <- ncol(Z)
             
             
             if(method != "CC") {
               include.zero = FALSE
             }
             
             if(ncol(XY) > 2) {
               message("XY slot with > 2 columns. The first two are taken as X-Y coordinates")
               XY <- XY[,1:2]
             } 
             
             if(latlon == TRUE) {
               XY <- SoDA::geoXY(XY[,2], XY[,1], unit=1)
             } 
             
             distancia <- dist(XY)
             
             
             if(is.null(smax) & is.null(nclass) & is.null(seqvec)) {
               smax <- max(distancia)
             }
             
             if(!is.null(int) & !is.null(smax)) {
               
               hmuch <- sum(distancia > 0 & distancia < int)
               if(hmuch < 5) {
                 stop("Scale not apropiated.Increase distance interval")
               }
               hlast <- sum(distancia > smax - int)
               if(hlast < 5) {
                 stop("Range not apropiated. Decrease smax value")
               }
             }
             
             #additional check when class(XY) == "dist"
             if(!is.null(smax)) {
               if(smax > max(distancia)) {
                 stop("scale not apropiated. Decrease smax")
               }
             }
             
             
             j <- 0
             
             
             # Funtion to estimate the stat in each iteration
             
             select_method <- function(u, ...) {
               
               if(method == "I") {
                 out <- int.moran(Z = u,  ...)
               } else if(method == "C") {
                 out <- int.geary(Z = u,  ...)
               } else if(method == "CC") {
                 out <- int.crosscor(Z = u, Y = Y, ...)
               }    
               out
             }
             
             
             # Iterating the latter with each individual variable
             
             lista<-list()
             listaw <- eco.lagweight(XY, 
                                     int = int, 
                                     smin = smin,
                                     smax = smax, 
                                     nclass = nclass,
                                     size = size,
                                     seqvec = seqvec,
                                     row.sd = row.sd,
                                     bin = bin,
                                     cummulative = cummulative)
             
             lag <- listaw@W
             breaks<- listaw@BREAKS
             
             d.max <- round(breaks[-1], 3)
             d.min <- round(breaks[-length(breaks)], 3)
             classint <- listaw@MEAN
             classint <- round(classint, 3)
             cardinal <- listaw@CARDINAL
             
             #output data frame/s construction
             
             #bootstrap case
             if(test == "bootstrap") {
               
               tabla <- data.frame(matrix(, length(d.min), 4))
               tabla[, 1] <- classint
               colnames(tabla) <- c("d.mean", "obs", "lwr", "uppr", "size")
               rownames(tabla) <- paste("d=", d.min, "-", d.max, sep = "")
               lista <- replicate(nvar, tabla, simplify = FALSE)
               names(lista) <- colnames(Z)
               
               #repetition of select_method for each run 
               
               
               for(j in 1:nvar) {
                 var.test <- Z[, j]
                 
                 for(i in 1:length(d.min))  {
                   cat("\r", "Computing",  
                       ceiling(i*100/length(d.min)), "%", "trait", j)
                   lag2 <- lag[[i]]
                   est <- select_method(u = var.test, 
                                        con = lag2, 
                                        nsim = nsim,
                                        alternative = alternative,
                                        test =  test,
                                        plotit = FALSE)
                   lista[[j]][i, 2] <- est$observation
                   lista[[j]][i, 3:4] <- est$quantile
                 }
                 lista[[j]][, 4] <- cardinal
                 
                 #when zero is included, use cor.test for d = 0
                 if(include.zero) {
                   cor.zero <- cor.test(var.test, Y)
                   lista[[j]] <- rbind(c(0, 0, 0, 0, 0), lista[[j]])
                   rownames(lista[[j]])[1] <- "d=0"
                   lista[[j]][1, 1] <- 0
                   lista[[j]][1, 2] <- cor.zero$estimate
                   lista[[j]][1, 3:4] <- cor.zero$conf.int
                   lista[[j]][1, 5] <- nrow(Z) 
                 }
                 cat("\n")
               }
               
               #permutation case
             } else if(test == "permutation") {
               
               
               tabla <- data.frame(matrix(0, length(d.min), 4))
               tabla[, 1] <- classint
               colnames(tabla) <- c("d.mean", "obs", "pval", "size")
               rownames(tabla) <- paste("d=", d.min, "-", d.max, sep = "")
               lista <- replicate(nvar, tabla, simplify = FALSE)
               names(lista) <- colnames(Z)
               
               #repetition of select_method for each run 
               
               
               for(j in 1:nvar) {
                 var.test <- Z[, j]
                 
                 for(i in 1:length(d.min))  {
                   cat("\r", "Computing",  
                       ceiling(i*100/length(d.min)), "%", "trait", j)
                   lag2 <- lag[[i]]
                   est <- select_method(u = var.test, 
                                        con = lag2, 
                                        nsim = nsim,
                                        alternative = alternative,
                                        test =  test,
                                        plotit = FALSE)
                   lista[[j]][i, 2] <- est$observation
                   lista[[j]][i, 3] <- est$p.value
                 }
                 lista[[j]][, 4] <- cardinal
                 #when zero is included, use cor.test for d = 0
                 if(include.zero) {
                   cor.zero <- cor.test(var.test, Y)
                   lista[[j]] <- rbind <- rbind(c(0, 0, 0, 0), lista[[j]])
                   rownames(lista[[j]])[1] <- "d=0"
                   lista[[j]][1, 1] <- 0
                   lista[[j]][1, 2] <- cor.zero$estimate
                   lista[[j]][1, 3] <- cor.test(var.test, Y)$p.value
                   lista[[j]][1, 4] <- nrow(Z)
                 }
                 
                 #sequential P correction
                 if(sequential) {
                   for(i in 1:length(d.min)) {
                     lista[[j]][i, 3] <- (p.adjust(lista[[j]][1:i, 3], 
                                                   method= adjust))[i]
                   }
                   
                 } else {
                   #standard-multiple P correction 
                   lista[[j]][ , 3] <- p.adjust(lista[[j]][ , 3], 
                                                method = adjust)
                 }
                 
                 cat("\n")
               }
               
             }
             
             
             # Configuring the output
             
             
             salida <- new("eco.correlog")
             
             
             if(method == "I") {
               outname <- "Moran's I"
             } else if(method == "C") {
               outname <- "Geary's C"
             } else if(method == "CC") {
               outname <- "Moran's Ixy"
             }
             
             
             salida@OUT <- lista
             salida@IN <- list(XY = XY, Z = Z,  Y =Y)
             salida@NAMES <- names(lista)
             salida@BREAKS <- breaks
             salida@CARDINAL <- cardinal
             salida@METHOD <- outname
             salida@DISTMETHOD <- listaw@METHOD
             salida@TEST <- test
             salida@NSIM <- nsim
             salida@PADJUST <- paste(adjust, "-sequential:", sequential)
             
             cat("\n")
             cat("done!")
             cat("\n\n")
             
             
             salida
             
           })
