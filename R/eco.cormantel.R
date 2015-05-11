# Mantel and partial Mantel correlograms

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015

setGeneric("eco.cormantel", 
           function(M, XY, MC = NULL, int = NULL, smin = 0,
                    smax =NULL, nclass = NULL, seqvec = NULL,
                    size = NULL,  bin = c("sturges", "FD"), 
                    nsim = 99, classM = c("dist", "simil"),
                    method = c("pearson", "spearman", "kendall"),
                    test = c("permutation", 
                             "bootstrap"),
                    alternative = c("auto", "two.sided", 
                                    "greater", "less"),
                    adjust = "holm", 
                    sequential = TRUE, 
                    latlon = FALSE,
                    ...) {
             
             
             alternative.i <- match.arg(alternative)
             test <- match.arg(test)
             bin <- match.arg(bin)
             classM <- match.arg(classM)
             method <- match.arg(method)
             
             #some check of the data
             
             #XY CHECKING
             if(ncol(XY)>2) {
               message(paste("XY with > 2 columns. The
                             first two are taken as X-Y coordinates"))
               XY <-XY[,1:2]
             } 
             
             if(latlon == TRUE) {
               XY <- SoDA::geoXY(XY[,2], XY[,1], unit=1)
             }
             distancia <- dist(XY)
             
             #####
             
             #LAG PARAMETERS
             
             
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
             
             #range parametrization
             
             listaw <- eco.lagweight(XY, 
                                     int = int, 
                                     smin = smin,
                                     smax = smax, 
                                     nclass = nclass,
                                     size = size,
                                     seqvec = seqvec,
                                     row.sd = FALSE,
                                     bin = bin,
                                     cummulative = FALSE)
             
             lag <- listaw@W
             
             d.mean <- listaw@MEAN
             d.mean <- round(d.mean, 3)
             cardinal <- listaw@CARDINAL
             
             breakpoints <- listaw@BREAKS
             
             d.max <- round(breakpoints[-1], 3)
             d.min <- round(breakpoints[-length(breakpoints)], 3)
             
             dist.dat<-paste("d=", d.min, "-", d.max)
             lengthbreak <- length(breakpoints) - 1
             
             if(is.null(MC)) {
               method.mantel <- "Mantel statistic"
             } else {
               method.mantel <- "Partial Mantel statistic"
             }
             
             cat("\r", "interval", 0,"/", lengthbreak , "completed")
             
             
             #starting the computation of the statistic
             
             ####bootstrap case ####
             if(test == "bootstrap") {
               tab <- data.frame(matrix(nrow = length(d.min), ncol=5))
               rownames(tab) <- dist.dat
               colnames(tab) <- c("d.mean","obs", "lwr", "uppr", "size")
               
               for(i in 1:length(lag)) {
                 
                 #mantel test
                 
                 result <- int.mantel(d1 = M, d2 = as.dist(lag[[i]]), 
                                      dc = MC, nsim = nsim, test = "bootstrap", 
                                      method = method, ...)
                 
                 obs <- result$obs
                 ext <- result$CI
                 ext1 <- ext[1]
                 ext2 <- ext[2]
                 
                 # change of sign for "dist" data
                 if(classM == "dist") {
                   obs <- - obs
                 }
                 
                 tab[i, ] <-c(d.mean[i],
                              round(obs, 4),
                              round(ext1, 4),
                              round(ext2, 4),
                              cardinal[i])
                 
                 cat("\r", "interval", i,"/", lengthbreak , "completed")
                 
               }
               
               ####permutation case ####
             } else if(test == "permutation") {
               tab <- data.frame(matrix(nrow = length(d.min), ncol= 5))
               rownames(tab) <- dist.dat
               colnames(tab) <- c("d.mean","obs", "exp", "pval", "cardinal")
               
               for(i in 1:length(lag)) {
                 
                 result <- int.mantel(d1 = M, d2 = as.dist(lag[[i]]), 
                                      dc = MC, nsim,  test = "permutation",
                                      alternative = alternative, 
                                      method = method, ...)
                 
                 obs <- result$obs
                 expected <- result$exp
                 p.val <- result$p.val
                 
                 # change of sign for "dist" data
                 if(classM == "dist") {
                   obs <- - obs
                   expected <- - expected
                 }
                 
                 tab[i,] <-c(round(d.mean[i], 3),
                             round(obs, 4),
                             round(expected, 4), 
                             round(p.val, 5),
                             cardinal[i])
                 
                 cat("\r", "interval", i,"/", lengthbreak , "completed")
                 
               }
               
               
               #sequential correction
               
               if(sequential) {
                 for(j in 1:nrow(tab)) {
                   tab[j, 4] <- (p.adjust(tab[1:j, 4], method= adjust))[j]
                 }
               } else {
                 tab[ , 4] <- p.adjust(tab[ , 4], method = adjust)
               }
               
             }
             
             salida <- new("eco.correlog")
             salida@OUT <- list (tab)
             salida@IN <- list(XY = XY, M = M, MC = MC)
             salida@BREAKS <- breakpoints
             salida@CARDINAL <- cardinal
             salida@METHOD <- c(method.mantel, method)
             salida@DISTMETHOD <- listaw@METHOD
             salida@TEST <- test
             salida@NSIM <- nsim
             salida@PADJUST <- paste(adjust, "-sequential:", sequential)
             
             salida
           })
