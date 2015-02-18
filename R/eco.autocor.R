# Moran's I, Geary's C, Join-Count, local Getis-Ord's G 
# and local Moran's I correlograms for an ecogen object
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015

setGeneric("eco.autocor",  
					 function(eco, int, smax,   
					 				 df = c("P", "G", "E", "GENIND", "C"),
					 				 select = c("moran", "geary", "joincount",
					 				 					 "getisord", "localmoran"),
					 				 nsim = 99,
					 				 indvar = NULL,
					 				 w = c("B", "W"),
					 				 latlon = FALSE,
					 				 fact = NA, grp = NA, 
					 				 ncod = NULL, ...) {
             
             
             # We start with some checks.
             
             cat("\n")
             cat(" analysis started at", as.character(Sys.time()),"\n")
             cat("\n")
             
             df <- match.arg(df)
             select <- match.arg(select)
             w <- match.arg(w)

             if(select == "joincount" & is.null(ncod) & df != "GENIND")  {
               stop("join count analysis requires a ncod argument when df is not 
                    the GENIND data frame")
             }
             
             if(select == "joincount" & df == "GENIND")  {
               ncod <- 1
             }
             
             
             xy <- eco@XY
             if(ncol(xy)>2) {
               cat("The XY data frame contains more that 2 columns.
                   (maybe altitude data, but it is ok). The program takes the 
                   first two columns as latitude -longitude information.", "\n\n")
               xy <- xy[, 1:2]
             }
             
             
             
             # Some specific checks and parameter settings 
             # for the case when a factor is included and when not.
             
             if(!is.na(fact)) {
               
               # the first mentioned case
               
               grupo <- eco@S
               fact <- match(fact, colnames(eco@S), nomatch = 0)
               fact <- fact[fact != 0]
               
               if(length(fact) == 0) {
                 stop("incorrect factor name")
               }
               
               if(is.na(grp)) {
                 stop("grp argument not provided")
               }
               
               where <- which(eco@S[, fact] == grp)
               xy <- eco@XY[where, ]
               if(latlon == FALSE) {
                 distancia <- dist(xy)
               } else {
                 distancia <- latlon2distm(xy)
               }
               
               hmuch <- sum(distancia < int)
               if(hmuch < 5) {
                 stop("Scale not apropiated.Increase distance interval")
               }
               
               hlast <- sum(distancia > smax - int)
               if(hlast <5 ) {
                 stop("Range not apropiated. Decrease smax value")
               }
               
               P <- eco@P[where, ]
               G <- eco@G[where, ]
               E <- eco@E[where, ]
               GENIND <- eco@GENIND$tab[where, ]
               C <- eco@C[where, ]
               
               
             } else {
               
               # the second case
               
               xy <- eco@XY
               if(latlon == FALSE) {
                 distancia <- dist(xy)
               } else {
                 distancia <- latlon2distm(xy)
               }
               
               hmuch <- sum(distancia < int)
               if(hmuch < 5) {
                 stop("Scale not apropiated.Increase distance interval")
               }
               hlast <- sum(distancia > smax - int)
               if(hlast < 5) {
                 stop("Range not apropiated. Decrease smax value")
               }
               
               
               P <- eco@P
               G <- eco@G
               E <- eco@E
               GENIND <- eco@GENIND$tab
               C <- eco@C
               
               
               }
             

             # Here we set the general parameters of the analysis 
             # and do some additional checks
             
             if(df == "P") {
               x <- P
             } else if(df == "G") {
               x <- G
             } else if(df == "E") {
               x <- E
             } else if(df == "GENIND") {
               x <- 2 * GENIND
             } else if(df == "C") {
               x <- C
             } else {
               stop("enter a valid data frame selection")
             }
             
             if(select == "joincount") {
               if(df == "G" && attributes(eco)$ploidy == 2) {
                 x <- eco.sort(x,  ncod / 2)
               }
             }
             
             
             if(!is.null(indvar)) {
               colx <- which(colnames(x) == indvar)
               x <- x[, colx, drop = FALSE]
             }
             
             nloc <- ncol(x)
             
             d.max<- seq(int, smax, int)
             d.min <- d.max - int
             d.min[1] <- 1e-10
             j <- 0
             
             # Computing distance intervals
             
             
             medint <- function(distancia, int, smax) {
               
               dist.range <- function(dism, int, smax)
               {
                 lista <- list()
                 j <- 1
                 for (i in seq(int, smax, int)) {
                   temp <- which((dism <= i) & (dism > i - int))
                   lista[[j]] <- temp
                   names(lista)[j] <- i
                   j <- j + 1
                 }
                 lista
               }
               rangos <- dist.range(distancia, int, smax)
               listamedias <- vector()
               dis <- distancia
               for(i in 1:length(rangos)) {
                 listamedias[i] <- mean(dis[rangos[[i]]])
               }
               listamedias
             }
             
             
             # Funtion to estimate the stat in each iteration
             
             select_method <- function(z, con, ...) {
               
               if(select == "moran") {
                 out <- eco.moran(z = z, con = con, nsim = nsim, 
                                  plotit =FALSE, test = "bootstrap",  ...)
               } else if(select == "geary") {
                 out <- eco.geary(z = z, con = con, nsim = nsim, 
                                  plotit =FALSE, test = "bootstrap", ...)
               } else if(select == "joincount") {
                 out <- eco.joincount(z = z, con = con, nsim = nsim, 
                                      ncod = ncod, test = "bootstrap", ...)
                 out <- out$results
               } else if(select == "getisord") {
               	out <- eco.getisord(z = z, con = con, nsim = nsim, 
               											test = "bootstrap", ...)
               	out <- out$results
               } else if(select == "localmoran") {
               	out <- eco.localmoran(z = z, con = con, nsim = nsim, 
               											test = "bootstrap", ...)
               	out <- out$results
               }
               out
             }
             
             
             # Iterating the latter with each individual variable
             
             lista<-list()
             
             lag <- eco.laglistw(xy, int, smax, w)
             
             #output data frame/s construction
             classint <- medint(distancia, int, smax)
             classint <- round(classint, 3)
             if(select == "moran" | select == "geary") {
               tabla <- data.frame(matrix(, length(d.min), 4))
               tabla[, 1] <- classint
               colnames(tabla) <- c("d.mean", "est", "lwr", "uppr")
               rownames(tabla) <- paste("d=", d.min, "-", d.max, sep = "")
               rownames(tabla)[1] <-paste("d=","0", "-", d.max[1], sep = "")
               
            } else if(select == "getisord" | select == "localmoran"){
             	tabla <- data.frame(matrix(, nrow(x), length(d.min)))
             	colnames(tabla) <- classint
             	rownames(tabla) <- rownames(x)
             	tabla <- replicate(3, tabla, simplify = FALSE)
             	names(tabla) <- c("observed", "lwr", "uppr")
             	lista <- replicate(ncol(x), tabla, simplify = FALSE)
             }
               
               for(j in 1:ncol(x)) {
                 var.test <- x[, j]
                 
                 #repetition of select_method for each run 
                 
                 if(select == "moran" | select == "geary") {
                   lista[[j]] <- tabla
                   for(i in 1:length(d.min))  {
                     cat("\r", "Computing", "trait", j," ", 
                         ceiling(i*100/length(d.min)), "%")
                     lag2 <- lag[[i]]
                     est <- select_method(var.test, lag2)
                     lista[[j]][i, 2] <- est$observation
                     lista[[j]][i, 3:4] <- est$quantile
                   }
                   names(lista)[j] <- colnames(x)[j]
                   cat("\n")
                   
                 } else if (select == "joincount") {
                 	outmat <- outer(var.test, var.test, FUN = "paste", sep = "")
                 	outmat <- (as.matrix(eco.sort(outmat, ncod)))
                 	outmat <- as.vector(outmat)
                 	outmat <- levels(as.factor(outmat))
                 	tabla <- data.frame(matrix(, length(outmat), length(d.min)))
                 	colnames(tabla) <- classint
                 	rownames(tabla) <- outmat
                 	tabla <- replicate(3, tabla, simplify = FALSE)
                 	names(tabla) <- c("observed", "lwr", "uppr")
                 	lista[[j]] <- tabla
                 	
                   for(i in 1:length(d.min)) {
                   
                     cat("\r", "Computing", "trait", j," ", 
                         ceiling(i*100/length(d.min)), "%")
                     lag2 <- lag[[i]]
                     est <- select_method(var.test, lag2)
                       lista[[j]][[1]][, i] <- est[ , 1]
                       lista[[j]][[2]][, i] <- est[ , 2]
                       lista[[j]][[3]][, i] <- est[ , 3]
                   }
                   cat("\n")
                 	
                 } else if(select == "getisord" | select == "localmoran") {
                 	for(i in 1:length(d.min))  {
                 		cat("\r", "Computing", "trait", j," ", 
                 				ceiling(i*100/length(d.min)), "%")
                 		lag2 <- lag[[i]]
                 		est <- select_method(var.test, lag2)
                 			lista[[j]][[1]][, i] <- est[, 1]
                 			lista[[j]][[2]][, i] <- est[, 2]
                 			lista[[j]][[3]][, i] <- est[, 3]
                 		}
                 	names(lista) <- colnames(x)
                 	}
                 cat("\n")
               }
                 	        
             # Configuring the output
             
            if (select == "moran" | select == "geary") {
            salida <- new("eco.boot")
            } else {
            salida <- new("eco.multiboot")
            }
            
             salida@OUT <- lista
             salida@NAMES <- colnames(x)
             salida@INTERVAL <- int
             salida@MAX <- smax
             salida@TYPE <- df
             salida@METHOD <- select
             salida@RANDTEST <- "bootstrap"

             cat("\n")
             cat("done!")
             cat("\n\n")
             
             ReturnVal <- tcltk::tkmessageBox(title = "Correlogram", 
                                              message = "process successful!",
                                              icon = "info", type = "ok")
             
             salida
             
               })
