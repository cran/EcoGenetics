# Mantel correlogram
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015


setGeneric("eco.cormantel", 
           function(z, xy, int, smax, 
           				  nsim = 99, 
           				  latlon = FALSE, 
           				  adjust = c("none", "holm", "hochberg", 
           				 					   "hommel", "bonferroni", "BH",
           				 					   "BY", "fdr"), 
                    alternative = c("auto", "two.sided", 
                                    "greater", "less"),
                    test = c("permutation", "bootstrap"), ...) {
             
             
             alternative.i <- match.arg(alternative)
             adjust <- match.arg(adjust)
             test <- match.arg(test)
             
             if(ncol(xy)>2) {
               cat("The XY data frame contains more that 2 columns.
                   (maybe altitude data, but it is ok). The program takes the 
                   first two columns as latitude -longitude information.", "\n\n")
               xy <- as.matrix(xy[,1:2])
             } 
             
             hmuch <- sum(dist(xy) < int)
             if(hmuch < 5) {
               stop("Scale not apropiated.Increase distance interval")
             }
             hlast <- sum(dist(xy) > smax - int)
             if(hlast < 5) {
               stop("Range not apropiated. Decrease smax value")
             }
             
             
             d.max<- seq(int, smax, int)
             d.min <- d.max - int
             d.min[1] <- 1e-10
             
             dist.dat<-paste("d=", d.min, "-", d.max)
             dist.dat[1]<-paste("d=","0", "-", d.max[1])
             
             if(class(z) != "dist") {
               m1 <- as.numeric(dist(as.data.frame(z), ...))
             } else { 
               m1 <- as.numeric(z)
             }
             
             if(latlon == FALSE) {
               distancia <- as.numeric(dist(xy))
             } else {
               distancia <- as.numeric(latlon2distm(xy))
             }
             
             listamedias <- vector()
             j <- 1
             
             cat("\n", "wait...","\n\n")
             
             if(test == "bootstrap") {
               tab <- data.frame(matrix(nrow = length(d.min), ncol=4))
               rownames(tab) <- dist.dat
               colnames(tab) <- c("d.mean","est", "lwr", "uppr")
               
               for (i in seq(int, smax, int)) {
                 temp <- which((distancia <= i) & (distancia > i - int))
                 d.mean <- mean(distancia[temp])
                 dummy <- distancia
                 dummy <- dummy - distancia
                 dummy[temp] <- 1
                 
                 obs <- -cor(m1, dummy)  #sign changed. Also for the replicates
                 
                 rep <-  -replicate(nsim, cor(sample(m1, replace = T), dummy))
                 ext <- quantile(rep, probs = c(0.05, 0.95),  na.rm = TRUE)
                 
                 tab[j,] <-c(d.mean, obs, ext)
                 j <- j + 1
                 
               }
               
               salida <- new("eco.boot")
               
             } else if(test == "permutation") {
               tab <- data.frame(matrix(nrow = length(d.min), ncol= 4))
               rownames(tab) <- dist.dat
               colnames(tab) <- c("d.mean","est", "exp", "pval")
               
               for (i in seq(int, smax, int)) {
                 temp <- which((distancia <= i) & (distancia > i - int))
                 d.mean <- mean(distancia[temp])
                 dummy <- distancia
                 dummy <- dummy - distancia
                 dummy[temp] <- 1
                 
                 obs <- -cor(m1, dummy)
                 
                 repsim <- -replicate(nsim, cor(sample(m1), dummy))
                 
                 alternative <- alternative.i
                 
                 expected <- mean(repsim)
                 if(alternative == "auto") {
                   alter <- expected - obs
                   if(alter > 0) {
                     alternative <- "less"
                   } else if (alter < 0) {
                     alternative <- "greater"
                   }
                 }
                 
                 if(alternative == "greater") {
                   howmany <- sum(repsim >= obs)
                   p.val <- (howmany + 1) / (nsim + 1)
                   
                 } else if(alternative == "less") {
                   howmany <- sum(repsim <= obs)
                   p.val <- (howmany + 1) / (nsim + 1)
                   
                 } else if(alternative == "two.sided") {
                   howmany <- sum(abs(repsim) >= abs(obs))
                   p.val <- (howmany + 1) / (nsim + 1)
                 }
                 tab[j,] <-c(round(d.mean, 3), round(obs, 4),
                 						round(expected, 4), p.val)
                 j <- j + 1
               }
               tab[, 3] <- round(p.adjust(tab[, 3], method = adjust), 4)
               
               salida <- new("eco.permut")
             }
             
             
             salida@OUT <- list (tab)
             salida@NAMES <- colnames(as.data.frame(as.matrix(z)))
             salida@INTERVAL <- int
             salida@MAX <- smax
             salida@TYPE <- "dfm"
             salida@METHOD <- "mantel"
             salida@RANDTEST <- test
           
             salida
           })
