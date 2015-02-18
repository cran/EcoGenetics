#' Empirical variogram
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015

setGeneric("eco.variogram",  
      function(z, xy, int, smax, w = c("B", "W"), nsim = 0, latlon = FALSE) {
 
 w <- match.arg(w)
 
  if(latlon == FALSE) {
 distancia <- dist(xy)
 } else {
   distancia <- latlon2distm(xy)
 }
       
 
 mat <- as.matrix(dist(z))
 
 d.max<- seq(int, smax, int)
 d.min <- d.max - int
 d.min[1] <- 1e-10
 classes <- length(d.min)
 
 dist.dat<-paste("d=", d.min, "-", d.max)
 dist.dat[1]<-paste("d=","0", "-", d.max[1])

 j <- 1
 d.mean <- numeric()
 dist2 <- dist(xy)
 for (i in seq(int, smax, int)) {
  temp <- which((distancia <= i) & (distancia > i - int))
  d.mean[j] <- round(mean(distancia[temp]), 3)
  j <- j+1
 }

 mat2 <- mat ^ 2
 wg <- eco.laglistw(xy, int, smax, w)
  wsub <- (2 * sapply(wg, sum))
 est <- sapply(wg, function(x) sum(x * mat2)) / wsub 
  
 
 if(nsim != 0) {
 boot.vario <- list()
 for(i in 1:nsim) {
 mat2 <- as.matrix(dist(sample(z, replace = TRUE)))
 mat2 <- mat2 ^ 2
 boot.vario[[i]] <- sapply(wg, function(x) sum(x * mat2)) / wsub
}
  boot.vario <- sapply(boot.vario, c)
 ext <- apply(boot.vario, 1, quantile, probs = c(0.05, 0.95),  na.rm = TRUE)
 tab <- data.frame(matrix(nrow = classes, ncol = 4))
 rownames(tab) <- dist.dat
 colnames(tab) <- c("d.mean","est", "lwr", "uppr")
 tab[, 1] <- d.mean
 tab[, 2] <- est
 tab[, 3:4] <- data.frame(t(ext))
 
 } else {
  tab <- data.frame(matrix(nrow = classes, ncol = 2))
  rownames(tab) <- dist.dat
  colnames(tab) <- c("d.mean","est")
  tab[, 1] <- d.mean
  tab[, 2] <- est
 }
 
 class(tab) <-"eco.variog"
 tab
})
  
