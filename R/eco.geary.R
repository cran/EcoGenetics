# Geary's C statistic with Monte-Carlo test
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015


setGeneric("eco.geary", 
      function(z, con, nsim = 99,
               alternative = c("auto", "two.sided", 
                               "greater", "less"),
               test = c("permutation", "bootstrap"), 
               adjust.n = TRUE, plotit = TRUE) {
 
  alternative <- match.arg(alternative)
  test <- match.arg(test)
  
 n <- length(z)

 ccon <- class(con)[1]
 
 if(ccon == "listw") {
  listwg <- sapply(con$neighbours, c, simplify = FALSE)
  weig <- sapply(con$weights, c, simplify = FALSE)
  wg <- outer(z, z)
  wg[] <- 0
  for(i in 1:nrow(wg)) {
   wg[i, ][listwg[[i]]] <- weig[[i]]
  }
 } else if(ccon == "matrix"){
  wg <- con
 } else {
  stop("weight object provided is not of class listw or matrix")
 }
 
  if(adjust.n == TRUE) {
   colTRUE <- apply(wg, 1, sum)
   colTRUE[colTRUE != 0] <- 1
   n <- sum(colTRUE)
  }
  
 wg <- as.numeric(wg)

 mat <- as.matrix(dist(z, upper = T))
 mat <- as.numeric(mat)
 
 num <- mat ^ 2
 denom <- 2 * sum (wg) 
  denom1 <- denom * var(z)
 
 obs <- sum(wg * num) /  denom1
 
repsim <- numeric()
   if(test == "permutation") {
    for(i in 1:nsim) {
    mat2 <- as.matrix(dist(sample(z), upper = T))
    mat2 <- as.numeric(mat2)
    num <- mat2 ^ 2
      repsim[i] <- sum(wg * num) /denom1
    
    expected <- median(repsim)
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
    res <- list("analysis" = "Geary's C",
          "alternative"= paste("alternative hypothesis:", 
                     alternative), 
          "observation" = round(obs, 4),
          "expectation" = round(expected, 4),
          "nsim" = nsim, "p.value" = round(p.val, 4))
    }
    
   } else {
    for(i in 1:nsim) {
     z.boot <- sample(z, replace = TRUE)
     mat <- as.matrix(dist(z.boot, upper = T))
     mat <- as.numeric(mat)
     num <- mat ^ 2
     denom1 <- denom * var(z.boot)
     repsim[i] <- sum(wg * num) /denom1
     qq <- quantile(repsim, probs = c(0.05, 0.95),  na.rm = TRUE)
     res <- list("analysis" = "Geary's C", 
           "observation" = round(obs, 4), 
           "nsim" = nsim,
           "quantile" = round(qq, 4))
    }
   }

if(plotit == TRUE) {
 hist(c(repsim, obs), xlab = "Geary's C", main = "Monte Carlo test", xlim = )
 abline(v = obs, col = "red")
 points(obs, 0, col = "green", pch = 15, cex = 3.6)
 text(obs, 0, "obs")
}
 
  res


 
 })
