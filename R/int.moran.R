# Moran internal

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015 

int.moran <- function(Z, con, nsim,
                      alternative, 
                      test = "permutation", adjust.n = FALSE, 
                      plotit) {
  
  N <- length(Z)
  wg <- int.check.con(con)
  
  #weight adjustment to number of connections
  if(adjust.n == TRUE) {
    colTRUE <- apply(wg, 1, sum)
    colTRUE[colTRUE != 0] <- 1
    Nc <- sum(colTRUE)
  } else  {
    Nc <- N
  }
  
  #Moran's I computation 
  z2 <- Z - mean(Z)
  SC <- drop(z2 %*% z2)
  VAR.Z <- SC / Nc
  W <- sum(wg)
  
  moranfun <- function(zc) {
    SCW <- wg %*% zc
    AUTOCOV.Z <- drop(zc %*% SCW) / W
    res <- AUTOCOV.Z / VAR.Z
    res
  }
  
  #observed value
  obs <- moranfun(z2)
  
  #Monte carlo replicates
  repsim <- numeric()
  for(i in 1:nsim) {
    samp <- sample(N)
    coefsup.mc <- z2[samp]
    repsim[i] <- moranfun(coefsup.mc)
  }
  
  
  #p value or CI computation
  random.m <- int.random.test(repsim = repsim, 
                              obs = obs,
                              nsim = nsim,
                              test = test, 
                              alternative = alternative)
  
  
  #listing results
  if(test == "permutation") {
    
    res <- list("analysis" = "Moran's I", 
                "alternative"= random.m$alter, 
                "observation" = round(random.m$obs, 4),
                "expectation" = round(random.m$exp, 4),
                "nsim" = nsim,
                "p.value" = round(random.m$p.val, 5),
                "quantile" = round(random.m$CI, 4))
    
  }else {
    
    res <- list("analysis" = "Moran's I",
                "observation" = round(random.m$obs, 4),
                "nsim" = nsim,
                "quantile" = round(random.m$CI, 4))
  }
  
  
  #plot
  if(plotit == TRUE) {
    hist(c(repsim, random.m$obs),
         xlab = "Moran's I", 
         main = "Monte Carlo test")
    abline(v = obs, col = "red")
    points(obs, 0, col = "green", 
           pch = 15, cex = 3.6)
    text(random.m$obs, 0, "obs")
  }
  
  res
  
} 
