#' Mantel and partial Mantel tests, internal.
#' 
#' @param d1 distance matrix.
#' @param d2 distance matrix.
#' @param dc conditional distance matrix (optional).
#' @param nsim Number of Monte-Carlo simulations. 
#' @param test If test = "bootstrap", the program generates a bootstrap 
#' resampling and the associated confidence intervals.
#' If test = "permutation" a permutation test is made and the P-values 
#' are computed. 
#' @param alternative The alternative hypothesis. If "auto" is selected (default) the
#' program determines the alternative hypothesis.
#' Other options are: "two.sided", "greater" and "less". 
#' @param plotit Should be generated a plot of the simulations?  
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal


int.mantel <- function(d1, d2, dc, nsim, 
                       test, alternative = "auto", 
                       plotit = FALSE, method = "pearson", 
                       na.omit = FALSE, ...) {
  
  m1 <- as.double(d1)
  m2 <- as.double(d2)
  m3 <- if(!is.null(dc)) as.double(dc)
  
  m1m <- as.matrix(d1)
  N <- nrow(m1m)
  lowerTri <- which(row(m1m) > col(m1m))
  
  if(is.null(m3)) {
    
    obs <- cor(m1, m2, method = method, ...)
    
    repsim <- double(nsim)
    for(i in seq_len(nsim)) {
      samp <- sample(N)
      repsim[i] <- cor((m1m[samp, samp])[lowerTri], m2, method = method, ...)
    }
    
    
  } else {
    
    x23 <- cor(m2, m3, method = method, ...)
    
    
    corpartial <- function(x1, ...) {
      x12 <- cor(x1, m2, method = method, ...)
      x13 <- cor(x1, m3, method = method, ...)
      num <- (x12 - x13 * x23) 
      denom <- sqrt(1 - (x13 ^ 2)) * sqrt(1 - (x23 ^ 2))
      num / denom
    }
    
    obs <- corpartial(m1, ...)
    
    repsim <- double(nsim)
    for(i in seq_len(nsim)) {
      samp <- sample(N)
      repsim[i] <- corpartial((m1m[samp, samp])[lowerTri], ...)
    }
    
    
  }
  
  res <- int.random.test(repsim = repsim, 
                         obs = obs,
                         nsim = nsim,
                         test = test,
                         alternative = alternative)
  
  
  if(is.null(dc)) { 
    xlab <- "Mantel statistic"
  } else {
    xlab <- "Partial Mantel statistic"
  }
  
  
  if(plotit) {
    
    hist(c(repsim, obs), xlab = xlab,
         main = "Monte Carlo test")
    abline(v = obs, col = "red")
    points(obs, 0, col = "green", pch = 15, cex = 3.6)
    text(obs, 0, "obs")
  }
  
  res
}
