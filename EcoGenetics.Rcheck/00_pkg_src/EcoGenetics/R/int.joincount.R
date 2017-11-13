#' Join-count statistic, internal.
#' 
#' @param Z Vector, matrix or data frame.
#' @param con Connection network.
#' @param nsim Number of Monte-Carlo simulations. 
#' @param test If test = "bootstrap", the program generates a bootstrap 
#' resampling and the associated confidence intervals of the null hypothesis.
#'  If test = "permutation" (default) a permutation test is made and the p value 
#'  is calculated.    
#' @param alternative The alternative hypothesis. If "auto" is selected (default) the
#' program determines the hypothesis by difference between the median of the simulations
#' and the observed value. Other options are: "two.sided", "greater" and "less".
#' if test == cross, for the first interval (d== 0) the p and CI are computed with cor.test.
#' @param adjust.n Should be adjusted the number of individuals? (warning, this would
#' change variances)
#' @param adjustjc Method for multiple correction of P-values 
#' passed to \code{\link[stats]{p.adjust}}.
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @keywords internal

int.joincount <- function(Z, con, nsim,
                          alternative, test = "permutation", 
                          adjust.n = FALSE, adjustjc = "none", 
                          plotit) {
  

  con <- int.check.con(con)
  con <- con[lower.tri(con)]
  con <- as.vector(con)
  
  if(test == "permutation") {
    replace <- FALSE
  } else {
    replace <- TRUE
  }
  
  
  #internal function 
  
  
  outerMatrix <- function(input) {
    input <- as.vector(as.matrix(input))
    outmat <- outer(input, input, paste0)
    outmat <- aue.sort(outmat, ploidy = 2)
    outmat <- outmat[lower.tri(outmat)]
    as.factor(outmat)
  }
  
    
  jcfun <- function(x) {
  xlev <- levels(x)
  temp <- lapply(seq(along = xlev), function(i) {
    as.integer(x == xlev[i])
  })
  out <- sapply(temp, function(x) sum(x * con)) 
  names(out) <- xlev
  out
  }
    
  myOuter <- outerMatrix(Z)
  obs <- jcfun(myOuter)
  
  #simulated datasets. permuting rows and columns, mantaining structure


  cat("computing randomization test ...\n")
  monte.c <- sapply(1:nsim, function(i) {
    sample(myOuter, length(myOuter), replace = replace)
  })
  monte.c <- apply(monte.c, 2, function(x) jcfun(as.factor(x)))
  #  outsamp <- Z[samp]
  #  jcfun(outsamp)
  #})
  monte.c <- t(monte.c)
  ran <- int.random.test(repsim = monte.c, obs = obs, 
                         nsim = nsim, test = test,
                         alternative = alternative,
                         adjust = adjustjc)
  
  
  #labeling rows
  #outmat <- outer(Z, Z, FUN = "paste", sep = "")
  #rownames(ran)<- levels(as.factor(aue.sort(outmat, ploidy = 2))) 
  ran <- data.frame(colnames(monte.c), ran)
  colnames(ran)[1] <- "pairs"
  res <- list("analysis" = "Join-count", 
              "nsim" = nsim,
              "results" = ran)
  
  if(plotit) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    graphics::layout(matrix(rep(c(1, 1, 2,2,2,2, 2), 7), 7,7, byrow=TRUE))
    mycol <- ran$p.val < 0.05
    mycol <- mycol + 1
    mycol <- c("blue", "red")[mycol]
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    legend("topright", legend = c("P < 0.05", "NS"), fill = c("red", "blue"), cex = 1.2)
    barplot(ran$obs, col = mycol, names.arg = ran$pairs, ylab = "Frequency",
            xlab = "Pairs", cex.axis = 1.5, cex.names = 1.5, cex.lab = 1.5)
  }
  
  res
  
}
