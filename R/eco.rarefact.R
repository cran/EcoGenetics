#' Multiple rarefaction estimates via bootstrap for ecogen genetic data frames.  
#' @param eco ecogen object.
#' @param x The name of the S slot column with the groups for
#' the rarefaction analysis.
#' @param nrep Number of repetitions for the bootstrap estimates.
#' @param mode Mode passed to the function \code{\link[gstudio]{rarefaction}}.
#' Default: "Ae".
#' @param type Argument passed to \code{\link[gstudio]{locus}}. If not
#' provided, the program set it as "separated" for data as microsatellites 
#' or individual haplotypes, and "aflp" for presence - absence data.
#' @description This program computes multiple-loci confidence intervals 
#' of 95\% by bootstrap,  passing arguments to the functions 
#' \code{\link{rarefaction}} and \code{\link{genetic_diversity}}.
#' The estimates are calculated by group. The metrics which could be 
#' estimated are specified in \code{\link{genetic_diversity}}. 
#' A bootstrap over groups also is performed. In the case of Ae, the mean 
#' is replaced by the harmonic mean.
#' @return ESTIMATES BY LOCI the estimates by loci and the 95\% CI.
#' @return OVERALL ESTIMATES the overall estimates and the 95\% CI.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @examples
#'  \dontrun{
#'  
#' data(eco.test)
#' Ae.raref <- eco.rarefact(eco,"structure", nrep = 50)
#' 
#' }
#' @export

setGeneric("eco.rarefact", 
					 
					 function(eco, x, 
					 				 nrep = 99, mode = "Ae", 
					 				 type = NA)  {
  
													 	
  
  cat("\n")
  cat(" analysis started at", as.character(Sys.time()),"\n")
  cat("\n")
  
  grupo <- eco@S
  fact <- match(x, colnames(eco@S), nomatch = 0)
  fact <- fact[fact != 0]
  
  if(length(fact) == 0) {
    stop("incorrect factor name")
  }
  
  
  npop <- max(as.numeric(eco@S[, fact]))
  
  datoshier <- eco.2hierfstat(eco, x)
  
  if(is.na(type)) {
  	if(eco@GENIND$type == "codom")
  	{
  		type <- "separated"
  	} else {
  		type <- "aflp"
  	}
  }
  
  dat <- eco.2gstudio(eco, type = type)
  
  raref <- intervalos <- list()
  nmar <- ncol(datoshier) - 1
  estimas <- matrix(ncol = npop, nrow = nmar)
  
  
  
  for(k in 1:nmar) {
    
    cat("\r","bootstrap over individuals ","|| ")
    cat("loci", k, "00 %", "||", as.character(Sys.time()), "", sep = " ")
    
    
    raref[[k]] <- as.data.frame(matrix(nrow = nrep, ncol = npop))
    for(i in 1:npop) {  
      ind <- rep(0, npop)
      for(j in 1:npop) {
        ind[j] <- sum(is.na(datoshier[grupo[, fact] == j, i + 1]))
        ind[j] <- sum(datoshier[, 1] == j) - ind[j]
      }
      ind <- min(ind)
      
      raref[[k]][,i] <- gstudio::rarefaction(dat[grupo[, fact] == i, k],
                                             mode = mode, size = ind,
                                             nperm = nrep)
      estimas[,i] <- gstudio::genetic_diversity(dat[grupo[, fact]==i, ],
                                                mode = mode)[, 2]
      
      cat ("\r", "bootstrap over individuals ", "|| ")
      cat("loci", k, formatC(ceiling((100 * i) / npop), width = 2,
                             flag = "0"), "%", "|| ")
      cat(as.character(Sys.time()), "         ", sep = " ") 
    }
    intervalos[[k]] <- sapply(raref[[k]], quantile, probs = c(0.05, 0.95))
  }
  
  cat ("\r", "bootstrap over individuals ", "|| ")
  cat("loci", k, "100", "%", "|| ")
  cat(as.character(Sys.time()), "         ", sep = " ")  
  cat("\n\n")
  
  
  
  if(nmar<4) {
    warning("warning, bootstrapping over less than 4 loci should be not
            correct!")
  } 
  
  
  temp <- estimas
  promedio <- data.frame(matrix(, nrep, npop))
  colnames(promedio) <- 1:npop
  
  for(i in 1:nrep) {
    temp <- data.frame(estimas[sample(1:nmar, nmar, TRUE), ])
    
    if(mode == "Ae") {
      fun <- function(a) 1 / mean(1 / a)
      promedio[i,] <- apply(temp, 2, fun)
    } else {
      promedio[i, ] <- apply(temp, 2, mean)
    }
    
    cat ("\r", " bootstrap over loci ", " ", 
         format(round(100 *(i / nrep), 1), trim = T, nsmall = 1),
         "% ", as.character(Sys.time()), "             ", sep = "")
  }
  
  intervalos.promedio <- apply(promedio, 2, quantile, probs = c(0.05, 0.95))
  medias <- rbind(apply(promedio, 2, function(a) 1 / mean(1 / a)), 
                  intervalos.promedio)
  rownames(medias) <- c("estimates", "5%", "95%")
  cat("\n\n")
  
  colnames(estimas) <- 1:npop
  rownames(estimas) <- colnames(eco@G) 
  
  
  intervalos2 <- as.data.frame(intervalos[[1]])
  temp <- rep(colnames(eco@G)[1], 2)
  rownames(intervalos2) <- paste(temp, ".", rownames(intervalos2), sep = "")
  
  
  for(i in 2:length(intervalos)) {
    temp <- rep(colnames(eco@G)[i], 2)
    rownames(intervalos[[i]]) <- paste(temp,".", rownames(intervalos[[i]]),
                                       sep = "")
    intervalos2 <- rbind(intervalos2, intervalos[[i]])
  }
  
  colnames(intervalos2) <- 1:npop
  
  final <- list("ESTIMATES_BY_LOCI" = estimas,"CI_BY_LOCI" = intervalos2, 
                "OVERALL_ESTIMATES" = medias)
  
  cat("\n")
  cat("done!")
  cat("\n\n")
  
  ReturnVal <- tcltk::tkmessageBox(title = "Bootstrap", 
                                   message = "process successful!",
                                   icon = "info", type = "ok")
  
  final
  })
