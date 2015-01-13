#' Creating an updated ecogen object by adding 
#' results to the slot OUT.
#' @param eco ecogen object.
#' @param ... objects to append to eco, typed without quotations. 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @examples
#' \dontrun{
#' 
#' data(eco.test)
#' variog <- eco.variogram(eco, eco$P$P1)
#' eco <- eco.append(eco, variog)
#' we.are.integers <- c(1:10)
#' we.are.characters <- c("John Coltrane", "Charlie Parker")
#' 
#' }
#' 
#' @export

setGeneric("eco.append", function(eco, ...) {
	

  res <- list(...)
  res.names <- (as.character(match.call()))
  res.names <- res.names[-c(1:2)]
  
  Z <- eco
  nob <- length(Z$OUT)
  nad <-  (nob+1):(nob + length(res))
  Z$OUT[nad] <- res
  names(Z$OUT)[nad] <- res.names
  orden <- order(names(Z$OUT))
  Z$OUT <- Z$OUT[orden]
  Z
})
