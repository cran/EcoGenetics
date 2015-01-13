#' Creating an updated ecogen object by removing
#' results of the slot OUT.
#' @param object Ecogen object.
#' @param ... Objects to remove, typed without quotations. 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @examples
#' \dontrun{
#' 
#' data(eco.test)
#' variog <- eco.variogram(eco, eco@@P$P1)
#' eco <- eco.append(eco, variog)
#' we.are.numbers <- c(1:10)
#' we.are.characters <- c("John Coltrane", "Charlie Parker")
#' eco <- eco.append(eco, we.are.numbers, we.are.characters)
#' eco
#' eco <- eco.remove(eco, we.are.numbers)
#' eco
#' 
#' }
#' @export


setGeneric("eco.remove", 
					 
					 function(object, ...) {
  
  res.names <- as.character(match.call())
  res.names <- res.names[-c(1:2)]
  del <- (names(object@OUT)) %in% res.names
  object@OUT <- object@OUT[!del]
  object
})
