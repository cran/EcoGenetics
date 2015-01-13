#' Converting an ecogen genetic data frame into a hierfstat 
#' data frame.
#' @param eco ecogen object.
#' @param fact The name of the S slot column with the groups 
#' for output data frame.
#' @description This program creates a hierfstat data frame 
#' with an ecogen object.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @examples
#' \dontrun{
#' 
#' data(eco.test)
#' hiereco <- eco.2hierfstat(eco, "structure")
#' require("hierfstat")
#' basic.stats(hiereco)
#' 
#' }
#' @export
#' 
setGeneric("eco.2hierfstat", 
					 function(eco, fact = NULL) {
	
u <- eco$G

grupo <- eco@S

if(is.null(fact))
{
  factord <- 	as.data.frame(rep(1, nrow(u)))
  cnom <- "pop"
  rnom <- rownames(eco@G)
  Gord <- u
} else {
	
fact <- match(fact, colnames(eco@S), nomatch = 0)
fact <- fact[fact != 0]
if(length(fact) == 0) {
  stop("incorrect factor name")
}
orden <- order(eco@S[, fact])
Gord <- u[orden,]
factord <- eco@S[orden, fact]
factord <- as.numeric(factord)
cnom <- colnames(eco@S[fact])
rnom <- rownames(eco@G)[orden]
}

datahier <- data.frame(factord, Gord)
colnames(datahier)[1] <- cnom
rownames(datahier) <- rnom
datahier

})
