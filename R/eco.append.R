# Creating an updated ecogen object by adding results to the slot OUT
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015

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
