# Subsetting an ecogen object by group

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015 

setGeneric("eco.subset",
					 
					 function(eco, fact, grp, missing = c(0, "NA",  "MEAN"), ...)  {
  
  grupo <- eco$S
  x <- match(fact, colnames(eco$S), nomatch = 0)
  x <- x[!x == 0]
  
  missing <- match.arg(missing)
  
  if(length(x) == 0) {
    stop("incorrect factor name")
  }
  
  if(grp > max(as.numeric(grupo[, x]))) {
    stop(sprintf("the number of groups (%d) exceeds the number of
                 groups in the data (%d)", grp,
                 max(as.numeric(grupo[, x]))))
  }
  
  grupo <- which(grupo[, x] == grp)
  z <- ecogen()
  z$P <- eco$P[grupo, ]
  z$G <- eco$G[grupo, ]
  z$E <- eco$E[grupo, ]
  z$XY <- eco$XY[grupo, ]
  z$S <- as.data.frame(eco$S[grupo, ])
  z$GENIND <- df2genind(eco$G[grupo, ], missing = missing, ...)
  
  colnames(z$S) <- colnames(eco$S)
  
  z$C <- eco$C[grupo, ]
  z$OUT <- list()
  
  attr(z, "format") <- attr(eco, "format")
  attr(z, "type") <-  attr(eco, "type")
  attr(z, "missing") <- attr(eco, "missing")
  attr(z, "ploidy") <- attr(eco, "ploidy")
  
  z
})
