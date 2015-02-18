# Subsetting an ecogen object by group
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015

setGeneric("eco.subset",
      
      function(object, fact, grp, missing = c(0, "NA",  "MEAN"), ...)  {
  
  grupo <- object$S
  x <- match(fact, colnames(object$S), nomatch = 0)
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
  z$P <- object$P[grupo, ]
  z$G <- object$G[grupo, ]
  z$E <- object$E[grupo, ]
  z$XY <- object$XY[grupo, ]
  z$S <- as.data.frame(object$S[grupo, ])
  z$GENIND <- df2genind(object$G[grupo, ], missing = missing, ...)
  
  colnames(z$S) <- colnames(object$S)
  
  z$C <- object$C[grupo, ]
  z$OUT <- list()
  z
})
