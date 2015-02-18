# Scaling a data frame or matrix to 0 - 1 range
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015

setGeneric("eco.rescale", 
					 
					 function(dfm) {
						
	dfm <- as.data.frame(dfm)
  col <- apply(dfm, 2, function(X) { 
  	(X - min(X, na.rm = TRUE)) / diff(range(X, na.rm = TRUE))
               })
  return(col)
})
