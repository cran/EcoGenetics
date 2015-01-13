#' Scaling a data frame or matrix to 0 - 1 range.
#' @param dfm Dataframe, matrix or vector to scale.
#' @description The program scales each column of a data frame or a matrix 
#' to 0 - 1 range, computing (X\emph{ij} - Xmin\emph{i}) / range(X)\emph{i} 
#' for each individual \emph{j} of the variable \emph{i}.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar} 
#' @examples
#' \dontrun{
#' 
#' data(eco.test)
#' require(adegenet)
#' pc <- dudi.pca(eco$P, scannf = FALSE, nf = 3)
#' pc <- pc$li
#' pc <- eco.rescale(pc)
#' plot(eco$XY[, 1], eco$XY[, 2], col = rgb(pc), pch = 16, cex = 1.5, xlab ="X", ylab= "Y")
#' 
#' }
#' @export

setGeneric("eco.rescale", 
					 
					 function(dfm) {
						
	dfm <- as.data.frame(dfm)
  col <- apply(dfm, 2, function(X) { 
  	(X - min(X, na.rm = TRUE)) / diff(range(X, na.rm = TRUE))
               })
  return(col)
})
