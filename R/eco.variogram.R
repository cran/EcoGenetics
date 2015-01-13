#'Variogram for variables of an ecogen object.
#' @param eco An ecogen object.
#' @param var Variable to analyze (see description).
#' @param ... Further arguments passed to \code{\link[gstat]{variogram}}.
#' @description This program calls  \code{\link[gstat]{variogram}} 
#' and plots a variogram of a selected variable contained in an ecogen object.
#' 
#' For example, for selecting the variable P2 of the data frame P,
#' it should be written eco$P$P2, and so on. 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @examples 
#' \dontrun{
#' 
#' data(eco.test)
#' variog<-eco.variogram(eco, eco$P$P2, cutoff = 3000)
#' variog
#' 
#' }
#' @export

setGeneric("eco.variogram",  
					 function(eco, var, ...)  {
  
  loc<-eco@XY
  colnames(loc)<-c("x", "y")
  sp::coordinates(loc) <- ~ x + y
  variog <- gstat::variogram(var ~ 1, 
  													 locations = sp::coordinates(loc), 
                             data = loc, ...)
  
  x <- variog$dist
  y <- variog$gamma
  
  ylab <- "Semivariance"
  xlab <- "Distance"
  
  title <- deparse(substitute(var))
  
  plotfun<-function(){
  	ggplot2::ggplot(data = variog) + 
  	ggplot2::geom_line(ggplot2::aes(x = dist,y = gamma),
  										 directions="hv",linetype = 1, 
  										 colour = "red") +
  	ggplot2::geom_point(ggplot2::aes(x = dist,y = gamma),colour = "black") +
  	ggplot2::ylab(ylab) + 
  	ggplot2::xlab(xlab) + 
  	ggplot2::labs(title = title) +
  	ggplot2::theme(axis.title = ggplot2::element_text(size = 14,face = "bold"),
  								 plot.title = ggplot2::element_text(size = 16,face = "bold"))     
  }
  
  print(plotfun())
  
  variog
  
})
