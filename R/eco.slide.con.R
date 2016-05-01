
#' Sliding window for connection networks
#' 
#' @description This program applies a function defined by the user, over the individuals
#' included in a connection network.
#' For a given variable, the program computes recursively a function for the individual 
#' of the network, using all the individuals connected to each.
#' The function uses a connection network product of the function eco.weight
#' @param var Input variable or matrix
#' @param r half a side for square, radius for circle, diagonal length for rhombus.
#' @param slide number of elements between two focal pixels, for column 
#' and row dimension
#' @param fun Function to apply in each focal pixel.
#' @param window window type. Default "square".
#' @param within should be computed the function in borders focal pixels only if 
#' the area is within the matrix? Default TRUE.
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @export


eco.slide.con <- function(x, con, fun) { 
 var <- data.frame(x)
 nonzero <- lapply(1:nrow(con@W), function(i) which(con@W[i, ] != 0))
 res <- apply(var, 2, function(y) unlist(lapply(nonzero, function(z) fun(y[z]))))
 res
}
  
