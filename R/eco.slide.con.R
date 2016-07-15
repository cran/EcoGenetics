
#' Sliding window along a network
#' 
#' @description This program applies a function defined by the user, over the individuals
#' included in a connection network.
#' For a given variable, the program computes recursively a function for the individuals 
#' of the network, using all the individuals connected to each.
#' The function uses a connection network generated with the function eco.weight.
#' @param x Input variable or matrix.
#' @param con Connection network.
#' @param fun Function to apply in each focal point.
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco2)
#' myMatrix <- eco2[["P"]]
#' con <- eco.weight(XY=eco2[["XY"]], method="knearest", k=5)
#' result <- eco.slide.con(myMatrix, con, function(x)mean(x, na.rm = TRUE))
#'
#' image(matrix(myMatrix[,1], 30, 30)) # original image
#' image(matrix(result[,1], 30, 30)) # smoothed image
#' 
#' data(eco3)
#' myMatrix2 <- eco3[["P"]]
#' con <- eco.weight(XY=eco3[["XY"]], method="knearest", k=5)
#' plot(con)
#' # smoothing values in myMatrix2 using the connection network:
#' result <- eco.slide.con(myMatrix2, con, function(x)mean(x, na.rm = TRUE))
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
  
