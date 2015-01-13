#' Clearing the working environment, maintaining only the specified objects.
#' @param ... Objects to retain.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @examples
#' \dontrun{
#' 
#' data(eco.test)
#' ls()
#' eco.clear(eco)
#' ls()
#' 
#' }
#' 
#' @export


eco.clear <- function(...) {
  
clean.names <- as.character(match.call())[-1]

env <- parent.frame()
cuales <- ls(envir = env) 
cuales  <- cuales %in% clean.names
rm(list = ls(envir = env)[!cuales], envir = env)

}
