# Clearing the working environment, maintaining only the specified objects
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015

eco.clear <- function(...) {
  
  clean.names <- as.character(match.call())[-1]
  
  env <- parent.frame()
  cuales <- ls(envir = env) 
  cuales  <- cuales %in% clean.names
  rm(list = ls(envir = env)[!cuales], envir = env)
  
}
