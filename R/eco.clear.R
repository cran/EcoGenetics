# Leandro Roser leandroroser@ege.fcen.uba.ar
# June 17, 2015 


# Clearing the working environment, maintaining only the specified objects


eco.clear <- function(..., all = FALSE) {
  
  clean.names <- as.character(match.call())[-1]
  
  env <- parent.frame()
  cuales <- ls(envir = env, all.names = all) 
  cuales  <- cuales %in% clean.names
  rm(list = ls(envir = env, all.names = all)[!cuales], envir = env)
  
}
