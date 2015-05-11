# Clearing the working environment, maintaining only the specified objects

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015

eco.clear <- function(..., all = FALSE) {
  
  clean.names <- as.character(match.call())[-1]
  
  env <- parent.frame()
  cuales <- ls(envir = env, all.names = all) 
  cuales  <- cuales %in% clean.names
  rm(list = ls(envir = env, all.names = all)[!cuales], envir = env)
  
}
