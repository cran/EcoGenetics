# Creating an updated ecogen object by removing
# results of the slot @OUT

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015

setGeneric("eco.remove", 
           
           function(eco, ...) {
             
             res.names <- as.character(match.call())
             res.names <- res.names[-c(1:2)]
             del <- (names(eco@OUT)) %in% res.names
             eco@OUT <- eco@OUT[!del]
             eco
           })
