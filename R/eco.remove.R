# Creating an updated ecogen object by removing results of the slot OUT
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015

setGeneric("eco.remove", 
           
           function(object, ...) {
             
             res.names <- as.character(match.call())
             res.names <- res.names[-c(1:2)]
             del <- (names(object@OUT)) %in% res.names
             object@OUT <- object@OUT[!del]
             object
           })
