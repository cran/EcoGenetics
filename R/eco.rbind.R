# Leandro Roser leandroroser@ege.fcen.uba.ar
# June 17, 2015 


# Combining the rows of two ecogen objects

setGeneric("eco.rbind", 
           function(e1, e2)  {
             
             nom <- c(rownames(e1$G), rownames(e2$G))
             
             if(any(duplicated(nom))) {
               stop("duplicated row names are not allowed")
             }
                 
             z <- ecogen(G = rbind(e1$G, e2$G))
             
             
             if(all(dim(z$G)) != 0) {
               type<-as.factor(as.vector(as.matrix(z$G)))
               
               if(length(levels(type)) != 2) {
                 if(e1$GENIND$ploidy == 1) {
                   
                   tempo <- df2genind(z$G, ploidy = 1)
                   
                 } else {
                   tempo <- df2genind(z$G)
                 } 
               } else {
                 tempo <- df2genind(z$G, type = "PA")
               }
             }
             
             z$GENIND$tab <- tempo$tab
             z$GENIND$ind.names <- tempo$ind.names
             z$GENIND$loc.names <- tempo$loc.names
             z$GENIND$loc.nall <- tempo$loc.nall
             z$GENIND$loc.fac <- tempo$loc.fac
             z$GENIND$all.names <- tempo$all.names
             z$GENIND$ploidy <- tempo$ploidy
             z$GENIND$type <- tempo$type
             
             z$XY <- rbind(e1$XY, e2$XY)
             z$P <- rbind(e1$P, e2$P)   
             z$E <- rbind(e1$E, e2$E)
             z$S <- rbind(e1$S, e2$S)
             z$C <- rbind(e1$C, e2$C)
             z$OUT <- list()
             
             attr(z, "format") <- attr(e1, "format")
             attr(z, "type") <-  attr(e1, "type")
             attr(z, "missing") <- attr(e1, "missing")
             attr(z, "ploidy") <- attr(e1, "ploidy")
             
             return(z)
           })
