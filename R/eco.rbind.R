# Combining the rows of two ecogen objects
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015

setGeneric("eco.rbind", 
           function(e1, e2)  {
             
             
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
             
             return(z)
           })
