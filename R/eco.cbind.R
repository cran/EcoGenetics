# Combining the columns of two ecogen object

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015

setGeneric("eco.cbind", 
           function(e1, e2, ..., 
                    missing = 0) {
             
             
             z <- ecogen()
             
             
             u <- unlist(list(...))
             vec <- c("P", "G", "E", "S", "C", "ALL")
             m <- vec %in% u
             
             if((m[6] == TRUE) | !any(m)) {
               m <- rep(TRUE, 5)
             }
             
             
             z1 <- list(e1@P, e1@G, e1@E, e1@S, e1@C)
             z2 <- list(e2@P, e2@G, e2@E, e2@S, e2@C)
             
             tem <- list()
             
             for(i in 1:5) {
               
               if(m[[i]]) {
                 a <- nrow(z1[[i]])
                 b <- nrow(z2[[i]])
                 
                 if(any(a,b) == 0) 
                 {
                   if(a == 0 && b == 0) {
                     tem[[i]] <-data.frame()
                   } else if(a == 0 && b != 0) {
                     tem[[i]] <- z2[[i]]
                   } else if(a != 0 && b == 0) {
                     tem[[i]] <- z1[[i]]
                   }
                 } else {
                   if(any(rownames(z1[[i]]) != rownames(z2[[i]]))) {
                     stop(cat("Individuals in vec[[i]]  
                              data frame do not have the same names"))
                   }
                   tem[[i]] <-cbind(z1[[i]], z2[[i]])
                   }
               }
               }
             
             
             if(m[1] == TRUE) {
               
               z@P <- tem[[1]]
             }
             if(m[2] == TRUE) {
               
               
               z@G<-tem[[2]]
               
               if(all(dim(z@G)) != 0) {
                 type<-as.factor(as.vector(as.matrix(z@G)))
                 
                 if(length(levels(type)) != 2) {
                   
                   if(e1@GENIND$ploidy == 1) {
                     tempo <- df2genind(z@G, ploidy = 1)
                     
                   } else {
                     tempo <- df2genind(z@G)
                   } 
                 } else {
                   tempo <- df2genind(z@G, type = "PA")
                 }
               }
               
               
               z@GENIND$tab <- tempo$tab
               z@GENIND$ind.names <- tempo$ind.names
               z@GENIND$loc.names <- tempo$loc.names
               z@GENIND$loc.nall <- tempo$loc.nall
               z@GENIND$loc.fac <- tempo$loc.fac
               z@GENIND$all.names <- tempo$all.names
               z@GENIND$ploidy <- tempo$ploidy
               z@GENIND$type <- tempo$type
             }
             if(m[3] == TRUE) {
               
               z@E<-tem[[3]]
             }
             if(m[4] == TRUE) {
               
               z@S <-tem[[4]]
             }
             if(m[5] == TRUE) {
               
               z@C <- tem[[5]]
               
             }
             
             if(nrow(e1@XY) != 0) {
               z@XY <- e1@XY
             } else if(nrow(e2@XY) != 0) {
               z@XY <- e2@XY
             } 
             
             attr(z, "format") <- attr(e1, "format")
             attr(z, "type") <-  attr(e1, "type")
             attr(z, "missing") <- attr(e1, "missing")
             attr(z, "ploidy") <- attr(e1, "ploidy")
             
             return(z)
             })
