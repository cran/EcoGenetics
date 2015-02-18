# Merging two ecogen objects. 
# Ordering the rows of an ecogen object according to the rows of another
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015


setGeneric("eco.merge",
           function(e1, e2, ...) {
             
             
             u <- unlist(list(...))
             vec <- c("P", "G", "E", "S", "C", "ALL")
             m <- vec %in% u
             
             if(!any(m)) {
               m <- rep(TRUE, 6)
             }
             
             
             if(m[6] == TRUE) {
               m <- rep(TRUE, 6)
             }
             
             z <- new("ecogen")
             
             if(all(dim(e1@XY) != 0) && 
                  all(dim(e2@XY) != 0)) {
               z@XY <- merge(data.frame(rownames(e1@XY),
                                        c(1:nrow(e1@XY)), e1@XY),
                             data.frame(rownames(e2@XY), e2@XY), by = 1)
               z@XY <- z@XY[order(z@XY[, 2]), ]
               rownames(z@XY) <- z@XY[, 1]
               z@XY <- z@XY[, -c(1:4)]
             }
             
             
             if((m[1] == TRUE) & all(dim(e1@P) != 0) & 
                  all(dim(e2@P) != 0)) {
               z@P <-  merge(data.frame(rownames(e1@P),
                                        c(1:nrow(e1@P)), e1@P),
                             data.frame(rownames(e2@P), e2@P), by = 1)
               z@P <- z@P[order(z@P[, 2]), ]
               rownames(z@P) <- z@P[, 1]
               z@P <- z@P[, -c(1, 2)]
             }
             
             
             if((m[2] == TRUE) & all(dim(e1@G) != 0) &
                  all(dim(e2@G) != 0)) {
               z@G <- merge(data.frame(rownames(e1@G), c(1:nrow(e1@G)),
                                       e1@G), data.frame(rownames(e2@G),
                                                         e2@G),by = 1)
               z@G <- z@G[order(z@G[, 2]), ]
               rownames(z@G) <- z@G[, 1]
               z@G <- z@G[, -c(1, 2)]
               
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
             
             
             if((m[3] == TRUE) & all(dim(e1@E) != 0) &
                  all(dim(e2@E) != 0)) {
               z@E <- merge(data.frame(rownames(e1@E), c(1:nrow(e1@E)),
                                       e1@E), data.frame(rownames(e2@E),
                                                         e2@E), by = 1)
               z@E <- z@E[order(z@E[, 2]), ]
               rownames(z@E) <- z@E[, 1]
               z@E <- z@E[, -c(1,2)]
             }
             
             if((m[4] == TRUE) & all(dim(e1@S) != 0) &
                  all(dim(e2@S) != 0)) {
               z@S <- merge(data.frame(rownames(e1@S),
                                       c(1:nrow(e1@S)), e1@S), 
                            data.frame(rownames(e2@S), e2@S), by = 1)
               z@S <- z@S[order(z@S[, 2]), ]
               rownames(z@S) <- z@S[, 1]
               z@S <- z@S[, -c(1, 2)]
             }
             
             if((m[5] == TRUE) & all(dim(e1@C) != 0) &
                  all(dim(e2@C) != 0)) {
               z@C <- merge(data.frame(rownames(e1@C),
                                       c(1:nrow(e1@C)), 
                                       e1@C),
                            data.frame(rownames(e2@C), 
                                       e2@C), by=1)
               z@C <- z@C[order(z@C[, 2]), ]
               rownames(z@C)<- z@C[, 1]
               z@C <- z@C[, -c(1,2)]
             }
             
             z
             
             
           } )
