# Converting an ecogen genetic data frame into a hierfstat data frame
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015

setGeneric("eco.2hierfstat", 
           function(eco, fact = NULL) {
             
             u <- eco$G
             
             grupo <- eco@S
             
             if(is.null(fact))
             {
               factord <- 	as.data.frame(rep(1, nrow(u)))
               cnom <- "pop"
               rnom <- rownames(eco@G)
               Gord <- u
             } else {
               
               fact <- match(fact, colnames(eco@S), nomatch = 0)
               fact <- fact[fact != 0]
               if(length(fact) == 0) {
                 stop("incorrect factor name")
               }
               orden <- order(eco@S[, fact])
               Gord <- u[orden,]
               factord <- eco@S[orden, fact]
               factord <- as.numeric(factord)
               cnom <- colnames(eco@S[fact])
               rnom <- rownames(eco@G)[orden]
             }
             
             datahier <- data.frame(factord, Gord)
             colnames(datahier)[1] <- cnom
             rownames(datahier) <- rnom
             datahier
             
           })
