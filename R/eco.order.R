# Leandro Roser leandroroser@ege.fcen.uba.ar
# June 17, 2015 


# Ordering the rows of the data frames contained in an ecogen object

setGeneric("eco.order", 
           function(eco) {
             
             
             posit0 <- seq(along=eco@XY[,1])
             
             posit <- seq(along = rownames(eco@P))
             temporalP <- merge(data.frame(rownames(eco@XY),
                                           posit0), 
                                data.frame(rownames(eco@P),
                                           rownames(eco@P),
                                           posit), 
                                by = 1)
             temporalP <- temporalP[order(temporalP$posit0), ]
             c.names <- colnames(eco@P)
             r.names <- temporalP[,3]
             eco@P <- data.frame(eco@P[temporalP$posit, ])
             colnames(eco@P) <- c.names
             rownames(eco@P) <- r.names  
             
             posit <- seq(along = rownames(eco@G))
             temporalG <- merge(data.frame(rownames(eco@XY),
                                           posit0), 
                                data.frame(rownames(eco@G), 
                                           rownames(eco@G),
                                           posit),
                                by = 1)
             temporalG <- temporalG[order(temporalG$posit0), ]
             c.names <- colnames(eco@G)
             r.names <- temporalG[,3]
             eco@G <- data.frame(eco@G[temporalG$posit, ])
             colnames(eco@G) <- c.names
             rownames(eco@G) <- r.names
             
             posit <- seq(along = rownames(eco@E))
             temporalE <- merge(data.frame(rownames(eco@XY),
                                           posit0),
                                data.frame(rownames(eco@E), 
                                           rownames(eco@E),
                                           posit),
                                by = 1)
             temporalE <- temporalE[order(temporalE$posit0), ]
             c.names <- colnames(eco@E)
             r.names <- temporalE[,3]
             eco@E <- data.frame(eco@E[temporalE$posit, ])
             colnames(eco@E) <- c.names
             rownames(eco@E) <- r.names
             
             
             posit <- seq(along = rownames(eco@S))
             temporalS <- merge(data.frame(rownames(eco@XY),
                                           posit0),
                                data.frame(rownames(eco@S), 
                                           rownames(eco@S), 
                                           posit),
                                by = 1)
             temporalS <- temporalS[order(temporalS$posit0), ]
             c.names <- colnames(eco@S)
             r.names <- temporalS[,3]
             eco@S <- data.frame(eco@S[temporalS$posit, ])
             colnames(eco@S) <- c.names
             rownames(eco@S) <- r.names
             
             
             posit <- seq(along = rownames(eco@C))
             temporalC <- merge(data.frame(rownames(eco@XY), 
                                           posit0), 
                                data.frame(rownames(eco@C),
                                           rownames(eco@C),
                                           posit),
                                by = 1)
             temporalC <- temporalC[order(temporalC$posit0), ]
             c.names <- colnames(eco@C)
             r.names <- temporalC[,3]
             eco@C <- data.frame(eco@C[temporalC$posit, ])
             colnames(eco@C) <- c.names
             rownames(eco@C) <- r.names
          
             eco
             
           })
