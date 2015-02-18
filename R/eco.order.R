# Ordering the rows of the data frames contained in an ecogen object
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015

setGeneric("eco.order", 
           function(object) {
             
             
             posit0 <- seq(along=object@XY[,1])
             
             posit <- seq(along = rownames(object@P))
             temporalP <- merge(data.frame(rownames(object@XY),
                                           posit0), 
                                data.frame(rownames(object@P),
                                           rownames(object@P),
                                           posit), 
                                by = 1)
             temporalP <- temporalP[order(temporalP$posit0), ]
             c.names <- colnames(object@P)
             r.names <- temporalP[,3]
             object@P <- data.frame(object@P[temporalP$posit, ])
             colnames(object@P) <- c.names
             rownames(object@P) <- r.names  
             
             posit <- seq(along = rownames(object@G))
             temporalG <- merge(data.frame(rownames(object@XY),
                                           posit0), 
                                data.frame(rownames(object@G), 
                                           rownames(object@G),
                                           posit),
                                by = 1)
             temporalG <- temporalG[order(temporalG$posit0), ]
             c.names <- colnames(object@G)
             r.names <- temporalG[,3]
             object@G <- data.frame(object@G[temporalG$posit, ])
             colnames(object@G) <- c.names
             rownames(object@G) <- r.names
             
             posit <- seq(along = rownames(object@E))
             temporalE <- merge(data.frame(rownames(object@XY),
                                           posit0),
                                data.frame(rownames(object@E), 
                                           rownames(object@E),
                                           posit),
                                by = 1)
             temporalE <- temporalE[order(temporalE$posit0), ]
             c.names <- colnames(object@E)
             r.names <- temporalE[,3]
             object@E <- data.frame(object@E[temporalE$posit, ])
             colnames(object@E) <- c.names
             rownames(object@E) <- r.names
             
             
             posit <- seq(along = rownames(object@S))
             temporalS <- merge(data.frame(rownames(object@XY),
                                           posit0),
                                data.frame(rownames(object@S), 
                                           rownames(object@S), 
                                           posit),
                                by = 1)
             temporalS <- temporalS[order(temporalS$posit0), ]
             c.names <- colnames(object@S)
             r.names <- temporalS[,3]
             object@S <- data.frame(object@S[temporalS$posit, ])
             colnames(object@S) <- c.names
             rownames(object@S) <- r.names
             
             
             posit <- seq(along = rownames(object@C))
             temporalC <- merge(data.frame(rownames(object@XY), 
                                           posit0), 
                                data.frame(rownames(object@C),
                                           rownames(object@C),
                                           posit),
                                by = 1)
             temporalC <- temporalC[order(temporalC$posit0), ]
             c.names <- colnames(object@C)
             r.names <- temporalC[,3]
             object@C <- data.frame(object@C[temporalC$posit, ])
             colnames(object@C) <- c.names
             rownames(object@C) <- r.names
             
             
             object
             
           })
