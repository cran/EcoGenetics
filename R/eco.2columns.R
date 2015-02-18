# Converting an ecogen genetic data frame with two alleles per column
# into a data frame with one allele per column 
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015

setGeneric("eco.2columns", 
           function(eco, ndig) {
             
             
             geno <- as.matrix(eco@G)
             
             m1 <- substr(geno, 1, ndig)
             
             m1 <- matrix(m1, ncol = ncol(geno),
                          nrow = nrow(geno)) 
             
             m2 <- substr(geno, ndig + 1, 2 * ndig)
             
             m2 <- matrix(m2, ncol = ncol(geno),
                          nrow = nrow(geno))
             
             m3 <- matrix(, ncol = 2 * ncol(geno),
                          nrow = nrow(geno))
             
             colnames(m3) <- c(1:ncol(m3))
             
             j <- 1
             for(i in 1:ncol(m3)) {
               
               if((i %% 2) == 0) {
                 
                 m3[, i] <- m2[, j]
                 j <- i/2
                 colnames(m3)[i] <- paste(colnames(geno)[j], ".2", sep = "")
               } else {
                 if(i != 1) {
                   j <- j + 1
                 }
                 m3[, i] <- m1[, j]
                 colnames(m3)[i] <- paste(colnames(geno)[j], ".1", sep = "")
               }
             }
             
             m3[m3 == "0"] <- "NA"
             m3[m3 == ""] <- "NA"
             m3 <- as.data.frame(m3)
             
             m3
             
           })
