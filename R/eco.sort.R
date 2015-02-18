# Ordering cells in a data frame or matrix
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015


setGeneric("eco.sort",
           function(x, ndig)  {
             
             geno <- as.matrix(x)
             
             m1 <- substr(geno, 1, ndig)
             m2 <- substr(geno, ndig + 1, 2 * ndig)
             m3 <- matrix(, ncol = ncol(geno), nrow = nrow(geno))
             m4 <- matrix(, ncol = ncol(geno), nrow = nrow(geno))
             
             m3[which(m1 < m2)] <- m1[which(m1 < m2)]
             m3[which(m2 < m1)] <- m2[which(m2 < m1)]
             m3[which(m2 == m1)] <- m2[which(m2 == m1)]
             
             m4[which(m1 > m2)] <- m1[which(m1 > m2)]
             m4[which(m2 > m1)] <- m2[which(m2 > m1)]
             m4[which(m2 == m1)] <- m2[which(m2 == m1)]
             
             m5 <- paste(m3, m4, sep = "")
             m5 <- as.data.frame(matrix(m5, ncol = ncol(geno), nrow = nrow(geno)))
             m5
             
           })
