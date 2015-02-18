# Creating input data for Geneland with an ecogen object
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015


setGeneric("eco.2geneland", 
           function(eco, ndig) {
             
             
             write.table(eco@XY, "XY.txt",
                         quote = FALSE,
                         row.names = FALSE,
                         col.names = FALSE)
             
             write.table(rownames(eco@XY),
                         "NAMES.txt",
                         quote = FALSE, 
                         row.names = FALSE,
                         col.names = FALSE)
             
             write.table(eco@P, 
                         "P.txt", 
                         quote = FALSE, 
                         row.names = FALSE,
                         col.names = FALSE)
             
             a <- eco@GENIND$type
             b <- eco@GENIND$ploidy
             if(a == "codom" && b != 1) {
               
               write.table(eco.2columns(eco, ndig),
                           "G.txt",quote = FALSE,
                           row.names = FALSE,
                           col.names = FALSE)
               
             } else if(a == "PA" | b == 1) {
               write.table(eco@G, 
                           "G.txt", 
                           quote = FALSE,
                           row.names = FALSE, 
                           col.names = FALSE)
             }
             
           })
