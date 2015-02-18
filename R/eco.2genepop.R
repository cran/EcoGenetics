# Exporting an ecogen genetic data frame into Genepop format
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015. 


setGeneric("eco.2genepop", 
           function(eco, x, ndig) {
             
             
             grupo <- eco$S
             fact <- match(x, colnames(eco$S), nomatch = 0)
             fact <- fact[!fact == 0]
             if(length(fact) == 0) {
               stop("incorrect factor name")
             }
             
             dat0 <- cbind(eco$S[, fact], eco$G)
             dat0 <- as.matrix(dat0)
             rownames(dat0) <- rownames(eco$G)
             
             if(sum(is.na(dat0)) != 0) {
               dat0[is.na(dat0)] <- 0
             }
             
             datos <- data.frame(matrix(nrow = nrow(dat0), ncol = ncol(dat0) - 1))
             
             n1 <- eco$GENIND$type
             n2 <- eco$GENIND$ploidy
             
             if((n1 == "codom") && (n2 != 1)) {
               
               if(is.na(ndig)) {
                 stop("incorrect ndig value")
               }
               a <- substr(dat0[, -1], 1, ndig)
               b<-substr(dat0[,-1],ndig+1,2*ndig)
               
             } else if((n1 == "PA") | (n2 > 1)) { 
               a <- as.matrix(eco$G)
               a[is.na(a)] <- 0
               
             } else if (n2 == 1) {
               a<-as.numeric(as.factor(as.matrix(eco$G)))
               a<-matrix(a, ncol = ncol(eco$G), nrow = nrow(eco$G))
               a[is.na(a)] <- 0
             }
             
             if(ndig == 1) {
               a <- paste("00", a, sep = "")
               b <- paste("00", b, sep = "")
               
             } else if(ndig == 2) {
               a <- paste("0", a, sep = "")
               b <- paste("0", b, sep = "")
               
             }
             
             if((n1 == "codom") && (n2 != 1)) {
               a[a == "00 "] <- "000"
               b[b == "00 "] <- "000"
               a <- paste(a, b, sep = "")
             }
             
             a <- matrix(a, nrow = nrow(dat0), ncol = (ncol(dat0) - 1))
             a <- cbind(rep(0, nrow(dat0)), a)
             a[, 1] <- paste(rownames(dat0), ",")
             
             lista <- list()
             pop <- rep(" ", ncol(dat0))
             pop <- as.matrix(t(pop))
             pop[1] <- "POP"
             maxf <- max(as.numeric(dat0[, 1]))
             matriz <- matrix(nrow = 0, ncol = ncol(pop))
             for(i in 1:maxf) {
               lista[[i]] <- a[dat0[, 1] == i, ]
               lista[[i]] <- rbind(pop, lista[[i]])
               matriz <- rbind(matriz, lista[[i]])
             }
             
             matriz[, 1] <- as.character(matriz[, 1])
             nombres <- rep("", (ncol(dat0)) ^ 2)
             nombres <- as.data.frame(matrix(nombres, ncol(dat0), ncol(dat0)))
             nombres[, 1] <- c("Data exported from EcoGenetics", colnames(dat0[, -1]))
             colnames(nombres) <- colnames(matriz)
             matriz <- rbind(as.matrix(nombres), matriz) 
             write.table(matriz, "infile.genepop.txt", row.names = FALSE,
                         col.names = FALSE, quote = FALSE)
           })
