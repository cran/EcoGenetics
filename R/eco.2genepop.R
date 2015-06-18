# Leandro Roser leandroroser@ege.fcen.uba.ar
# June 17, 2015 


# Exporting an ecogen genetic data frame into Genepop format

setGeneric("eco.2genepop", 
           function(eco, grp = NULL, ndig, 
                    name = "infile.genepop.txt") {
             
             
             grupo <- eco$S
             if(!is.null(grp)) {
             fact <- match(grp, colnames(eco$S), nomatch = -999)
             if(fact == -999) {
               stop("incorrect factor name")
             }
             structures <- as.factor(as.numeric(eco$S[, fact])) #recoding levels
             dat0 <- cbind(structures, eco$G)
             } else {
               dat0 <- cbind(rep(1, nrow(eco$XY)), eco$G)
             }
             
            
             dat0 <- as.matrix(dat0)
             rownames(dat0) <- rownames(eco$G)
             
             if(sum(is.na(dat0)) != 0) {
               dat0[is.na(dat0)] <- 0
             }
             
             datos <- data.frame(matrix(nrow = nrow(dat0), 
                                        ncol = ncol(dat0) - 1))
             
             n1 <- eco$GENIND$type
             n2 <- eco$GENIND$ploidy
             
             if((n1 == "codom") && (n2 != 1)) {
               
               if(is.na(ndig)) {
                 stop("incorrect ndig value")
               }
               a <- substr(dat0[, -1], 1, ndig)
               b<-substr(dat0[,-1],ndig+1, 2*ndig)
               a[a == " "] <- 0
               b[b == " "] <- 0
               
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
               a[a == "00"] <- "000"
               a[a == "0"] <- "000"
               b[b == "00"] <- "000"
               a[a == "0"] <- "000"
               a <- paste(a, b, sep = "")
             }
             
             a <- matrix(a, nrow = nrow(dat0), ncol = (ncol(dat0) - 1))
             a <- cbind(rep(0, nrow(dat0)), a)
             a[, 1] <- paste(rownames(dat0), ",")
             
             lista <- list()
             grp <- rep(" ", ncol(dat0))
             grp <- as.matrix(t(grp))
             grp[1] <- "POP"
             maxf <- max(as.numeric(dat0[, 1]))
             matriz <- matrix(nrow = 0, ncol = ncol(grp))
             for(i in 1:maxf) {
               lista[[i]] <- a[dat0[, 1] == i, ]
               lista[[i]] <- rbind(grp, lista[[i]])
               matriz <- rbind(matriz, lista[[i]])
             }
             matriz <- rbind(matriz, rep("", ncol(matriz)))
             
             matriz[, 1] <- as.character(matriz[, 1])
             nombres <- rep("", (ncol(dat0)) ^ 2)
             nombres <- as.data.frame(matrix(nombres, ncol(dat0), ncol(dat0)))
             nombres[, 1] <- c("Data exported from EcoGenetics", colnames(dat0[, -1]))
             colnames(nombres) <- colnames(matriz)
             matriz <- rbind(as.matrix(nombres), matriz) 

             write.table(matriz, name, row.names = FALSE,
                         col.names = FALSE, quote = FALSE)
             	
           })
