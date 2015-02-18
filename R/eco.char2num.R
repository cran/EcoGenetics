# Converting data coded with characters into data coded by numbers
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015

setGeneric("eco.char2num", 
           
           function(data, ndig = 0) {
             
             
             singlechar <- function(x) {
               y <- as.vector(as.matrix(x))
               y <- as.factor(y)
               original.code <- levels(y)
               y <- as.numeric(y)
               max.char <- max(nchar(y))
               y <- formatC(y, width=max.char, flag="0")
               y <- as.factor(y)
               new.code <- levels(y)[levels(y) != "NA"]
               y <- as.character(y)
               y[y == "NA"] <- paste(rep(0,max.char), collapse = "")
               y.tab <- data.frame(original.code, new.code)
               res <- list(y, y.tab)
               res
               
             }
             
             if(ndig == 0) {
               res <- singlechar(data)
               res[[1]] <- data.frame(matrix(res[[1]], ncol  = ncol(data),
                                             nrow = nrow(data)))
               for(i in 1:ncol(res[[1]])) {
                 res[[1]][,i] <- as.character(res[[1]][,i])
               }
               colnames(res[[1]]) <- colnames(data)
               rownames(res[[1]]) <- rownames(data)
               names(res) <- c("recoded_data", "conversion_table")
               res
               
             } else if(ndig != 0) {
               
               geno <- as.matrix(data)
               
               m1 <- substr(geno, 1, ndig)
               m2 <- substr(geno, ndig + 1, 2 * ndig)
               m2 <-c(m1, m2)
               m2[m2 == ""] <- NA
               recoding <- singlechar(m2)
               cadena <- recoding[[1]]
               m1 <- cadena[1:length(m1)]
               m2 <- cadena[(length(m1)+1): length(cadena)]
               m3 <- paste(m1, m2, sep ="")
               m3 <- matrix(m3, ncol  = ncol(data), nrow = nrow(data))
               m3 <-data.frame(m3)
               colnames(m3) <- colnames(data)
               rownames(m3) <- rownames(data)
               res <- list("recoded_data" = m3, "conversion_table" = recoding[[2]])
               for(i in 1:ncol(res[[1]])) {
                 res[[1]][,i]<-as.character(res[[1]][,i])
               }
               res
             }
             
             
           })

