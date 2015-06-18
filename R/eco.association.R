# Leandro Roser leandroroser@ege.fcen.uba.ar
# June 17, 2015 


# Chi-square and Fisher's exact test for association of loci and alleles 
# with a factor

setGeneric("eco.association",  
           function(eco,
                    assoc = c("within", "between"),
                    x, 
                    method = c("fisher.test", "chisq.test"), 
                    nrep = 99,
                    adjust = "none", 
                    ndig = NA) {
             
             
             cat("\n","#### Test of association <", assoc,">, with", nrep,
                 "repetitions","####", "\n")
             
             assoc  <- match.arg(assoc)
             method <- match.arg(method)
             
             
               
               cat(paste("\n","Method: ", method,"\n",
                         "P-adjust method: ", adjust, "\n\n"))
           
             fact <- match(x, colnames(eco@S), nomatch = 0)
             fact <- fact[fact != 0]
             if(length(fact) == 0) {
               stop("incorrect factor name")
             }
             
             
             if(assoc == "within") {
               
               grp <- eco@S[,fact]
               met <- c("chisq.test", "fisher.test")
               nom <- match(method,met, nomatch = 0)
               nom <- nom[nom != 0]
               if(length(nom) == 0) {
                 stop(paste("invalid", method))
               }
               marcadores <- eco@GENIND$tab
               marcadores[marcadores == 1.0] <- 2
               marcadores[marcadores == 0.5] <- 1
               ifelse(method == "chisq.test", 
                      chi <- apply(marcadores, 2, 
                                   function(u){ 
                                     chisq.test(u, grp, simulate.p.value = T, 
                                                B = nrep)
                                   }),
                      chi <- apply(marcadores, 2,  
                                   function(u){
                                     fisher.test(u, 	grp, simulate.p.value = T,
                                                 B = nrep)
                                   }))
               
               if(method == "chisq.test") {
                 estadistico <- sapply(chi, function(u) c(u$statistic))
                 estadistico<- round((data.frame(estadistico)), 3)
                 pvalor <- sapply(chi, function(u) c(u$p.value))
                 pvalor <- round(data.frame(p.adjust(pvalor, method = adjust)), 3)
                 pvalor <- data.frame(estadistico, pvalor)
                 colnames(pvalor) <- c("statistic", "P value")
                 rownames(pvalor) <- colnames(marcadores)
                 return(pvalor)
               } else if(method == "fisher.test") {
                 
                 pvalor <- sapply(chi, function(u) c(u$p.value))
                 pvalor <- round(data.frame(p.adjust(pvalor, method = adjust)), 3)
                 colnames(pvalor) <- "P value"
                 rownames(pvalor) <- colnames(marcadores)
                 return(pvalor)
               }
               
             } else if(assoc == "between") {
               
               if(is.na(ndig)) {
                 stop("please provide a value for the number digits per allele (ndig)")
               }
               
               x <- aue.sort(eco@G, ndig)
               x <-as.matrix(x)
               x[x == 0] <-NA
               pop <- eco@S[, fact]
               nloc <- ncol(x)
               y <-list()
               for(i in 1:nloc) {
                 y[[i]] <-as.factor(x[,i])
               }
               
               tab <- list()
               for(i in 1:nloc) {
                 tab[[i]] <- table(pop,y[[i]], exclude =NA)[,-1]
               }
               
               chi <- data.frame(matrix(nrow = 2, ncol = nloc))
               rownames(chi) <- c("statistic", "P value")
               colnames(chi) <- colnames(x)
               
               
               
               if(method == "chisq.test") { 
                 chi <- data.frame(matrix(nrow = 2, ncol = nloc))
                 rownames(chi) <- c("statistic", "P value")
                 colnames(chi) <- colnames(x)
                 
                 for(i in 1:nloc) {
                   temp <- chisq.test(tab[[i]], 
                                      simulate.p.value =TRUE,
                                      B = nrep)
                   chi[, i] <- c(round(temp$statistic, 3), 
                                 round(temp$p.value, 3))
                   
                 }
                 if(adjust != "none") {
                   chi[2,] <- round(p.adjust(chi[2,], method = adjust), 3)
                   return(t(chi))
                   
                 } else {
                   return(t(chi))
                 }
               } else if(method == "fisher.test")  {
                 chi <- data.frame(matrix(nrow = 1, ncol = nloc))
                 rownames(chi) <- "P value"
                 colnames(chi) <- colnames(eco@G)
                 
                 
                 for(i in 1:nloc) {
                   temp <- fisher.test(tab[[i]], 
                                       simulate.p.value = TRUE,
                                       B = nrep)
                   chi[, i] <- temp$p.value
                 }
                 
                 if(adjust != "none") {
                   chi<- round(p.adjust(chi, method = adjust), 3)
                   return(t(chi))
                   
                 } else {
                   chi <- round(chi,3)
                   return(t(chi))
                 }
               }
               
             }
             
           })
