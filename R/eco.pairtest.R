# Kruskall - Wallis + Wilcoxon (Mann-Whitney U) and aov + Tukey-HSD tests 
# for an ecogen object

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015

setGeneric("eco.pairtest", 
           
           function(eco, df = c("P", "E","GENIND", "C"),
                    x, test = c("wilcoxon", "tukey"),
                    adjust = "fdr",
                    only.p = TRUE, ...) {
             
             test <- match.arg(test)
             
             grupo <- eco@S
             fact <- match(x, colnames(eco@S), nomatch = 0)
             fact <- fact[!fact == 0]
             if(length(fact) == 0) {
               stop("incorrect factor name")
             }
             
             grupos <- as.factor(eco@S[, fact])
             
             P <- match.arg(df)
             
             if(P == "P") {
               P<-eco$P
             } else if(P == "E") {
               P<-eco$E
             } else if(P == "GENIND") {
               P<- eco$GENIND$tab
             } else if(P == "C") {
               P<-eco$C
             }
             
             if(test == "wilcoxon") {
               
               
               niveles <- as.numeric(max(levels(grupos)))
               lev <- list()
               
               tabla <- table(1:length(grupos), grupos)
               
               a <- matrix(0, nrow = niveles, ncol = niveles)
               a <- upper.tri(a)
               index <- which(a == TRUE, arr.ind = TRUE)
               index <- data.frame(index, rep(0, nrow(index)), rep(0, nrow(index)))
               colnames(index) <- c("pop1", "pop2", "est", "P")
               res <- list()
               krus <- list()
               for(i in 1:ncol(P)) {
                 tabla <- index
                 temp <- kruskal.test(P[,i]~grupos)
                 krus[[i]] <- c(temp$statistic, temp$p.value)
                 for(j in 1:nrow(index)) {
                   wil <- wilcox.test(P[grupos == index[j, 1], i], 
                                      P[grupos == index[j, 2],i], ...)
                   tabla[j, 3] <- wil$statistic
                   tabla[j, 4] <- round(p.adjust(wil$p.value, method = adjust), 4)
                 }
                 res[[i]] <- tabla
               }
               krus <- sapply(krus, c)
               krus <- round(krus, 4)
               colnames(krus) <- colnames(P)
               rownames(krus) <- c("statistic", "P")
               names(res) <- names(P)
               
               if(only.p == TRUE) {
                 pvalues <- sapply(res, function(y) {y[, 4]})
                 colnames(pvalues) <- colnames(P)
                 rownames(pvalues) <- paste(index[, 1], index[, 2], sep = "-")
                 return(list(kruskall.test = krus, wilcoxon.test = pvalues, 
                             p.correction = adjust))
               } else {
                 return(list(kruskall.test = krus, wilcoxon.test = res, 
                             p.correction = adjust))
               } 
               
             } else if(test == "tukey") {
               
               m.tukey <- function(x, y) {
                 general <- listafac <- tukeyfac <- list()
                 for(i in 1:ncol(x)) {
                   listafac[[i]] <- aov(x[, i] ~ grupos)
                   temporal <- anova(listafac[[i]])
                   general[[i]] <- c(temporal[,4][1], temporal[, 5][1])
                   tukeyfac[[i]] <- post <- TukeyHSD(listafac[[i]], ...)
                 }
                 general <- sapply(general, c)
                 rownames(general) <- c("F", "P")
                 colnames(general) <- colnames(P)
                 names(tukeyfac) <- names(P)
                 res <- list(general, tukeyfac)
               }
               resultados <- m.tukey(P, grupos)
               tuk <- resultados[[2]]
               general <- round(resultados[[1]], 4)
               p.values <- sapply(tuk, function(y) round(y$grupos[, 4], 4))
               rownames(p.values) <- rownames(tuk[[1]]$grupos)
               colnames(p.values) <- colnames(P)
               
               if(only.p == FALSE) {
                 return(list(aov = general, pairtest = tuk))
               } else  {
                 return(list(aov = general, pairtest = p.values))
               }
               
             }
           })
