#' Wilcoxon (Mann-Whitney U) and Tukey-HSD tests for an ecogen object.
#' @param eco Ecogen object.
#' @param df The data frame for the analysis. Could be "P", "E" or "C".
#' @param x The name of the S slot column with the groups for the analysis.
#' @param test Test to perform ("wilcoxon", "tukey").
#' @param  only.p  Should be only returned a matrix with P values? 
#' Default TRUE.
#' @param adjust The p.value adjustment method passed to \code{\link{p.adjust}}. 
#' for the Wilkoxon tests. Defalut is "fdr".
#' @param ... Additional arguments passed to \code{\link{wilcox.test}} 
#' or  \code{\link{TukeyHSD}}.
#' @details This program returns the Wilcoxon (Mann-Whitney U) or Tukey-HSD 
#' statistics and p values for the multiple comparisons of the variables contained
#' in the selected data frame,  among the levels of a factor of the slot "S".
#' @seealso \code{\link{wilcox.test}} \code{\link{TukeyHSD}}
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @examples 
#' \dontrun{
#' 
#' data(eco.test)
#' wil <- eco.pairtest(eco = eco, df = "P", x = "structure")
#' wil <- eco.pairtest(eco = eco,df = "E", x = "structure")
#' wil <- eco.pairtest(eco = eco, df = "P", x = "structure", only.p = FALSE)
#' wil <- eco.pairtest(eco = eco,df = "P", x = "structure", test = "tukey")
#' 
#' }
#' @export

setGeneric("eco.pairtest", 
					 
					 function(eco, df = c("P", "E", "C"),
					 				 x, test = c("wilcoxon", "tukey"),
					 				 adjust = c("fdr", "holm", "hochberg", 
					 				 					 "hommel", "bonferroni", "BH",
					 				 					 "BY", "none"),
					 				 only.p = TRUE, ...) {
					 	
	test <- match.arg(test)
	adjust <- match.arg(adjust)

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
	}
	else if(P == "C") {
		P<-eco$C
	}
  
  if(test == "wilcoxon") {
    
    
    niveles <- as.numeric(max(levels(grupos)))
    lev <- list()
    
    tabla <- table(1:length(grupos), grupos)
    
    a <- matrix(0, nrow = niveles, ncol = niveles)
    a <- upper.tri(a)
    index <- which(a == TRUE, arr.ind = TRUE)
    index <- data.frame(index, rep(0, 3), rep(0, 3))
    colnames(index) <- c("pop1", "pop2", "est", "P")
    res <- list()
    for(i in 1:ncol(P)) {
      tabla <- index
      for(j in 1:nrow(index)) {
        wil <- wilcox.test(P[grupos == index[j, 1], i], 
                           P[grupos == index[j, 2],i], ...)
        tabla[j, 3] <- wil$statistic
        tabla[j, 4] <- round(p.adjust(wil$p.value, method = adjust), 4)
      }
      res[[i]] <- tabla
    }
    names(res) <- names(P)
    
    if(only.p == TRUE) {
    pvalues <- sapply(res, function(y) {y[, 4]})
    colnames(pvalues) <- names(P)
    rownames(pvalues) <- paste(index[, 1], index[, 2], sep = "-")
    return(pvalues)
    } else if(only.p == FALSE) {
      return(res)
    } 
    
  } else if(test == "tukey") {
    
    m.tukey <- function(x, y) {
      listafac <- tukeyfac <- list()
      for(i in 1:ncol(x)) {
        listafac[[i]] <- aov(x[[i]] ~ grupos)
        tukeyfac[[i]] <- post <- TukeyHSD(listafac[[i]], ...)
      }
      names(tukeyfac) <- names(P)
      tukeyfac
    }
    tuk <- m.tukey(P, grupos)
    
    p.values <- sapply(tuk, function(y) round(y$grupos[, 4], 4))
    rownames(p.values) <- rownames(tuk[[1]]$grupos)
    
    if(only.p == FALSE) {
      return(tuk)
    }
    return(p.values)
    
  }
})
