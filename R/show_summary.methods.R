# Leandro Roser leandroroser@ege.fcen.uba.ar
# June 17, 2015 


###############################################################################
#                             SHOW METHODS
################################################################################


# show eco.gsa

setMethod("show", "eco.gsa", function(object) {
  method <- (object@METHOD)[1]
  cormethod <- (object@METHOD)[2]
  
  if(length(object@MULTI) != 0) {
    cat("\n", 
        "###########################","\n",
        paste(" ", method), "\n", 
        "###########################","\n\n")
    if(method == "Mantel test" | method == "Partial Mantel test") {
      cat("  Correlation coefficent used ->", cormethod, "\n")
    }
    cat(
        paste("  Number of simulations ->", object@NSIM), "\n",
        paste(" P-values correction method ->", object@ADJUST), "\n\n",
        paste(" Results:"), "\n\n")
    print(object@MULTI)
    cat("\n")
    
  } else {
    cat("\n", 
        "############################","\n",
        paste(" ", method), "\n", 
        "############################","\n\n")
    if(method == "Mantel test" | method == "Partial Mantel test") {
      cat("  Correlation coefficent used ->", cormethod, "\n")
    }
    cat(
        paste("  Number of simulations ->", object@NSIM), "\n",
        paste(" Alternative ->", object@ALTER), "\n", 
        paste(" P-value ->", object@PVAL), "\n", 
        paste(" Observed value ->", object@OBS), "\n",
        paste(" Expected value ->", object@EXP), "\n")
    cat("\n")
  }
})


#################################

# show eco.lsa

setMethod("show", "eco.lsa", function(object)  {
  cat("\n", 
      "#########################","\n",
      paste(" ", object@METHOD), "\n", 
      "#########################","\n\n",
      paste(" Test ->", object@TEST), "\n",
      paste(" Number of simulations ->", object@NSIM), "\n",
      paste(" Conditional ->", object@COND), "\n")
  if(object@TEST == "permutation") {
    cat(paste("  P-adjust method ->", object@PADJ))
  }
  cat("\n\n", 
      paste(" Results:"), "\n\n")
  print(object@OUT)
  cat("\n")
})

#################################

# show eco.correlog

setMethod("show", "eco.correlog", function(object) {
  
  randtest <- object@TEST
  method <- (object@METHOD)[1]
  cormethod <- (object@METHOD)[2]
  
  if(length(randtest != 0) & randtest != "none") {
  if(randtest == "permutation") {
    cat("\n", 
        "############################","\n",
        paste(" ", method), "\n", 
        "############################","\n\n")
    if(method == "Mantel statistic" | method == "Partial Mantel statistic") {
      cat("Correlation coefficent used ->", cormethod, "\n")
    }
    cat(paste(" Number of simulations ->", object@NSIM), "\n",
        paste(" Random test ->", object@TEST), "\n",
        paste(" P-adjust method ->", object@PADJUST), "\n\n", 
        paste(" Results:","\n\n"))
    print(object@OUT)
    cat("\n")
 
     } else if (randtest == "bootstrap") {
 
  cat("\n", 
      "############################","\n",
      paste(" ", method), "\n", 
      "############################","\n\n")
       if(method == "Mantel statistic" | method == "Partial Mantel statistic") {
         cat("Correlation coefficent used ->", cormethod, "\n")
       }
       cat(paste(" Number of simulations ->", object@NSIM), "\n",
      paste(" Random test ->", object@TEST), "\n",
      paste(" Results:","\n\n"))
  print(object@OUT)
       cat("\n")
     }
  } else {
    cat("\n", 
        "############################","\n",
        paste(" ", method), "\n", 
        "############################","\n\n")
    if(method == "Mantel statistic" | method == "Partial Mantel statistic") {
      cat("Correlation coefficent used ->", cormethod, "\n")
    }
    cat(paste(" Results:","\n\n"))
    print(object@OUT)
    cat("\n")
  }
  
  })


#################################

# show eco.weight

setMethod("show", "eco.weight", function(object)  {
  cat("\n", 
      "###################","\n",
      paste(" spatial weights"), "\n", 
      "###################","\n\n",
      paste(" Method ->", object@METHOD), "\n",
      " Parameters ->", paste("(", object@PAR, " = ",  
                              object@PAR.VAL, ")", sep =""), "\n",
      paste(" Row-standardization ->", object@ROW.SD), "\n")
  if(object@METHOD == "circle" | object@METHOD == "knearest") {
    cat(paste("  Self-included ->", object@SELF), "\n",
        paste(" Number of individuals ->", nrow(object@XY)), "\n",
        paste(" Non-zero (non-self) links ->", object@NONZERO, "%"), "\n",
        paste(" Individuals with non-zero (non-self) links ->", 
              object@NONZEROIND, "%"), "\n",
        paste(" Average (non-self) links per individual ->", object@AVG), "\n")
  } else {
    cat(paste("  Non-zero (non-self) links ->", object@NONZERO, "%"), "\n",
        paste(" Number of individuals ->", nrow(object@XY)), "\n",
        paste(" Individuals with non-zero (non-self) links ->", 
              object@NONZEROIND, "%"), "\n",
        paste(" Average (non-self) links per individual ->", object@AVG), "\n\n")
  }
})

#################################

# show eco.detrend

setMethod("show", "eco.detrend", function(object)  {
  cat("\n", 
      "###################","\n",
      paste(" ", "Data detrending"), "\n", 
      "###################","\n\n",
      paste(" Polynomial degree ->", object@POLY.DEG), "\n",
      paste(" Residuals - detrended data (@RES) ->"), "\n\n")
  print(object@RES)
  cat("\n\n", "Other data: 
      @XY -> projected coordinates
      @MODEL-> models
      @ANALYSIS -> eco.mlm object with details", 
      "\n\n")
})

#################################

# show eco.lagweight

setMethod("show", "eco.lagweight", function(object)  {
  cat("\n", 
      "#####################","\n",
      paste("spatial weights list"), "\n", 
      "#####################","\n\n",
      " Parameters ->", paste("(", object@PAR, " = ",  
                              round(object@PAR.VAL, 3), ")", sep =""), "\n",
      paste(" Row-standardization ->", object@ROW.SD), "\n",
      paste(" Self-included ->", object@SELF), "\n",
      paste(" Cummulative ->", object@CUMMUL), "\n",
      paste(" Number of classes ->", length(object@MEAN)), "\n",
      paste(" Method ->", object@METHOD), "\n\n")
  
})


#################################

# show eco.mlm

setMethod("show", 	"eco.mlm", function(object) {
  
  cat("\n\n", paste(" @MLM:", "multiple model results"), "\n", 
      paste(" @SUMMARY.MLM:", "summary of the results"), "\n", 
      paste(" @ANOVA.MLM:", "analysis of variance tables"), "\n",
      paste(" @PREDICTED:", "predicted values"), "\n",			
      paste(" @RESIDUALS:", "residuals of the analysis"), "\n\n")
})

#################################

# show eco.mctree

setMethod("show", 	"eco.mctree", function(object) {
  
  cat("\n\n", paste(" @TREES:", "trees"), "\n", 
      paste(" @CLASSPREDICT:", "predictions of the analysis"), "\n", 
      paste(" @FREQUENCIES:", "number of individuals predicted in each node"), "\n",
      paste(" @PREDICTED:", "predicted values"), "\n",			
      paste(" @RESIDUALS:", "residuals of the analysis"), "\n\n")
})



###############################################################################
#                             SUMMARY METHODS
################################################################################



#################################

# Summary for eco.lmtree output

setMethod("summary", 	"eco.mlm", 	function(object) {
  
  anova.pres <- function(mod) {
    
    asteriscos <- function(x, y) {
      
      if((x < 0.05) && (x >= 0.01)) {
        y <- paste(y, "*")
      } else if((x < 0.01) && (x >= 0.001)) {
        y <- paste(y, "**")
      } else if(x < 0.001) {
        y <- paste(y, "***")
      } else if(x >= 0.05) {
        y <- paste(y, "ns")
      }
    }
    
    
    final <- as.data.frame(matrix(ncol = length(object@MLM), 
                                  nrow = ncol(object@df2) + 1))
    rownames(final) <- c("Corrected model", colnames(object@df2))
    colnames(final) <- colnames(object@df1)
    pvalor <- rep(0, ncol(object@df1))
    
    for(i in 1:length(object@MLM)) { 
      pvalor[i] <- object@ANOVA.MLM[[i]][, 5][1]   
    }
    pvalor <- simplify2array(pvalor)
    pvalor <- p.adjust(pvalor, method = "fdr")
    
    
    for(i in 1:length(object@MLM)) {
      
      an <- object@ANOVA.MLM[[i]]
      an <- as.data.frame(an)
      mm <- object@SUMMARY.MLM[[i]]
      if(is.null(mm$fstatistic)) {
        final[1, i] <- "ns"
        
      } else {
        
        final[1, i] <- paste("F", "(", mm$fstatistic[2],
                             ", ", mm$fstatistic[3], ")","",
                             "=", round(mm$fstatistic[1],
                                        3), sep = "")
        final[1, i] <- asteriscos(pvalor[i], final[1, i])
        
      }
      
      for(j in 1:(nrow(an))) {
        
        solapa <- match(rownames(an), rownames(final))
        if(is.na(solapa[j]) == FALSE) {
          final[solapa[j], i] <- paste("F", "(", an[, 1][j],
                                       ", ", an[, 1][nrow(an)], ")",
                                       "= ", round(an[, 4][j], 3),
                                       sep= "")
          final[solapa[j], i] <- asteriscos(an[, 5][j], 
                                            final[solapa[j], i])
        }
      }
    }
    final[is.na(final) == TRUE] = "-"
    final
  }
  
  cat("\n", "ANOVA's F statistics, degrees of freedom and P-VALUES", "\n\n")
  print(anova.pres(object@MLM))
  cat("ns: non significative, *P<0.05, **P<0.01, ***P<0.001, ",
      "- null model")
  
})

#################################

#  summary,eco.mctree-method

setMethod("summary", 	"eco.mctree", 	function(object) {
  
  
  count <-0
  for(i in seq(along = object@TREES))
  {
    where.lev <- object@TREES[[i]]@where
    where.lev <- as.factor(where.lev)
    plot.tre<- max(levels(where.lev[[i]])) > 1
    if(plot.tre) {
      count <- count +1
    }
    if(plot.tre) {
      plot(object@TREES[[i]], main = names(object@TREES)[i])
    }
  }
  if(count == 0) {
    cat("\n", "There are not trees with significant splits to plot", "\n\n")
  }
  
})


####Show eco.IBD #####################################################

setMethod("show", "eco.IBD", function(object) {
  
  randtest <- (object@TEST)[1]
  metodo <- (object@METHOD)[2]
  
  cat("\n", 
      "###############","\n",
      " Fij analysis", "\n", 
      "###############","\n\n",
      paste(" Method ->", metodo), "\n",
      paste(" Number of simulations ->", object@NSIM), "\n")
      if(randtest == "permutation") {
      cat(paste("  P adjust method ->", object@PADJUST), "\n\n")
      } 
      if(metodo == "local") {
      cat(paste(" Conditional ->", (object@TEST)[2]), "\n\n")
        } 
      cat(
      "  Available information: ", "\n",
      paste(" @SP ->", "SP analysis"), "\n", 
      paste(" @OUT ->", "table with results"),"\n\n")
})



