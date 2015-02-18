########--RelDist-Class--#########################################################
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015

setClass("RelDist", 
         
         representation(OUT = "list",
                        INTERVAL = "numeric", 
                        MAX = "numeric",
                        TYPE = "character", 
                        NAMES = "character",
                        METHOD = "character", 
                        RANDTEST = "character"), 
         
         prototype( OUT = list(), 
                    INTERVAL = 0,
                    MAX = 0, 
                    TYPE = character(), 
                    NAMES = character(), 
                    METHOD = character(),
                    RANDTEST = character())
)

##################################################################################
# Other classes definitions

setClass("eco.multiboot", contains = "RelDist")
setClass("eco.boot", contains = "RelDist")
setClass("eco.permut", contains = "RelDist")
setClass("eco.variog")
setClass("eco.gm")
setClassUnion("dataframeORmatrix",
							c("data.frame", "matrix"))

########--Plot---RelDist---#######################################################

# Plot for lag objects


setMethod("plot", "eco.boot", function(x, var = NULL) {
  
	if(length(x@OUT) == 1) {
		var2 <- 1
	} else if(is.null(var)) {
		stop("a variable (var) argument must be selected")
	} else {
		var2 <- which(names(x@OUT) %in% var)
	}
  
  var <- as.numeric(var2)
  datos <- as.data.frame(x@OUT[[var2]])
  
  ylabel <- unlist(strsplit(x@METHOD, ""))
  ylabel[1] <- toupper(ylabel[1])
  ylabel <- paste(ylabel, collapse = "")
  xlabel <- "Mean lag distance (m)"
  title <- names(x@OUT[var2])
  
  method <- x@METHOD
  
  method2 <- pmatch(method, c("moran", "geary", "mantel"))
  if(is.na(method2)) {
    stop("not available plot method yet for", method, "objects")
  }
  
  z <- ggplot2::ggplot(data = datos) + 
    ggplot2::geom_line(ggplot2::aes(x = d.mean, y = uppr), 
                       directions = "hv", 
                       linetype = 2, 
                       colour = "red") +  
    ggplot2::geom_line(ggplot2::aes(x = d.mean, y = lwr),
                       direction = "hv", 
                       linetype = 2, 
                       colour = "red") + 
    ggplot2::geom_point(ggplot2::aes(x = d.mean, y = est)) + 
    ggplot2::geom_line(ggplot2::aes(x = d.mean, y = est)) + 
    ggplot2::xlab(xlabel) + 
    ggplot2::ylab(ylabel) + 
    ggplot2::labs(title = title) + 
    ggplot2::geom_ribbon(ggplot2::aes(x = d.mean, ymax = lwr, 
    																	ymin = uppr), 
                         fill = "red",
                         alpha = 0.05) + 
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
                   axis.title = ggplot2::element_text(size = 14, 
                   																	 face = "bold"), 
                   plot.title = ggplot2::element_text(size = 16,
                   																	 face = "bold")) 

  
  return(z)
})

# Plot for eco.multiboot objects

setMethod("plot", "eco.multiboot", function(x, var = NULL) {
	
	if(length(x@OUT) == 1) {
		var <- x@NAMES
		var2 <- 1
	} else if(is.null(var)) {
		stop("a variable (var) argument must be selected")
	} else {
		var2 <- which(names(x@OUT) %in% var)
	}
	
	method <- x@METHOD
	ylabel <- "Mean lag distance (m)"
	datos <- x@OUT[[var2]][[1]]
	title <- names(x@OUT[var2])
	
	if(method == "joincount") { 
		xlabel <- "Join class"
		legend <- "Join-Count"
	} else {
		xlabel <- "Individual"
		legend <- "Getis-Ord's"
	}
	
	
	
	x <- rownames(datos)
	y <- colnames(datos)
	df<-expand.grid(x,y)
	statistic <- as.vector(t(datos))
	df$statistic <- statistic
	colnames(df) <- c("x", "y", "statistic")
	envir <- environment()
	
	z <- ggplot2::ggplot(df, ggplot2::aes(x = x,y = y), environment = envir) +
		ggplot2::geom_tile(ggplot2::aes(fill = statistic))+ 
		ggplot2::scale_fill_gradient2(legend, low = "blue", high = "red", space = "Lab")+
		ggplot2::xlab(xlabel)+ 
		ggplot2::ylab(ylabel)+
		ggplot2::labs(title = paste(var)) +
		ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
									 axis.title = ggplot2::element_text(size = 14, face = "bold"), 
									 plot.title = ggplot2::element_text(size = 16, face = "bold"))
	
	return(z)
  
})

# Plot for eco.permut objects

 
setMethod("plot", "eco.permut", function(x, var = NULL) {
	
	if(length(x@OUT) == 1) {
		var2 <- 1
	} else if(is.null(var)) {
		stop("a variable (var) argument must be selected")
	} else {
		var2 <- which(names(x@OUT) %in% var)
	}
	
	datos <- as.data.frame(x@OUT[[var2]])
	
	ylabel <- unlist(strsplit(x@METHOD, ""))
	ylabel[1] <- toupper(ylabel[1])
	ylabel <- paste(ylabel, collapse = "")
	xlabel <- "Distance (m)"
	title <- names(x@OUT[var2])
	
	method <- x@METHOD
	
	method2 <- pmatch(method, "mantel")
	if(is.na(method2)) {
		stop("not available plot method for this object")
	}
	
	pval2 <- as.numeric(datos$pval < 0.05)
	pval2[pval2 == 1] <- "P<0.05"
	pval2[pval2 == 0] <- "NS"
	
	localenv <- environment()
	
	z <- ggplot2::ggplot(datos, environment = localenv) + 
		ggplot2::geom_point(ggplot2::aes(x = d.mean, y = est, colour = pval2),  
												size = 5) + 
		ggplot2::geom_line(ggplot2::aes(x = d.mean, y = est)) + 
		ggplot2::xlab(xlabel) + 
		ggplot2::ylab(ylabel) + 
		ggplot2::labs(title = title) + 
		ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
									 axis.title = ggplot2::element_text(size = 14, face = "bold"), 
									 plot.title = ggplot2::element_text(size = 16, face = "bold")) + 
		ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
									 axis.title = ggplot2::element_text(size = 14, face = "bold"))+
		ggplot2::scale_colour_discrete(name  ="P value",
																	 labels=c("NS", "P < 0.05"))
	
	return(z)
})

# Plot for eco.variog objects

setMethod("plot", "eco.variog", function(x) {
  
  ylab <- "Semivariance"
  xlab <- "Distance"
  
  tab <- data.frame(sapply(x,c))
  
  if(ncol(tab) == 2) {
    na.dum <- rep(NaN, nrow(tab))
    tab <- cbind(tab, na.dum, na.dum)
    names(tab)[3:4] <- c("lwr", "uppr")
  }
  
  
  z <- ggplot2::ggplot(data = tab) + 
    ggplot2::geom_line(ggplot2::aes(x = d.mean, y = est),
                       directions="hv",linetype = 1, 
                       colour = "red") +
    ggplot2::geom_point(ggplot2::aes(x = d.mean,y = est),colour = "black") +
    ggplot2::ylab(ylab) + 
    ggplot2::xlab(xlab) + 
    ggplot2::theme(axis.title = ggplot2::element_text(size = 14,face = "bold"),
                   plot.title = ggplot2::element_text(size = 16,face = "bold")) +
    ggplot2::geom_line(ggplot2::aes(x = d.mean, y = uppr), 
                       directions = "hv", 
                       linetype = 2, 
                       colour = "red") +  
    ggplot2::geom_line(ggplot2::aes(x = d.mean, y = lwr),
                       direction = "hv", 
                       linetype = 2, 
                       colour = "red") + 
    ggplot2::geom_ribbon(ggplot2::aes(x = d.mean, ymax = lwr, ymin = uppr), 
                         fill = "red",
                         alpha = 0.05) 
  
  return(z)
  
})



####################---MlmModel-Class----#########################################

# eco.mlm-class


setClass("eco.mlm",
         
         representation( MLM = "list",
                         SUMMARY.MLM = "list", 
                         ANOVA.MLM = "list",
                         df1 = "data.frame",
                         df2 = "data.frame"),
)

# eco.mctree-class


setClass("eco.mctree",
         
         representation( TREES = "list",
                         PREDICTIONS = "list", 
                         FREQUENCIES = "list",
                         df1 = "data.frame",
                         df2 = "data.frame"), 
)

####Show lmtree #################################################################

# Show eco.mlm


setMethod("show", 	"eco.mlm", function(object) {
  
  cat("\n\n", paste(" @MLM:", "multiple model results"), "\n", 
      paste(" @SUMMARY.MLM:", "summary of the results"), "\n", 
      paste(" @ANOVA.MLM:", "analysis of variance tables"), "\n\n")
})


# Show eco.mctree

setMethod("show", 	"eco.mctree", function(object) {
  
  cat("\n\n", paste(" @TREES:", "trees"), "\n", 
      paste(" @PREDICTIONS:", "predictions of the analysis"), "\n", 
      paste(" @FREQUENCIES:", "Number of individuals predicted in each node", "\n\n"))
})


#################### Summary eco.lmtree #########################################

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

# Summary for eco.mctree output



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
		
