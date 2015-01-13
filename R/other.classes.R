########--RelDist-Class--#####################################################
#' RelDist-class
#' @name RelDist-class
#' @keywords internal
#' @slot OUT output from eco.autocor analysis
#' @slot INTERVAL interval
#' @slot MAX max distance
#' @slot TYPE type
#' @slot NAMES names
#' @slot METHOD method 
#' @include ecogen.definition.R
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @aliases RelDist-class

setClass("RelDist", 
				 
				 representation(OUT = "list",
				 							 INTERVAL = "numeric", 
				 							 MAX = "numeric",
				 							 TYPE = "character", 
				 							 NAMES = "character",
				 							 METHOD = "character"), 
				 
				 prototype( OUT = list(), 
				 					 INTERVAL = 0,
				 					 MAX = 0, 
				 					 TYPE = character(), 
				 					 NAMES = character(), 
				 					 METHOD = character())
)


########--Plot---RelDist---###################################################


#' Plot for RelDist objects.
#' @param x RelDist object.
#' @param var Variable to plot (numeric, see examples).
#' @seealso  \code{\link{eco.autocor}}
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @examples
#' \dontrun{
#' 
#' eco.ac <- eco.autocor(eco, int = 50, smax = 1000,
#' d.class = "phenotype")
#' plot(eco.ac, var = 6)
#' plot(eco.ac, var = 3)
#' 
#' }
#' @rdname RelDist-methods
#' @aliases plot,RelDist-method
#' @exportMethod plot 


setMethod("plot", "RelDist", function(x, var) {
	
	var <- as.numeric(var)
	datos <- as.data.frame(x@OUT[[var]])
	
	ylabel <- unlist(strsplit(x@METHOD, ""))
	ylabel[1] <- toupper(ylabel[1])
	ylabel <- paste(ylabel, collapse = "")
	xlabel <- "Distance"
	title <- names(x@OUT[var])
	
	z <- ggplot2::ggplot(data = datos) + 
		ggplot2::geom_line(ggplot2::aes(x = d.mean, y = upr), 
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
		ggplot2::geom_ribbon(ggplot2::aes(x = d.mean, ymax = lwr, ymin = upr), 
								fill = "red",
								alpha = 0.05) + 
		ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
					axis.title = ggplot2::element_text(size = 14, face = "bold"), 
					plot.title = ggplot2::element_text(size = 16, face = "bold")) + 
		ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
					axis.title = ggplot2::element_text(size = 14, face = "bold"))
	
	#geom_errorbar(aes(x = d.mean, ymax = est + sd, ymin= est-sd))
	
	return(z)
})


####################---MlmModel-Class----#####################################
#' eco.mlm-class
#' @name eco.mlm-class
#' @keywords internal
#' @slot MLM mlm results
#' @slot SUMMARY.MLM summary of the results
#' @slot ANOVA.MLM anovas for the results
#' @slot df1 data frame
#' @slot df2 data frame
#' @include ecogen.definition.R
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @aliases eco.mlm-class

setClass("eco.mlm",
				 
				 representation( MLM = "list",
				 								SUMMARY.MLM = "list", 
				 								ANOVA.MLM = "list",
				 								df1 = "data.frame",
				 								df2 = "data.frame"),
				 )

#' eco.mctree-class
#' @name eco.mctree-class
#' @keywords internal
#' @slot TREES trees obtained
#' @slot PREDICTIONS predictions of the analysis
#' @slot FREQUENCIES frequencies of individuals per class in nodes
#' @slot df1 data frame
#' @slot df2 data frame
#' @include ecogen.definition.R
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @aliases eco.mctree-class


setClass("eco.mctree",
				 
				 representation( TREES = "list",
				 								 PREDICTIONS = "list", 
				 								 FREQUENCIES = "list",
				 								df1 = "data.frame",
				 								df2 = "data.frame"), 
				 )

####Show lmtree #####################################################
#' show eco.mlm
#' @keywords internal
#' @rdname eco.lmtree-methods
#' @aliases show,eco.mlm-method


setMethod("show", 	"eco.mlm", function(object) {
	
	cat("\n\n", paste(" @MLM:", "multiple model results"), "\n", 
			paste(" @SUMMARY.MLM:", "summary of the results"), "\n", 
			paste(" @ANOVA.MLM:", "analysis of variance tables"), "\n\n")
})

#' show eco.mctree
#' @keywords internal
#' @rdname eco.lmtree-methods
#' @aliases show,eco.mctree-method


setMethod("show", 	"eco.mctree", function(object) {
	
	cat("\n\n", paste(" @TREES:", "trees"), "\n", 
			paste(" @PREDICTIONS:", "predictions of the analysis"), "\n", 
			paste(" @FREQUENCIES:", "Number of individuals predicted in each node", "\n\n"))
})


#################### Summary eco.lmtree #########################################


#' Summary for eco.lmtree output
#' @param object output object of \code{\link{eco.lmtree}}
#' @return A table with a summary of the analysis for "mlm" analysis, 
#' the plot of the trees with significant splits for "mctree" analysis.
#' @seealso \code{\link{eco.lmtree}}
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @examples
#' \dontrun{
#' 
#' data(eco.test)
#' #' mod <- eco.lmtree(df1 = eco$P, df2 = eco$E, 
#' analysis = "mlm")                                    
#' summary(mod)                                    #summary for "mlm" analysis
#' 
#' mod <- eco.lmtree(df1 = eco$P, df2 = eco$E,
#' analysis = "mctree", fact = eco$S$structure)               
#' summary(mod)                                    #summary for "mctree" analysis
#' 
#' }
#' @rdname eco.lmtree-summary
#' @aliases summary,eco.lmtree-method
#' @exportMethod summary

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

#' @rdname eco.lmtree-summary
#' @aliases summary,eco.lmtree-method
#' @exportMethod summary

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
		
