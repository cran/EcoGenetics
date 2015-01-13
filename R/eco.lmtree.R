#' Fitting Multiple Linear Regression models by stepwise AIC selection and
#' Multiple Classification and Regression Trees via party.
#' @param df1 Data frame with dependent variables as columns.
#' @param df2 Data frame with independent variables as columns.
#' @param analysis Class of analysis to perform. "mlm" 
#' for multiple linear regression analysis, 
#' or "mctree" for a multiple classification tree analysis.
#' @param mod.class "+" for additive model, "*" for model with 
#' interaction, in both cases, these models will include all terms
#' in the dependent data frame. If other model than these two 
#' is desired, it could be specified the model as a string 
#' with the names of those columns of the independent variable 
#' that should be used as terms. This string corresponds
#' to the right side x of a formula y ~ x (see examples).
#' @param fact Optional factor for estimating the frequencies
#' of individuals from different levels in each node, when the
#' analysis performed is "mctree".
#' @param ... Further arguments passed to \code{\link[stats]{lm}} or
#' \code{\link[party]{ctree}}
#' @description This program fits for each dependent variable, a Multiple Linear
#' Regression model calling the function \code{\link[stats]{step}} for choosing 
#' the best model by AIC criterion, or a Multiple Classification and Regression Trees
#' model, using the package party. 
#' The summary of the model returns information about 
#' the significance of the models, F statistics and degrees of freedom,
#' when is fitted a "mlm"; otherwise, when the model fitted is a "mctree", the summary
#' returns the plots of those trees with significant splits. 
#' @return 
#' When the analysis selected is "mlm", the output object 
#' has three main slots:
#' @@MLM: the results of the model, @@SUMMARY.MLM with the 
#' summary for each variable returned by the \code{\link{lm}} 
#' function and @@ANOVA.MLM with the ANOVAs results.
#' When the analysis selected is "mctree", the output object 
#' has also three main slots:
#' @@TREES: Trees returned by the \code{\link[party]{ctree}} analysis.
#' @@PREDICTIONS: Predictions of the analysis.
#' @@FREQUENCIES: Number of individuals predicted in each node.
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @examples
#' \dontrun{
#' 
#' data(eco.test)
#' mod <- eco.lmtree(df1 = eco$P, df2 = eco$E, 
#' analysis = "mlm")                                #mlm aditive model 
#' mod
#' summary(mod)
#' 
#' mod <- eco.lmtree(df1 = eco$P, df2 = eco$E, 
#' analysis = "mctree", fact = eco$S$structure)     # mctree aditive model 
#' mod@@FREQUENCIES
#' summary(mod)
#' 
#' mymod <- "E1+E2*E3"
#' mod <- eco.lmtree(df1 = eco$P, df2 = eco$E, 
#' analysis = "mlm", mod.class = mymod)             #mlm custom model
#' summary(mod)
#' 
#' mod <- eco.lmtree(df1 = eco$P, df2 = eco$E, 
#' analysis = "mctree", mod.class = mymod, 
#' fact = eco$S$structure)                          # mctree custom model
#' summary(mod)
#' 
#' }
#' @export

setGeneric("eco.lmtree", 
					 function(df1, df2, 
										analysis = c("mlm", "mctree"), mod.class = "+", 
										fact = NULL, ...) 	{

	           
  analysis <- match.arg(analysis)
  
  if((mod.class == "+") || (mod.class == "*")) {
  indep <- paste(colnames(df2), collapse = mod.class, sep = "")
  } else {
  	indep <- mod.class
  }
  
  data <- data.frame(df1, df2)
  
  if(analysis == "mlm") {
  	
  	mlm.mod <- new("eco.mlm")
  	mod <- results <- anovas <- list()
  	
  for(i in 1:ncol(df1)) {
    dep <- paste(colnames(df1)[i], "~", sep = "")
    smod <- paste(dep, indep, sep = "")
    smod <- as.formula(smod)
    mod[[i]] <- lm(smod, data = data, ...)
    mod[[i]] <- step(mod[[i]], scope = list(colnames(df1)[i] ~ 1,
                                            upper = mod[[i]]))
  }
  
  for(i in 1:ncol(df1)) {
    results[[i]] <- summary(mod[[i]])
    anovas[[i]] <- anova(mod[[i]])
  }
  
  names(mod) <- colnames(df1)
  names(results) <- colnames(df1)
  names(anovas) <- colnames(df1)
  
  mlm.mod@MLM <- mod
  mlm.mod@SUMMARY.MLM <- results
  mlm.mod@ANOVA.MLM <- anovas
  mlm.mod@df1 <- df1
  mlm.mod@df2 <- df2
  
  mlm.mod
  return(mlm.mod)
}

else if(analysis == "mctree") {
	
	tre.new <- new("eco.mctree")
	
	tre <- list()
	prediction <- list()
	freq <- list()
	
	for(i in 1:ncol(df1))
	{
		dep <- paste(colnames(df1)[i], "~", sep = "")
		smod <- paste(dep, indep, sep = "")
		smod <- as.formula(smod)
		tre[[i]] <- party::ctree(smod, data = data, ...)
		prediction[[i]] <- party::where(tre[[i]])
		if(!is.null(fact)) {
		freq[[i]] <- table(fact, prediction[[i]])
		}
	}
	names(tre) <- colnames(df1)
	names(prediction) <- colnames(df1)
	if(!is.null(fact)) {
	names(freq) <- colnames(df1)
	} 
	if(is.null(fact)){
		freq[[1]] <- "factor not provided for computing frequencies"
	}
	
	tre.new@TREES <- tre
	tre.new@PREDICTIONS <- prediction
	tre.new@FREQUENCIES <- freq
	return(tre.new)
}
})

