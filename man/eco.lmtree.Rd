% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eco.lmtree.R
\name{eco.lmtree}
\alias{eco.lmtree}
\title{Fitting Multiple Linear Regression models by stepwise AIC selection and
Multiple Classification and Regression Trees via party}
\usage{
eco.lmtree(
  df1,
  df2,
  analysis = c("mlm", "mctree"),
  mod.class = "+",
  fact = NULL,
  ...
)
}
\arguments{
\item{df1}{Data frame with dependent variables as columns.}

\item{df2}{Data frame with independent variables as columns.}

\item{analysis}{Class of analysis to perform. "mlm" 
for multiple linear regression analysis, 
or "mctree" for a multiple classification tree analysis.}

\item{mod.class}{"+" for additive model, "*" for model with 
interaction, in both cases, these models will include all terms
in the dependent data frame. If other model than these two 
is desired, it could be specified as a string 
with the names of those columns of the independent variable 
that should be used as terms. This string corresponds
to the right side "x" of a formula y ~ x (see examples).}

\item{fact}{Optional factor for estimating the frequencies
of individuals from different levels in each node, when the
analysis performed is "mctree".}

\item{...}{Further arguments passed to \code{\link[stats]{lm}} or
\code{\link[party]{ctree}}}
}
\value{
When the analysis selected is "mlm", the output object 
has three main slots:

> MLM: the results of the model

> SUMMARY.MLM the summary for each variable returned by the \code{\link{lm}} 
function 

> ANOVA.MLM with the ANOVAs results.


When the analysis selected is "mctree", the output object 
has also three main slots:

> TREES: Trees returned by the multiple \code{\link[party]{ctree}} analysis.

> PREDICTIONS: Predictions of the analysis.

> FREQUENCIES: Number of individuals predicted in each node.


\strong{ACCESS TO THE SLOTS}
The content of the slots can be accessed 
with the corresponding accessors, using
the generic notation of EcoGenetics 
(<ecoslot.> + <name of the slot> + <name of the object>).
See help("EcoGenetics accessors") and the Examples
section below
}
\description{
This program fits for each dependent variable, a Multiple Linear
Regression model calling the function \code{\link[stats]{step}} for choosing 
the best model by AIC criterion, or a Multiple Classification and Regression Trees
model, using the package party. 
The summary of the model returns information about 
the significance of the models, F-statistics and degrees of freedom,
when is fitted a "mlm"; otherwise, when the model fitted is a "mctree", the summary
returns the plots of those trees with significant splits.
}
\examples{

\dontrun{

data(eco.test)

# mlm additive model
mod <- eco.lmtree(df1 = eco3[["P"]], df2 = eco3[["E"]], 
analysis = "mlm")                                 
mod
summary(mod)

# mctree additive model
mod <- eco.lmtree(df1 = eco3[["P"]], df2 = eco3[["E"]], 
analysis = "mctree", fact = eco3[["S"]]$pop) 


#-----------------------
# ACCESSORS USE EXAMPLE
#-----------------------

# the slots are accessed with the generic format 
# (ecoslot. + name of the slot + name of the object). 
# See help("EcoGenetics accessors")

summary(mod)

ecoslot.FREQUENCIES(mod)        # slot FREQUENCIES

# frequency table with counts of individuals in populations x terminal nodes
tabfreq <- do.call(cbind, ecoslot.FREQUENCIES(mod))
namestab <- lapply(ecoslot.FREQUENCIES(mod), ncol)
namestab <- lapply(namestab, rep)
namestab <- rep(names(namestab), namestab)
colnames(tabfreq) <- namestab
tabfreq

# mlm custom model
mymod <- "E1+E2*E3"
mod <- eco.lmtree(df1 = eco[["P"]], df2 = eco[["E"]], 
analysis = "mlm", mod.class = mymod)            
summary(mod)

# mctree custom model
mod <- eco.lmtree(df1 = eco[["P"]], df2 = eco[["E"]], 
analysis = "mctree", mod.class = mymod, fact = eco[["S"]]$pop)   
 
                      
summary(mod)

}

}
\references{
Hothorn T., K. Hornik, and A. Zeileis. 2006. Unbiased Recursive Partitioning: 
A Conditional Inference Framework. Journal of Computational and Graphical Statistics, 
15: 651-674.
}
\author{
Leandro Roser \email{learoser@gmail.com}
}
