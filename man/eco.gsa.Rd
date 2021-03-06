% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eco.gsa.R
\name{eco.gsa}
\alias{eco.gsa}
\title{Global spatial analysis}
\usage{
eco.gsa(
  Z,
  con,
  Y = NULL,
  method = c("I", "C", "CC", "JC"),
  nsim = 99,
  alternative = c("auto", "two.sided", "greater", "less"),
  ploidy = 1,
  adjust = "fdr",
  plotit = TRUE
)
}
\arguments{
\item{Z}{Vector with a variable, or matrix/data frame with variables in columns.}

\item{con}{An object of class eco.weight obtained with the function \code{\link{eco.weight}},
a listw object or a matrix, giving the spatial weights for the analysis. If "con" is a matrix,
an attribute "xy" including the projected coordinates is required.}

\item{Y}{Vector with the second variable for Moran's Ixy.
If Z has multiple variables, the program will compute the coefficent 
for each with Y.}

\item{method}{Method of analysis: "I" for Moran's I, "C" for Geary's C, "CC" for
the Bivariate Moran's or "JC" for Join-count.}

\item{nsim}{Number of Monte-Carlo simulations.}

\item{alternative}{The alternative hypothesis. If "auto" is selected (default) the
program determines the alternative hypothesis.
Other options are: "two.sided", "greater" and "less".}

\item{ploidy}{For join count analysis: number of elements for the values of the vector passed, given value: for example,
if ploidy=1, "13" and "31" are considered a same level ("31" is sorted by the program as "13"); 
if ploidy = 1, "13" and "31" represent two different levels.}

\item{adjust}{Correction method of P-values for multiple tests, 
passed to \code{\link[stats]{p.adjust}}. Defalut is "fdr".}

\item{plotit}{should be generated a plot for univariate results?}
}
\value{
The program returns an object of class "eco.gsa" with the following slots:

> METHOD method used in the analysis

> OBS observed value when a single variable is tested

> EXP expected value when a single variable is tested

> PVAL P-value when a single variable is tested

> ALTER alternative hypotesis when a single variable is tested

> NSIM number of simulations

> MULTI table with observed and expected values, P-values and alternative
hypoteses when multiple variables are tested


\strong{ACCESS TO THE SLOTS}
The content of the slots can be accessed 
with the corresponding accessors, using
the generic notation of EcoGenetics 
(<ecoslot.> + <name of the slot> + <name of the object>).
See help("EcoGenetics accessors") and the Examples
section below.
}
\description{
Univariate and multivariable global spatial analysis.This program computes Moran's I, Geary's C,  Bivariate Moran's I or
Join-count statistics with P-values.

The program allows the analysis of a single variable or multiple variables. 
in this latter last case, the variables must be in columns and the individuals in rows.

In join-count analysis, a ploidy argument must be supplied. The data is then ordered with the 
function \code{\link{aue.sort}}. This step is required for the analysis of
genotypic data. An individual with the alleles A and B, coded as AB, is identical 
to other coded as BA. The order step ensures the AB and BA 
will be considered the same genotype. For the analysis of frequencies of single alleles, the input 
is count data (ploidy times the frequency, as provided by the slot
A of an ecogen object: the count data A' can be obtained as A' <- ploidy * A),
using the function with the arguments ploidy = 1.
}
\examples{

\dontrun{

data(eco.test)

# Moran's I 

### one test
con <- eco.weight(eco[["XY"]], method = "circle", d1 = 0, d2 = 2)
global <- eco.gsa(Z = eco[["P"]][, 1], con = con, method = "I", nsim = 200)
global

require(adegenet)
con2<-chooseCN(eco[["XY"]], type = 1, result.type = "listw", plot.nb = FALSE)
global <- eco.gsa(Z = eco[["P"]][, 1], con = con2, method = "I", nsim = 200)
global

#-----------------------
# ACCESSORS USE EXAMPLE
#-----------------------

# the slots are accesed with the generic format 
# (ecoslot. + name of the slot + name of the object). 
# See help("EcoGenetics accessors")

# observed value
ecoslot.OBS(global)

# p-value
ecoslot.PVAL(global)

#----------------
# multiple tests
#----------------
data(eco3)
con <- eco.weight(eco3[["XY"]], method = "circle", d1 = 0, d2 = 500)
global <- eco.gsa(Z = eco3[["P"]], con = con, method = "I", nsim = 200)
global 
# Plot method for multivariable eco.gsa objects:
eco.plotGlobal(global)

#--------------------------------
# accessor use in multiple tests
#--------------------------------

ecoslot.MULTI(global)


#----------------------------------------

# Gearys's C 

con <- eco.weight(eco[["XY"]], method = "circle", d1 = 0, d2 = 2)
global.C <- eco.gsa(Z = eco[["P"]][, 1], con = con, method = "C", nsim = 200)
global.C

#----------------------------------------

# Bivariate's Moran's Ixy

con <- eco.weight(eco[["XY"]], method = "circle", d1 = 0, d2 = 2)
global.Ixy <- eco.gsa(Z = eco[["P"]][, 1], Y = eco[["E"]][, 1], 
con = con, method = "CC", nsim = 200)
global.Ixy

#----------------------------------------

# Join-count

## using the allelic frequency matrix of an ecogen object. 
## The data is diploid. Frequencies are transformed into counts
## as ploidy * frequency_matrix:

Z = 2* eco[["A"]]


jc <- eco.gsa(Z[, 1], con = con, method = "JC")
eco.plotGlobal(jc)

# multiple tests
# using the first ten alleles of the matrix
global.JC <- eco.gsa(Z[, 1:10], con = con, method = "JC", nsim = 99)
global.JC
# plot method for multivariable join-count
eco.plotGlobal(global.JC)

# counting joins between genotypes in the first locus the G matrix:
global.JC <- eco.gsa(Z = eco[["G"]][, 1], ploidy = 2, con = con, method = "JC", nsim = 99)
global.JC
eco.plotGlobal(global.JC)

}

}
\references{
Geary R. 1954. The contiguity ratio and statistical mapping. 
The incorporated statistician, 115-146.

Moran P. 1950. Notes on continuous stochastic phenomena. Biometrika, 17-23. 

Reich R., R. Czaplewski and W. Bechtold. 1994. 
Spatial cross-correlation of undisturbed, natural shortleaf pine stands 
in northern Georgia. Environmental and Ecological Statistics, 1: 201-217.

Sokal R. and N. Oden 1978. Spatial autocorrelation in biology: 
1. Methodology. Biological journal of the Linnean Society, 10: 199-228.

Sokal R. and N. Oden. 1978. Spatial autocorrelation in biology. 
2. Some biological implications and four applications of evolutionary and 
ecological interest. Biological Journal of the Linnean Society, 10: 229-49.
}
\author{
Leandro Roser \email{learoser@gmail.com}
}
