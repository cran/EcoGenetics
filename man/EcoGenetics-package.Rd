\name{EcoGenetics-package}
\alias{EcoGenetics-package}
\alias{EcoGenetics}
\docType{package}
\title{
 Spatial Analysis of Phenotypic, Genotypic and Environmental Data
}
\description{
EcoGenetics provides geostatistical tools for the spatial analysis of phenotypic, genotypic and environmental data. The package has been designed to make easy the storing, handling and integration of the available information from multiple sources, under an S4 philosophy. 
}


\details{
 \tabular{ll}{
 Package: \tab EcoGenetics\cr
 Type: \tab Package\cr
 Version: \tab 1.2.0-3\cr
 Date: \tab 2016-07-15\cr
 License: \tab GPL (>=2) \cr
 }
 
 
 
 
 ----------------------- o-----------------------
 
  \bold{\emph{I. STRUCTURE OF THE PACKAGE}}

EcoGenetics has four basic modules. The \strong{base module} is composed by general functions (\link[=eco.lmtree]{multiple lm}, \link[=eco.detrend]{detrending spatial data utility}, etc.). 
The \strong{general spatial module} computes global (\link[=eco.gsa]{Moran's I}, \link[=eco.mantel]{Mantel test}, etc.) and local (\link[=eco.lsa]{Getis-Ord's G*}, \link[=eco.lsa]{local Moran's I}, etc.) spatial tests. These analyses use a \link[=eco.weight]{spatial weights matrix}, provided by the  \strong{spatial weights module}.  The \strong{lag analysis module} performs two basic analyses: the obtention  of \link[=eco.variogram]{variograms} and \link[=eco.correlog]{correlograms} (see also \link[=eco.cormantel]{this link}). This module uses the other tool provided by the spatial weights module: a sequence of \link[=eco.lagweight]{spatial weights matrices}. 

The package have also special plot methods, as \link[=eco.rankplot]{rankplot},  
\link[=eco.forestplot]{forestplot}, \link[=eco.rasterplot]{rasterplot} (implemented the last two for 
\link[=eco.lsa]{local spatial analysis}). Several conversor of data to other programs are available, (as to \link[=eco.2genepop]{genepop} - an \link[=eco.genepop2df]{importer} tool is also defined for genepop-, \link[=eco.2spagedi]{SPAGeDi}, etc.). Basic manipulation of genetic matrices is allowed by \link{eco.convert} and \link{eco.format}. Tools for computation of NDVI in Landsat imaginery, 
post-process of rasters and temporal analysis can be found in \link{eco.NDVI}, \link{eco.NDVI.post} and \link{eco.theilsen}. Other useful functions are \link{aue.sort} (for ordering alleles), \link{eco.alfreq} (to comput histograms of allelic frequencies for detection of bottlenecks), \link{eco.post.geneland} and \link{eco.pairtest}. 

The results obtained with the main functions defined in EcoGenetics are
object of class \link[=Methods]{S4}. As a default characteristic of the package design,
these objects have a "show"" method for a general overview of the results, and methods to extract the results stored in slots (generic \link[=ecoslot.XY]{accessors} and double square brackets (\link[=Extract]{"[["}) definitions). 

For storing and pre-processing the data to analyze, the package 
defines a special class: the class \code{\link{ecogen}}.
 
 \bold{\emph{II. STRUCTURE OF ECOGEN OBJECTS: HANDLING AND INTEGRATING INFORMATION}}


Ecological genetics research requires the integration of data originated in different sources. The class \code{\link{ecogen}} has been designed for handling multidimensional data. Its basic structure is the following:

- An \strong{XY} slot, storing a data frame with geographic coordinates.  

- A \strong{P} slot, storing a phenotypic data frame.  

- A \strong{G} slot, storing a genotypic data frame. 

- An \strong{A} slot containing as allelic frequencies the information of G (only
available for codominant markers.)

- An \strong{E} slot, storing an environmental data frame.  

- A \strong{S} slot, storing a data frame with classes assigned to the 
individuals. 


- A \strong{C} slot, for a custom data frame. 

- An \strong{OUT} slot, containing a list for the storage of the results.


For dominant (presence-absence) markers, the slot A is empty.


\bold{\emph{III. A BRIEF OVERVIEW OF THE DATA HANDLING METHODS DEFINED FOR ECOGEN OBJECTS}}


The construction of a new "ecogen" object from a data frame is made with the homonymous function.

\emph{library("EcoGenetics")}

\emph{data(eco.test)}

\emph{eco <- ecogen(XY = coordinates, P = phenotype, G = genotype, E = environment, S = structure, order.G = TRUE)}


# The following methods can be used with ecogen objects:

# - - - - - - 

\bold{# 1. subsetting/ordering by row-method, using single square brackets ("[")}

\emph{eco.sub <- eco[1:50]}

\emph{eco.ord <- eco[c(50:1)]}

\emph{eco.ord[["XY"]] # object ordered and subsetted rows.
This is a powerful method for handling the stack of ecogen data frames}

\emph{eco <- eco[order(eco[["S"]][, 1])] #ordering S[, 1] and then the object}

\bold{# 2. \link[=eco.merge]{merging-method}, for two objects}

\emph{eco1 <- eco}

\emph{merged <- eco.merge(eco, eco1)}

\bold{# 3. \link[=eco.subset]{subsetting-method}, in reference to a group in the S slot (in this case, "1")}

\emph{eco.subS <- eco.subset(eco,"pop", 1)}


\bold{# 4. \link[=eco.rbind]{binding by row-method} (duplicated row names not allowed)}

\emph{eco2 <- eco; rownames(eco2[["P"]]) <-226:450}

\emph{eco.r <- eco.rbind(eco, eco2)}

\bold{# 5. \link[=eco.cbind]{binding by column-method}}

\emph{eco.c <- eco.cbind(eco, eco1)}

\bold{# 6. \link[=eco.order]{ordering by row-method}, using the rows in the XY data frame as reference}

\emph{ecoslot.P(eco1) <- eco[["P"]][sample(1:173), ]} ## object with unordered rows in P data frame

\emph{ordered <- eco.order(eco1)}

\bold{# 7. get-method using generic \link[=ecoslot.XY]{accessors} and "[[" (equivalent methods):}

\emph{ecoslot.XY(eco) ; eco[["XY"]]}

\emph{ecoslot.P(eco); eco[["P"]]}

\emph{ecoslot.G(eco); eco[["G"]]}

\emph{ecoslot.A(eco); eco[["A"]]  (only available for codominant data}

\emph{ecoslot.E(eco); eco[["E"]]}

\emph{ecoslot.S(eco); eco[["S"]]}

\emph{ecoslot.C(eco); eco[["C"]]}

\emph{ecoslot.OUT(eco); eco[["OUT"]]}


\bold{# 8. set-method using \link[=ecoslot.XY]{accessors} and "[[" (equivalent methods):}

\emph{eco.temp <- ecogen(XY = coordinates, P = phenotype)}

\emph{eco.temp}

\emph{ecoslot.G(eco.temp, order.G = TRUE) <- genotype}

## this is equivalent, in square brackets notation, to:

\emph{eco[["G", order.G=TRUE]] <- genotype}

\emph{ecoslot.E(eco.temp) <- environment}

## identical to eco[["E"]] <- environment

\emph{ecoslot.S(eco.temp) <- structure}

## identical to eco[["S"]] <- structure

\emph{eco.temp}

\bold{# 9. appending-method (storing information generated via \link[=ecoslot.XY]{accessors)}}

## fitting a multiple linear regression model:

\emph{linear.analysis <- eco.lmtree(eco[["P"]], eco[["E"]], "mlm")}

## storing the results:

\emph{ecoslot.OUT(eco) <- linear.analysis}

\emph{eco}

## Storing multiple result at once:
 
\emph{a <- c(1:10)}

\emph{b <- c(2:30)}
  
\emph{ecoslot.OUT(eco) <- list(eco, a, b)}

\emph{eco}

## the use of the accessor OUT has its equivalence in double
square brackets notation:

\emph{eco[["OUT"]] <- list(eco, a, b)}


## summary table

\emph{ecoslot.OUT(eco)}
  
## the data frame shows the results stored in alphabetical order and their classes. The specification of a second name, return the corresponding stored object:

\emph{ecoslot.OUT(eco, "a")} 

## note that the append method is a particular case of the
## set method with accessors / "[[" for the slot OUT.

\bold{# 10. \link[=eco.remove]{removing-method}}

## removing objects a and b from eco

\emph{eco <- eco.remove(eco, a, b)}




- - - - - - 

Finally, the workspace can be \link[=eco.clear]{cleared}, only storing the desired object:

\emph{ls()}

\emph{eco.clear(eco)}

\emph{ls()}

 ----------------------- o -----------------------

}


\author{

 Leandro Roser, Juan Vilardi, Beatriz Saidman and Laura Ferreyra
 
 Maintainer: Leandro Roser <leandroroser@ege.fcen.uba.ar>
}

\references{
 
Anselin L. 1995. Local indicators of spatial association-LISA. Geographical analysis. 27: 93-115.

Borcard D., F. Gillet, and P. Legendre. 2011. Numerical ecology with R. Springer Science & Business Media.

Chander G., B. Markham, and D. Helder. 2009. Summary of current radiometric calibration coefficients for Landsat MSS, TM, ETM+, and EO-1 ALI sensors. Remote sensing of environment, 113: 893-903.

Chavez P. 1989. Radiometric calibration of Landsat Thematic Mapper multispectral images. Photogrammetric Engineering and Remote Sensing, 55: 1285-1294.

Chavez P. 1996. Image-based atmospheric corrections-revisited and improved. Photogrammetric engineering and remote sensing, 62: 1025-1035.

Double M., R. Peakall, N. Beck, and Y. Cockburn. 2005. Dispersal, philopatry, and infidelity: dissecting local genetic structure in superb fairy-wrens (Malurs cyaneus). Evolution 59: 625-635.

Dray S., and A., Dufour. 2007. The ade4 package: implementing the duality diagram for ecologists. Journal of statistical software, 22: 1-20.
 
Freedman D., and P. Diaconis. 1981. On the histogram as a density estimator: L 2 theory. Probability theory and related fields, 57: 453-476.
  
Geary R. 1954. The contiguity ratio and statistical mapping. The incorporated statistician, 115-146.

Getis A., and J. Ord. 1992. The analysis of spatial association by use of distance statistics. Geographical analysis, 24: 189-206.

Goslee S. 2011. Analyzing remote sensing data in R: the landsat package. Journal of Statistical Software, 43: 1-25.

Goudet J. 2005. Hierfstat, a package for R to compute and test hierarchical F-statistics. Molecular Ecology Notes, 5: 184-186.
 
Guillot G., F. Mortier and A. Estoup. 2005. GENELAND: a computer package for landscape genetics. Molecular Ecology Notes, 5: 712-715.
 
Jombart T. 2008. adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24: 1403-1405.

Kalisz S., J. Nason, F.M. Handazawa, and S. Tonsor. 2001. Spatial population genetic structure in Trillium grandiflorum:  the roles of dispersal, mating, history, and selection.  Evolution 55: 1560-1568.

Legendre P., and L. Legendre. 2012. Numerical ecology. Third English edition .Elsevier Science, Amsterdam, Netherlands.

Lichstein J., T. Simons, S. Shriner, and K. Franzreb. 2002. Spatial autocorrelation and autoregressive models in ecology. Ecological monographs, 72: 445-463.

Loiselle B., V. Sork, J. Nason, and C. Graham. 1995. Spatial genetic structure of a tropical understory shrub, Psychotria officinalis (Rubiaceae). American Journal of Botany 1420-1425.

Moran P. 1950. Notes on continuous stochastic phenomena. Biometrika, 17-23.

Oden N., and R. Sokal. 1986. Directional autocorrelation: an extension of spatial correlograms to two dimensions. Systematic Zoology, 35: 608-617.

Ord J., and A. Getis. 1995. Local spatial autocorrelation statistics: distributional issues and an application. Geographical analysis, 27: 286-306.

Reich R., R. Czaplewski and W. Bechtold. 1994. Spatial cross-correlation of undisturbed, natural shortleaf pine stands in northern Georgia. Environmental and Ecological Statistics, 1: 201-217.

Sokal R., and N. Oden 1978. Spatial autocorrelation in biology: 1. Methodology. Biological journal of the Linnean Society, 10: 199-228.

Sokal R., and N. Oden. 1978. Spatial autocorrelation in biology. 2. Some biological implications and four applications of evolutionary and ecological interest. Biological Journal of the Linnean Society, 10: 229-49.

Sokal R. 1979. Ecological parameters inferred from spatial correlograms. In: G. Patil and M. Rosenzweig, editors. Contemporary Quantitative Ecology and elated Ecometrics. International Co-operative Publishing House: Fairland, MD, pp. 167-96.

Sokal R. 1986. Spatial data analysis and historical processes. In: E. Diday, Y. Escoufier, L. Lebart, J. Pages, Y. Schektman, and R. Tomassone, editors. Data analysis and informatics, IV. North-Holland, Amsterdam, The Netherlands, pp. 29-43.

Sokal R., N. Oden and B. Thomson. 1998. Local spatial autocorrelation in a biological model. Geographical Analysis, 30: 331-354.

Sokal R., and B. Thomson. 2006. Population structure inferred by local spatial autocorrelation: an example from an Amerindian tribal population. American journal of physical anthropology, 129: 121-131.

Song C., C. Woodcock, K. Seto, M. Lenney and S. Macomber. 2001. Classification and change detection using Landsat TM data: when and how to correct atmospheric effects?. Remote sensing of Environment, 75: 230-244.

Sturges H. 1926. The choice of a class interval. Journal of the American Statistical Association, 21: 65-66.

Tucker C. 1979. Red and photographic infrared linear combinations for monitoring vegetation. Remote sensing of Environment, 8: 127-150.

Vekemans, X., and O. Hardy. 2004. New insights from fine-scale spatial genetic structure analyses in plant populations. Molecular Ecology, 
13: 921-935.

Wu. 1986. Jackknife, bootstrap and other resampling methods in regression analysis. the Annals of Statistics, 1261-1295.
 
}
 
\examples{

 \dontrun{
 
 
#---Detrending spatial data with polynomial interpolation---#
 
 
data(eco2)

## original data
data1 <- matrix(eco2[["P"]][,1], 30, 30)
image(data1)

## original data + trend
data2 <- matrix(eco2[["P"]][,2], 30, 30)
image(data2)

## data detrending
data2.det <- eco.detrend(Z = eco2[["P"]][,2], XY =  eco2[["XY"]], degree =  1)
data2.det <- ecoslot.RES(data2.det)
data2.det <- matrix(data2.det$df1, 30, 30)
image(data2.det)


#---Multiple Linear Regression fit---#
 
 
data(eco.test)
mymod <- "E1+E2*E3"
mod <- eco.lmtree(df1 = eco[["XY"]], df2 = eco[["E"]], 
analysis = "mlm", mod.class = mymod) 
summary(mod)
 
 
#---Multiple Conditional Inference Trees---#


data(eco.test)
mymod <- "E1+E2*E3"
mod <- eco.lmtree(df1 = eco[["P"]], df2 = eco[["E"]], 
analysis = "mctree", mod.class = mymod, fact = eco[["S"]]$structure)
summary(mod)


#---Global spatial analysis---#
 
 
## Moran's I 

### one test
data(eco.test)
con <- eco.weight(eco[["XY"]], method = "circle", d1 = 0, d2 = 2)
global <- eco.gsa(Z = eco[["P"]][, 1], con = con, , method = "I", nsim = 200)
global

require(adegenet)
con2<-chooseCN(eco[["XY"]], type = 1, result.type = "listw", plot.nb = FALSE)
global <- eco.gsa(Z = eco[["P"]][, 1], con = con2, , method = "I", nsim = 200)
global

### multiple tests
con <- eco.weight(eco[["XY"]], method = "circle", d1 = 0, d2 = 2)
global <- eco.gsa(Z = eco[["P"]], con = con, , method = "I", nsim = 200)
global


## Geary's C 

data(eco.test)
global.C <- eco.gsa(Z = eco[["P"]][, 1], con = con, method = "C", nsim = 200)
global.C


## Bivariate's Moran's Ixy

data(eco.test)
global.Ixy <- eco.gsa(Z = eco[["P"]][, 1], Y = eco[["E"]][, 1],
con = con, method = "CC", nsim = 200)
global.Ixy


## Join-Count

data(eco.test)
global.JC <- eco.gsa(Z = 2* eco[["A"]][, 1], ncod = 1,
con = con, method = "JC", nsim = 5)
global.JC


# Mantel test

data(eco.test)
eco.mantel(d1 = dist(eco[["P"]]), d2 = dist(eco[["E"]]), nsim = 99)  


## Partial Mantel test

data(eco.test)
eco.mantel(d1 = dist(eco[["P"]]), d2 = dist(eco[["E"]]),
dc = dist(eco[["XY"]]), nsim = 99)  


#---Local spatial analysis---#


## Getis-Ord's G*

data(eco.test)
require(ggplot2)
con<- eco.weight(eco[["XY"]], method = "knearest",  k = 4, self = TRUE) 
### self = TRUE for G*
getis.ak <- eco.lsa(eco[["P"]][, 1], con, method = "G*", nsim = 99, 
adjust = "none")
getis.ak

### to plot the results, the function "eco.lsa" calls "eco.rankplot"
### (see ?eco.rankplot) when test = "permutation" and "eco.forestplot"
### (see ?eco.forestplot) when test = "bootstrap"

p <- plot(getis.ak)   ###  rankplot graph
p                     ###  points with colors of the color-scale:
                      ###  points with P < 0.05. Yellow points : 
                      ###  points with P > 0.05
p <- plot(getis.ak, significant = FALSE)
p                     ###  all points have a color of the color-scale

### bootstrap example
getis.akb <- eco.lsa(eco[["P"]][, 1], con, method = "G*", nsim = 99,
test = "bootstrap")
p <- plot(getis.akb)     ### forestplot graph
p + ggplot2::theme_bw()  ### the plot can be modified with ggplot2
                         ### In this case, the background color is modified 


## Getis-Ord's G

data(eco.test)
require(ggplot2)
con<- eco.weight(eco[["XY"]], method = "knearest", k = 4)
### self = FALSE for G
getis <- eco.lsa(eco[["P"]][, 1], con, method = "G", nsim = 99, adjust = "none")
plot(getis)


## Local Moran's I

 #-------------------------
 # TESTING PHENOTYPIC DATA-
 #-------------------------
 
 con <- eco.weight(eco[["XY"]], method = "knearest",  k = 4, row.sd = TRUE) 
 # row standardized weights = TRUE
 
 # test for the first trait of the data frame P 
 localmoran <- eco.lsa(eco[["P"]][, 1], con, method = "I", nsim = 99)     
 
 plot(localmoran)
 
 # test for several variables
 
 all.traits <- apply(eco[["P"]], 2, eco.lsa,  con, method = "I", nsim = 99)
 
 # Observed statistic and P-values tables (individuals x traits)
 stat.P <- sapply(all.traits, function(x) return(ecoslot.OUT(x)[,1]))
 pval.P <- sapply(all.traits, function(x) return(ecoslot.OUT(x)[,4]))
 
 # Plot of the phenotypic spatial patterns
 
 par(mfrow = c(2,4))
 for(i in 1:8) {
 image(matrix(stat.P[,i], 15,15))
 }
 
 par(mfrow = c(2,4))
 for(i in 1:8) {
 image(matrix(pval.P[,i], 15,15))
 }
 
 
 #-------------------------
 # TESTING GENOTYPIC DATA-
 #-------------------------
 
 # eco[["A"]] is a matrix with the genetic data of "eco"
 # as frequencies for each allele in each individual
 # (use [["G"]] for presence-absence data).
 
 head(eco[["A"]])      # head of the matrix - 40 alleles
 
 con <- eco.weight(eco[["XY"]], method = "knearest",  k = 4, row.sd = TRUE) 
 # row standardized weights = TRUE
 
 # test for a single allele
 localmoran.geno <-  eco.lsa(eco[["A"]][, 32], con, method = "I", nsim = 99)
 
 # test for several alleles -  40 alleles (it runs in less than 1 min 
 # for 99 simulations per allele;  999 simulations takes ~ 11 s per allele, 
 # less than 8 min in total.) 
 all.alleles <- apply(eco[["A"]], 2, eco.lsa,  con, method = "I", nsim = 99)
 
 # plot all alleles to get an overview of the spatial patterns
 lapply(all.alleles, plot)
 
 # Observed statistic and P-values tables (individuals x loci)
 stat.G <- sapply(all.alleles, function(x) return(ecoslot.OUT(x)[,1]))
 pval.G <- sapply(all.alleles, function(x) return(ecoslot.OUT(x)[,4]))
 
 # counting individuals with P < 0.05 for each allele 
 # (5 * 225 /100 ~  12   significant tests by random)
 signif <- lapply(all.alleles, function(x) sum(ecoslot.OUT(x)[,4] < 0.05))
 signif <- unlist(signif)
 
 # filtering alleles, loci with > 12 significant individual tests
 
 A.local <- eco[["A"]][, signif > 12]     #filtered matrix
 stat.G.f <- stat.G[, signif > 12] 
 pval.G.f <- pval.G[, signif > 12]
 
 # Plot of the genotypic spatial patterns
 
 # one plot possibility, using the EcoGenetics method <rankplot>
 all.local <- all.alleles[signif > 12] 
 lapply(all.local, plot)
 
 # other plot possibility, using <image>
 par(mfrow = c(3,4))
 for(i in 1:12) {
 image(matrix(stat.G[,i], 15,15))
 }
 
 par(mfrow = c(3,4))
 for(i in 1:12) {
 image(matrix(pval.G[,i], 15,15))
 }


## Local Geary's C

data(eco.test)
require(ggplot2)
con<- eco.weight(eco[["XY"]], method = "knearest",  k = 4, row.sd = TRUE)
### row standardized weights = TRUE
localgeary <- eco.lsa(eco[["P"]][, 1], con, method = "C", nsim = 99, 
adjust = "none")
plot(localgeary)


#---Moran's I, Geary's C and bivariate Moran's Ixy correlograms---#


## Moran's I correlogram

### single test
data(eco.test)
require(ggplot2)
moran <- eco.correlog(Z=eco[["P"]][,1], XY = eco[["XY"]], method = "I",
smax=10, size=1000)
plot(moran)

### multiple tests
moran2 <- eco.correlog(Z=eco[["P"]], XY = eco[["XY"]], method = "I",
smax=10, size=1000)
plot(moran2, var ="P2") ## single plots
plot(moran2, var ="P3") ## single plots

graf <- plot(moran2, meanplot = TRUE)  ## multiple plot with mean correlogram
                                       ## and jackknifed confidence intervals.

plot(graf[[1]])
plot(graf[[2]])

### correlogram plots support the use of ggplot2 syntax
moranplot <- plot(moran2, var ="P3") + theme_bw() + theme(legend.position="none")
moranplot


## Geary's C correlogram

data(eco.test)
require(ggplot2)
geary <- eco.correlog(Z = eco[["P"]][,1], XY = eco[["XY"]], method = "C",
smax=10, size=1000)
plot(geary)


## Moran's Ixy cross-correlogram

data(eco.test)
require(ggplot2)
cross <- eco.correlog(Z = eco[["P"]][,1], XY = eco[["XY"]], Y = eco[["E"]][, 1],
method = "CC", int = 2, smax = 15)
plot(cross)


#---Mantel and partial Mantel correlograms---#


data(eco.test)
require(ggplot2)

## Mantel correlogram

corm <- eco.cormantel(M = dist(eco[["P"]]), size=1000,smax=7, XY = eco[["XY"]],
nsim = 99)
plot(corm)

corm <- eco.cormantel(M = dist(eco[["P"]]), size=1000,smax=7, XY = eco[["XY"]],
nsim = 99, test = "bootstrap")
plot(corm)



## Partial Mantel correlogram

corm <- eco.cormantel(M = dist(eco[["P"]]), MC = dist(eco[["E"]]),
size=1000, smax=7, XY = eco[["XY"]], nsim = 99)
plot(corm)

### correlogram plots support the use of ggplot2 syntax
mantelplot <- plot(corm) + theme_bw() + theme(legend.position="none")
mantelplot


#---Empirical variogram---#


data(eco.test)
require(ggplot2)
variog <- eco.variogram(Z = eco[["P"]][, 2],XY =  eco[["XY"]])
plot(variog)

### variogram plots support the use of ggplot2 syntax
variogplot <- plot(variog) + theme_bw() + theme(legend.position="none")
variogplot

 
#---Computing NDVI with atmospheric correction
#---over a time series, and estimation of the temporal mean---#

require(raster)
data(tab)
data(eco3)
set.seed(6)

temp <- list()

## we create 4 simulated rasters for the data included in the object tab:

for(i in 1:4) {
temp[[i]] <- runif(19800, 0, 254)
temp[[i]] <- matrix(temp[[i]], 180, 110)
temp[[i]] <- raster(temp[[i]], crs="+proj=utm")
extent(temp[[i]])<-c(3770000, 3950000, 6810000, 6920000)
}

writeRaster(temp[[1]], "20040719b4.tif", overwrite = T)
writeRaster(temp[[2]], "20040719b3.tif", overwrite = T)
writeRaster(temp[[3]], "20091106b4.tif", overwrite = T)
writeRaster(temp[[4]], "20091106b3.tif", overwrite = T)

## Computing NDVI images:

eco.NDVI(tab, "COST", "NDVI", "LT5")

## Mean NDVI image computed over the NDVI images that we calculated:

eco.NDVI.post(tab, "COST", "NDVI", what = c("mean", "var"))
mean.ndvi <- raster("NDVI.COST.mean.tif")
plot(mean.ndvi)

## Extraction of the mean NDVI for each point in the object eco and plot of the data:

ndvi <- extract(mean.ndvi, eco3[["XY"]])
ndvi<- aue.rescale(ndvi)
plot(eco3[["XY"]][, 1], eco3[["XY"]][, 2], col=rgb(ndvi, 0, 0),
pch=15, main = "Mean NDVI", xlab = "X", ylab  = "Y")
 
 
#---Theil-sen estimation for a raster---#

require("raster")
set.seed(6)

temp <- list()

for(i in 1:100) {
temp[[i]] <- runif(36,-1, 1)
temp[[i]] <- matrix(temp[[i]], 6, 6)
temp[[i]] <- raster::raster(temp[[i]])
}

temp <- brick(temp)

writeRaster(temp,"temporal.tif", overwrite=T)
rm(temp)
ndvisim <- brick("temporal.tif")
date <- seq(from = 1990.1, length.out = 100, by = 0.2)
eco.theilsen(ndvisim, date)
pvalue <- raster("pvalue.tif")
slope <- raster("slope.tif")
par(mfrow = c(1, 2))
plot(pvalue, main = "p-value")
plot(slope, main = "slope")


#---Conversion of ecogen data to several formats and format tool---#


## Interconverting an ecogen genetic data frame among several formats 

data(eco3)

### One allele per column
loc2al <- eco.convert(eco3[["G"]], "matrix", "alleles.matrix", ploidy = 2)
loc2al

### Inverse operation (collapse alleles into locus)
al2loc <- eco.convert(loc2al, "alleles.matrix", "matrix", ploidy = 2)
al2loc

### Separating alleles with a character string
loc2loc <- eco.convert(eco3[["G"]], "matrix", "matrix", ploidy = 2, sep.out = "/")
loc2loc

### Inverse operation (removing separator)
loc2loc.nosep <- eco.convert(loc2loc, "matrix", "matrix", ploidy = 2, sep.in = "/", sep.out = "")
loc2loc.nosep


## Format tool

data(eco.test)

### Adding zeros

example <- as.matrix(genotype[1:10,])
mode(example) <- "character"
### example data
example
recoded <- eco.format(example, ncod = 1, ploidy = 2, nout = 3)
### recoded data
recoded
### leading zeros
recoded2 <- eco.format(example, ncod = 1, ploidy = 2, nout = 3, 
fill.mode = "first")
### recoded data
recoded2


### Tetraploid data, separating alleles with a "/"
tetrap <- as.matrix(example)
### simulated tetraploid example data
tetrap <- matrix(paste(example,example, sep = ""), ncol = ncol(example))
recoded <- eco.format(tetrap, ncod = 1, ploidy = 4, sep.out = "/")
### recoded data
recoded


## Creating input data for Geneland with an ecogen object

data(eco.test)
eco.2geneland(eco, 1)


## Exporting an ecogen genetic data frame into Genepop format

data(eco.test)
eco.2genepop(eco, grp = "pop", name = "infile.genepop.txt")
### an output file "infile.genepop.txt" is generated in the working directory


## Converting a diploid ecogen genetic data frame into a gstudio object

data(eco.test)
gsteco <- eco.2gstudio(eco, "separated")
gsteco


## Converting an ecogen genetic data frame into a hierfstat data frame

data(eco.test)
hiereco <- eco.2hierfstat(eco, "pop")


## Exporting an ecogen genetic data frame into SPAGeDI format

data(eco.test)
eco.2spagedi(eco, "pop", ndig = 1,int=2, smax=6, name="infile.spagedi.txt")

}
}
