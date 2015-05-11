\name{EcoGenetics-package}
\alias{EcoGenetics-package}
\alias{EcoGenetics}
\docType{package}
\title{
 Analysis of Phenotypic, Genotypic and Environmental Data
}
\description{
This package allows users to analyze phenotypic, genotypic and environmental data with a single object of class "ecogen". EcoGenetics has been designed to make easy the storing, handling and integration of the available information from multiple sources. It also provides geostatistical tools for analyzing spatial patterns, and extra functions. 
}

\details{
 \tabular{ll}{
 Package: \tab EcoGenetics\cr
 Type: \tab Package\cr
 Version: \tab 1.2.0\cr
 Date: \tab 2015-05-11\cr
 License: \tab GPL (>=2) \cr
 }
 
 ----------------------- o
 -----------------------
 
 \bold{\emph{Structure of an ecogen object: handling and integrating information}}


Ecological genetics research requires the integration of data originated in different sources. The object class "ecogen" has been designed for handling multidimensional data. Its basic structure is the following:

- A \strong{P} slot, storing a phenotypic data frame.  

- A \strong{G} slot, storing a genotypic data frame.  

- An \strong{E} slot, storing an environmental data frame.  

- A \strong{S} slot, storing a data frame with classes assigned to the 
individuals. 

- A \strong{GENIND} slot,containing in \emph{genind} format the information of @G.

- A \strong{C} slot, for a custom data frame. 

- An \strong{OUT} slot, containing a list for the storage of the results.


The construction of a new "ecogen" object is made with the homonymous function.This object class allows the importation from different genetic formats (codominant-dominant data of class "data.frame", "DNAbin" and "allignment"):

\emph{library("EcoGenetics")}

\emph{data(eco.test)}

\emph{# "data.frame" data}

\emph{eco <- ecogen(XY = coordinates, P = phenotype, G = genotype, E = environment, S = structure)}

\emph{# "DNAbin" data}

\emph{require("ape")}

\emph{data(woodmouse)}

\emph{G <- woodmouse}

\emph{eco.dnabin <- ecogen(G = G)}

\emph{# "allignment" data}

\emph{require("seqinr")}

\emph{data(mase)}

\emph{G <- mase}

\emph{eco.allign <- ecogen( G = G)}


Continuing with the first example, the following available methods can be used with ecogen objects:

- - - - - - 

\bold{# 1. subsetting-method}

\emph{eco.sub <- eco[1:50]}

\bold{# 2. merging-method, for two objects}

\emph{eco1 <- eco}

\emph{merged <- eco.merge(eco, eco1)}

\bold{# 3. subsetting-method, in reference to a group in the S slot (in this case, "1")}

\emph{eco.subS <- eco.subset(eco,"pop", 1)}


\bold{# 4. binding by row-method}

\emph{eco.r <- eco.rbind(eco, eco1)}

\bold{# 5. binding by column-method}

\emph{eco.c <- eco.cbind(eco, eco1)}

\bold{# 6. ordering by row-method, using the rows in the XY data frame as reference}

\emph{eco1$P <- eco$P[sample(1:173), ] ## object with unordered rows in P data frame}

\emph{ordered <- eco.order(eco1)}

\bold{# 7. appending-method (storing information generated)}

\emph{## fitting a multiple linear regression model:}
 
\emph{linear.analysis <- eco.lmtree(eco$P, eco$E, "mlm")}
 
\emph{## storing the results:}
 
\emph{eco <- eco.append(eco, linear.analysis)}
 
\emph{eco}

\emph{## Storing multiple result at once:}
 
\emph{a <- c(1:10)}

\emph{b <- c(2:30)}
  
\emph{eco <- eco.append(eco, a, b)}

\emph{eco}
  
\emph{# the panel shows the object with the results in alphabetical order and their classes in parentheses}

\bold{# 8. removing-method}

\emph{eco <- eco.remove(eco, a, b) #removing objects a and b eco}

- - - - - - 

Finally, the workspace can be cleared, only storing the ecogen object:

\emph{eco.clear(eco)}

 ----------------------- o
 -----------------------
 
 \bold{\emph{Structure of the package}}

The package has four basic modules. The \strong{base module} is composed by general functions (multiple lm, conversion utilities, etc.). 
The \strong{general spatial module} computes global (Moran's I, Mantel test, etc.) and local (Getis-Ord's G*, local Moran's I, etc.) spatial tests. These analyses use a spatial weights matrix, provided by the  \strong{spatial weights module}.  The \strong{lag analysis module} performs two basic analyses: the obtention  of variograms and correlograms. This module uses the other tool provided by the spatial weights module: a sequence of spatial weights matrices. 
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

Dray S., and A., Dufour. 2007. The ade4 package: implementing the duality diagram for ecologists. Journal of statistical software, 22: 1-20.
 
Freedman D., and P. Diaconis. 1981. On the histogram as a density estimator: L 2 theory. Probability theory and related fields, 57: 453-476.
  
Geary R. 1954. The contiguity ratio and statistical mapping. The incorporated statistician, 115-146.

Getis A., and J. Ord. 1992. The analysis of spatial association by use of distance statistics. Geographical analysis, 24: 189-206.

Goslee S. 2011. Analyzing remote sensing data in R: the landsat package. Journal of Statistical Software, 43: 1-25.

Goudet J. 2005. Hierfstat, a package for R to compute and test hierarchical F-statistics. Molecular Ecology Notes, 5: 184-186.
 
Guillot G., F. Mortier and A. Estoup. 2005. GENELAND: a computer package for landscape genetics. Molecular Ecology Notes, 5: 712-715.
 
Jombart T. 2008. adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24: 1403-1405.

Legendre P., and L. Legendre. 2012. Numerical ecology. Third English edition .Elsevier Science, Amsterdam, Netherlands.

Lichstein J., T. Simons, S. Shriner, and K. Franzreb. 2002. Spatial autocorrelation and autoregressive models in ecology. Ecological monographs, 72: 445-463.

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

Wu. 1986. Jackknife, bootstrap and other resampling methods in regression analysis. the Annals of Statistics, 1261-1295.
 
}
 
\examples{

 \dontrun{
 
 
#---Detrending spatial data with polynomial interpolation---#
 
 
data(eco2)

## original data
data1 <- matrix(eco2$P[,1], 30, 30)
image(data1)

## original data + trend
data2 <- matrix(eco2$P[,2], 30, 30)
image(data2)

## data detrending
data2.det <- eco.detrend(Z = eco2$P[,2], XY =  eco2$XY, degree =  1)
data2.det <- data2.det@RES
data2.det <- matrix(data2.det$df1, 30, 30)
image(data2.det)


#---Multiple Linear Regression fit---#
 
 
data(eco.test)
mymod <- "E1+E2*E3"
mod <- eco.lmtree(df1 = eco$P, df2 = eco$E, 
analysis = "mlm", mod.class = mymod) 
summary(mod)
 
 
#---Multiple Conditional Inference Trees---#


data(eco.test)
mymod <- "E1+E2*E3"
mod <- eco.lmtree(df1 = eco$P, df2 = eco$E, 
analysis = "mctree", mod.class = mymod, fact = eco$S$structure)
summary(mod)


#---Global spatial analysis---#
 
 
## Moran's I 

### one test
data(eco.test)
con <- eco.weight(eco$XY, method = "circle", d1 = 0, d2 = 2)
global <- eco.gsa(Z = eco$P[, 1], con = con, , method = "I", nsim = 200)
global

require(adegenet)
con2<-chooseCN(eco$XY, type = 1, result.type = "listw", plot.nb = FALSE)
global <- eco.gsa(Z = eco$P[, 1], con = con2, , method = "I", nsim = 200)
global

### multiple tests
con <- eco.weight(eco$XY, method = "circle", d1 = 0, d2 = 2)
global <- eco.gsa(Z = eco$P, con = con, , method = "I", nsim = 200)
global


## Geary's C 

data(eco.test)
global.C <- eco.gsa(Z = eco$P[, 1], con = con, method = "C", nsim = 200)
global.C


## Bivariate's Moran's Ixy

data(eco.test)
global.Ixy <- eco.gsa(Z = eco$P[, 1], Y = eco$E[, 1],
con = con, method = "CC", nsim = 200)
global.Ixy


## Join-Count

data(eco.test)
global.JC <- eco.gsa(Z = 2* eco$GENIND$tab[, 1], ncod = 1,
con = con, method = "JC", nsim = 10)
global.JC


# Mantel test

data(eco.test)
eco.mantel(d1 = dist(eco$P), d2 = dist(eco$E), nsim = 99)  


## Partial Mantel test

data(eco.test)
eco.mantel(d1 = dist(eco$P), d2 = dist(eco$E),
dc = dist(eco$XY), nsim = 99)  


#---Local spatial analysis---#


## Getis-Ord's G*

data(eco.test)
require(ggplot2)
con<- eco.weight(eco$XY, method = "knearest",  k = 4, self = TRUE) 
### self = TRUE for G*
getis.ak <- eco.lsa(eco$P[, 1], con, method = "G*", nsim = 99, 
adjust = "none")
getis.ak

### to plot the results, the function "eco.lsa" calls "eco.rankplot"
(see ?eco.rankplot) when test = "permutation" and "eco.forestplot"
(see ?eco.forestplot) when test = "bootstrap"

p <- plot(getis.ak)   ###  rankplot graph
p                     ###  points with colors of the color-scale:
                      ###  points with P < 0.05. Yellow points : 
                      ###  points with P > 0.05
p <- plot(getis.ak, significant = FALSE)
p                     ###  all points have a color of the color-scale

### bootstrap example
getis.akb <- eco.lsa(eco$P[, 1], con, method = "G*", nsim = 99,
test = "bootstrap")
p <- plot(getis.akb)     ### forestplot graph
p + ggplot2::theme_bw()  ### the plot can be modified with ggplot2
                         ### In this case, the background color is modified 


## Getis-Ord's G

data(eco.test)
require(ggplot2)
con<- eco.weight(eco$XY, method = "knearest", k = 4)
### self = FALSE for G
getis <- eco.lsa(eco$P[, 1], con, method = "G", nsim = 99, adjust = "none")
plot(getis)


## Local Moran's I

data(eco.test)
require(ggplot2)
con<- eco.weight(eco$XY, method = "knearest",  k = 4, row.sd = TRUE)
### row standardized weights = TRUE
localmoran <- eco.lsa(eco$P[, 1], con, method = "I", nsim = 99, 
adjust = "none")
plot(localmoran)


## Local Geary's C

data(eco.test)
require(ggplot2)
con<- eco.weight(eco$XY, method = "knearest",  k = 4, row.sd = TRUE)
### row standardized weights = TRUE
localgeary <- eco.lsa(eco$P[, 1], con, method = "C", nsim = 99, 
adjust = "none")
plot(localgeary)


#---Moran's I, Geary's C and bivariate Moran's Ixy correlograms---#


## Moran's I correlogram

### single test
data(eco.test)
require(ggplot2)
moran <- eco.correlog(Z=eco$P[,1], XY = eco$XY, method = "I",
smax=10, size=1000)
plot(moran)

### multiple tests
moran2 <- eco.correlog(Z=eco$P, XY = eco$XY, method = "I",
smax=10, size=1000)
plot(moran2, var ="P2")
plot(moran2, var ="P3")

### correlogram plots support the use of ggplot2 syntax
moranplot <- plot(moran2, var ="P3") + theme_bw() + theme(legend.position="none")
moranplot


## Geary's C correlogram

data(eco.test)
require(ggplot2)
geary <- eco.correlog(Z=eco$P[,1], XY = eco$XY, method = "C",
smax=10, size=1000)
plot(geary)


## Moran's Ixy cross-correlogram

data(eco.test)
require(ggplot2)
cross <- eco.correlog(Z=eco$P[,1], XY = eco$XY, Y = eco$E[, 1],
method = "CC", int= 2, smax=15)
plot(cross)


#---Mantel and partial Mantel correlograms---#


data(eco.test)
require(ggplot2)
corm <- eco.cormantel(M = dist(eco$P), size=1000,smax=7, XY = eco$XY,
nsim = 99)
plot(corm)

corm <- eco.cormantel(M = dist(eco$P), size=1000,smax=7, XY = eco$XY,
nsim = 99, test = "bootstrap")
plot(corm)


## partial Mantel correlogram

corm <- eco.cormantel(M = dist(eco$P), MC = dist(eco$E),
size=1000, smax=7, XY = eco$XY, nsim = 99)
plot(corm)

### correlogram plots support the use of ggplot2 syntax
mantelplot <- plot(corm) + theme_bw() + theme(legend.position="none")
mantelplot


#---Empirical variogram---#


data(eco.test)
require(ggplot2)
variog <- eco.variogram(Z = eco$P[, 2],XY =  eco$XY)
plot(variog)

## variogram plots support the use of ggplot2 syntax
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

ndvi <- extract(mean.ndvi, eco3$XY)
ndvi<- aue.rescale(ndvi)
plot(eco3$XY[, 1], eco3$XY[, 2], col=rgb(ndvi, 0, 0),
pch=15, main = "Mean NDVI", xlab = "X", ylab  = "Y")
 
 
#---Theil-sen estimation for a raster---#

require("raster")
require("animation")
set.seed(6)

temp <- list()

for(i in 1:100) {
temp[[i]] <- runif(36,-1, 1)
temp[[i]] <- matrix(temp[[i]], 6, 6)
temp[[i]] <- raster::raster(temp[[i]])
}

temp <- brick(temp)

oopt <- ani.options(interval = 0.01)
for (i in 1:ani.options("nmax")) {
plot(temp[[i]])
ani.pause()
}
ani.options(oopt)

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


#---Conversion of ecogen data to several formats---#


## Converting a diploid ecogen genetic data frame with two alleles 
## per column into a data frame with one allele per column

data(eco3)
eco.2al <- eco.2columns(eco3, 1)
eco.2al


## Creating input data for Geneland with an ecogen object

data(eco.test)
eco.2geneland(eco, 1)


## Exporting an ecogen genetic data frame into Genepop format

data(eco.test)
eco.2genepop(eco, grp = "pop", ndig = 1, name = "infile.genepop.txt")
### an output file "infile.genepop.txt" is generated in the working directory


## Converting a diploid ecogen genetic data frame into a gstudio object

data(eco.test)
gsteco <- eco.2gstudio(eco, "separated")
gsteco


## Converting an ecogen genetic data frame into a hierfstat data frame

data(eco.test)
hiereco <- eco.2hierfstat(eco, "pop")
require("hierfstat")
basic.stats(hiereco)


## Exporting an ecogen genetic data frame into SPAGeDI format

data(eco.test)
eco.2spagedi(eco, "pop", ndig = 1,int=2, smax=6, name="infile.spagedi.txt")

}
}
