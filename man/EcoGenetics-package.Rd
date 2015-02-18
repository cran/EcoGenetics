\name{EcoGenetics-package}
\alias{EcoGenetics-package}
\alias{EcoGenetics}
\docType{package}
\title{
 Analysis of Phenotypic, Genotypic and Environmental Data
}
\description{
 Analysis of phenotypic, genotypic and environmental data 
 with an object of class "ecogen", designed to make easy the integration
 of the available information from multiple sources. Also including 
 geostatistical tools and extra functions.
}

\details{
 \tabular{ll}{
 Package: \tab EcoGenetics\cr
 Type: \tab Package\cr
 Version: \tab 1.1.0\cr
 Date: \tab 2015-02-18\cr
 License: \tab GPL (>=2) \cr
 }
 
In the ecogenetics discipline, when working with a set of individuals, we have different kind of data that we wish to integrate and analyze. The ecogen class is a powerful unit of information storing and processing, designed to make comfortable and faster that work.

The basic structure of an ecogen object is the following:

\enumerate{
\item A \strong{P} slot, storing a phenotypic data frame.
\item A \strong{G} slot, storing a genotypic data frame.
\item An \strong{E} slot, storing an environmental data frame.
\item A \strong{S} slot, storing a data frame with classes assigned to the individuals. For example, an S data frame could consist of a first column with demes and a second with region information, both numeric. The only requirement that this data frame must fulfil, is to provide classification rules. 
\item A \strong{GENIND} slot, with a genind object, containing in genind
format the G information.
\item A \strong{C} slot, for a custom data frame. For example, it could be epigenetic information or other class of G data that should be stored apart.
\item An \strong{OUT} slot, containing a list for the storage of the results. 
}

Each time we wish to store new results, we use the function \code{\link{eco.append}}, assigning the new result to our ecogen object. In the examples below we explain better how this work. First, with the data preloaded into our workspace, we create a new ecogen object. There are two option. We create an empty object and then we add the data,  or we add all the data at once:

In the first one option:

\strong{empty <- ecogen()}

\strong{empty}


Next, we include into our object, for example, the P data frame:

\strong{empty$P <- phenotype}


In the second one option:

\strong{eco <- ecogen(XY = coordinates, P = phenotype, G = genotype, E = environment, S = as.data.frame(structure))}


We are free of including the data desired when needed.

With all the information stored in our ecogen object, we clear out the workspace with the function \code{\link{eco.clear}}, removing all the older data:

\strong{eco.clear(eco)}


Then, we can for example merge two objects of this class.  For simplicity, we create a new object shuffling eco rows:

 \strong{eco1 <- eco}
 
 \strong{eco1$XY <- eco$XY[sample(1:173), ]}
 
 
We merge both data sets:

 \strong{merged <- eco.merge(eco, eco1)}
 
 
New object with the data of both eco and eco1:

 \strong{merged} 
 
Subsetting by rows, individuals 1 to 50: 

\strong{merged[1:50]}

Subsetting by groups, in this case by group 1:

\strong{eco.subset(eco,"structure", 1)}


We can analyze the data and store the new information generated. For example, we fit a multiple linear regression model:
 
 \strong{linear.analysis <- eco.lmtree(merged$P, merged$E, "mlm")}
 
 \strong{ linear.analysis}
 
 \strong{summary(linear.analysis)}
 
 
 We store the results:
 
 \strong{merged <- eco.append(merged, linear.analysis)} 
 
 \strong{merged} 
 
 
 We store multiple result at once:
 
  \strong{a <- c(1:10)}
  
  \strong{b <- c(2:30)}
  
  \strong{eco <- eco.append(eco, a, b)} 
  
  
As it can be seen, the panel show us the object with the results and their classes in parentheses. Always the results are stored in alphabetical order for easy seach. If needed, as $OUT is a list we may write                      \emph{eco$OUT$linear.analysis} to recover this data.
 
 We can also remove results from our object:
 
 \strong{eco <- eco.remove(eco, a)}
 
 
 
 Finally, we clear out the workspace, only storing the ecogen object:
 
 \strong{eco.clear(eco)}    
 
 
 }
 
\author{

 Leandro Roser, Juan Vilardi, Beatriz Saidman and Laura Ferreyra
 
 Maintainer: Leandro Roser <leandroroser@ege.fcen.uba.ar>
}

\references{

Anselin, L. (1995). Local indicators of spatial association - LISA. Geographical analysis, 27(2), 93-115.

 Chander, G., Markham, B. L., & Helder, D. L. (2009). Summary of current radiometric calibration coefficients for Landsat MSS, TM, ETM+, and EO-1 ALI sensors. Remote sensing of environment, 113(5), 893-903.
 
 Chavez, P. S. (1996). Image-based atmospheric corrections-revisited and improved. Photogrammetric engineering and remote sensing, 62(9), 1025-1035.
 
 Cliff, A. D., & Ord, J. K. (1981). Spatial processes: models & applications (Vol. 44). London: Pion.
 
 Dray, S., & Dufour, A. B. (2007). The ade4 package: implementing the duality diagram for ecologists. Journal of statistical software, 22(4), 1-20.
 
 Goudet, J. (2005). Hierfstat, a package for R to compute and test hierarchical F-statistics. Molecular Ecology Notes, 5(1), 184-186.
 
 Guillot, G., Mortier, F., & Estoup, A. (2005). GENELAND: a computer package for landscape genetics. Molecular Ecology Notes, 5(3), 712-715.
 
 Hothorn, T., Hornik, K., & Zeileis, A. (2006). Unbiased recursive partitioning: A conditional inference framework. Journal of Computational and Graphical statistics, 15(3), 651-674.
 
 Jombart, T. (2008). adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24(11), 1403-1405.
 
  Legendre, P., & Fortin, M. J. (1989). Spatial pattern and ecological analysis. Vegetatio, 80(2), 107-138.
  
 Legendre, P., & Legendre, L. F. (2012). Numerical ecology (Vol. 24). Elsevier.
 
 Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. Cancer research, 27(2 Part 1), 209-220.
 
 Moran, P. A. (1948). The interpretation of statistical maps. Journal of the Royal Statistical Society. Series B (Methodological), 10(2), 243-251.
 
 Nei, M. (1978). Estimation of average heterozygosity and genetic distance from a small number of individuals. Genetics, 89(3), 583-590.
 
 Ord, J. K., & Getis, A. (1995). Local spatial autocorrelation statistics: distributional issues and an application. Geographical analysis, 27(4), 286-306.
 
 Pebesma, E. J. (2004). Multivariable geostatistics in S: the gstat package. Computers & Geosciences, 30(7), 683-691.
 
 Petit, R. J., & Pons, O. (1998). Bootstrap variance of diversity and differentiation estimators in a subdivided population. Heredity, 80(1), 56-61.
 
 Sen, P. K. (1968). Estimates of the regression coefficient based on Kendalls tau. Journal of the American Statistical Association, 63(324), 1379-1389.
 
Sokal, R. R., & Oden, N. L. (1978). Spatial autocorrelation in biology: 1. Methodology. Biological journal of the Linnean Society, 10(2), 199-228.

 Tucker, C. J. (1979). Red and photographic infrared linear combinations for monitoring vegetation. Remote sensing of Environment, 8(2), 127-150.

 Wu, C. F. J. (1986). Jackknife, bootstrap and other resampling methods in regression analysis. the Annals of Statistics, 1261-1295.

}
 
\examples{

 \dontrun{
 
 \# Creating an empty ecogen object:
 
 empty <- ecogen()
 
 
 \# Creating a filled ecogen object with microsatellite data:
 
 data(eco.test)
 eco <- ecogen(XY = coordinates, P = phenotype, G = genotype, 
 E = environment, S = as.data.frame(structure))
 eco <- eco.sortalleles(eco, 1)           \# ordering the alleles in a codominant marker
 
 
 \# Creating an ecogen object with "DNAbin" imported data:
 
 library("ape")
 data(woodmouse)
 G <- woodmouse
 eco <- ecogen( G = G, format = "DNAbin" ,missing = 0)
 
 
 \# Creating an ecogen object with "allignment" imported data:
 
 library("seqinr")
 data(mase)
 G <- mase
 eco <- ecogen( G = G, format = "alignment" ,missing = 0)
 
 
 \# Operations over ecogen objects:
 
 data(eco.test)
 
 eco2 <- eco.cbind(eco, eco, "P")         \# binding rows
 
 eco2 <- eco.rbind(eco, eco)              \#binding columns
 
 eco2 <- eco[1:20]                        \# subsetting from rows 1 to 20

 eco2 <- eco.subset(eco, "structure", 1)  \# subsetting by groups
 

 variog <- eco.variogram(eco$P[, 2], eco$XY, int=100, smax=1000)
 
 eco <- eco.append(eco, variog)           \#updating by adding results
 eco <- eco.append(eco, variog)           \#updating by removing results
 
 \# Ordering rows 
 
 data(eco.test)
 eco1 <- eco
 eco1$P <- eco$P[sample(1:173), ]         \#object with shuffled rows
 eco1$E <- eco$E[sample(1:173), ]
 ordered <- eco.order(eco1)
 ordered$XY; ordered$P; ordered$E; ordered$S
 
\# Merging two ecogen objects:

 data(eco.test)
 eco1 <- eco
 eco1$XY <- eco$XY[sample(1:173), ]       \# object with shuffled rows
 eco1$P <- eco$P[sample(1:173), ]
 eco$XY
 eco1$XY
 merged <- eco.merge(eco, eco1)
 eco2 <- eco.cbind(eco, eco, "ALL")
 eco2$XY
 
 
 \# Multiple Linear Regression fit:
 
 data(eco.test)
 mymod <- "E1+E2*E3"
 mod <- eco.lmtree(df1 = eco$P, df2 = eco$E, 
 analysis = "mlm", mod.class = mymod) 
 summary(mod)
 
 
 \# Multiple Conditional Inference Trees:

data(eco.test)
mymod <- "E1+E2*E3"
mod <- eco.lmtree(df1 = eco$P, df2 = eco$E, 
analysis = "ctree", mod.class = mymod, fact = eco$S$structure)
summary(mod)

 
 \# Bootstrap estimates of the effective 
 \# number of alleles per locus and group, 
 \# and harmonic overall mean:
 
 data(eco.test)
 Ae.raref <- eco.rarefact(eco, "structure", nrep = 50)
 
 
\# Bootstrap estimates of heterozygocity:

 data(eco.test)
 boot.test <- eco.boothet(eco, "structure", nrep = 20, pres = T)
 
 
\# Moran I statistic with Monte-Carlo test:

data(eco.test)
require(adegenet)
con <- chooseCN(eco$XY, type = 1, result.type = "listw",
plot.nb = FALSE)
eco.moran(eco$P[, 3], con, nsim = 99)


\# Geary C statistic with Monte-Carlo test:

data(eco.test)
require(adegenet)
con<-chooseCN(eco$XY, type = 1, result.type = "listw", 
plot.nb = FALSE)
eco.geary(eco$P[, 3], con, nsim = 99)


\# Join-Count statistic with Monte-Carlo test:

data(eco.test)
require(adegenet)
con<-chooseCN(eco$XY, type = 1, result.type = "listw",
plot.nb = FALSE)
eco.joincount(2 * eco$GENIND$tab[, 3], con,
nsim = 99, ncod = 1)


\# Local Getis-Ord G and G* statistics with Monte-Carlo test:

data(eco.test)
require(adegenet)
con<-chooseCN(eco$XY, type = 1, result.type = "listw", plot.nb = FALSE)
eco.getisord(eco$P[, 3], con, nsim = 99,
test = "permutation", classG = "G")
getis <- eco.getisord(eco$P[, 3], con, nsim = 99, classG = "G*")
eco.rankplot(getis, eco$XY)
getis <- eco.getisord(eco$P[, 3], con, nsim = 99, classG = "G*",
test = "bootstrap")
eco.forestplot(getis)

\# Local Moran I test

data(eco.test)
require(adegenet)
con <- chooseCN(eco$XY, type = 1, result.type = "listw", plot.nb = FALSE)
local <- eco.localmoran(eco$P[, 2], con, nsim = 99)
eco.rankplot(local, eco$XY)
local <- eco.localmoran(eco$P[, 2], con, nsim = 99, test = "bootstrap")
eco.forestplot(local)


\# Mantel test:

data(eco.test)
eco.mantel(dist(eco$P), dist(eco$XY), 99)


\# Mantel correlogram:

data(eco.test)
corm <- eco.cormantel(z = eco$P, xy = eco$XY,
nsim = 99,int=50, smax=1000, test = "permutation")
plot(corm)
corm <- eco.cormantel(z = eco$P, xy = eco$XY,
nsim = 99,int=50, smax=1000, test = "bootstrap")
plot(corm)

 
\# Moran I, Geary C, Join-Count, local Getis-Ord G and local Moran I correlograms   

 eco.ac <- eco.autocor(eco, int = 200, smax = 2000,
df = "P", select = "moran", nsim = 10)
plot(eco.ac, var = "P2")

eco.ac <- eco.autocor(eco, int = 200, smax = 2000,
df = "P", select = "geary", nsim = 10)
plot(eco.ac, var = "P2")

eco.ac <- eco.autocor(eco, int = 200, smax = 2000,
df = "GENIND", select = "joincount", nsim = 10,
indvar="L1.1")
plot(eco.ac)

eco.ac <- eco.autocor(eco, int = 200, smax = 2000,
df = "P", select = "getisord", nsim = 10)
plot(eco.ac, var = "P2")
eco.rankplot(eco.ac, eco$XY, "P2")
eco.forestplot(eco.ac, "P2")

eco.ac <- eco.autocor(eco, int = 200, smax = 2000,
df = "P", select = "localmoran", nsim = 10)
plot(eco.ac, var = "P2")
eco.rankplot(eco.ac, eco$XY, "P2")
eco.forestplot(eco.ac, "P2")


\# Semivariogram for a selected variable:

data(eco.test)
variog <- eco.variogram(eco$P[, 2], eco$XY, int=100, smax=1000)
plot(variog)
variog <- eco.variogram(eco$P[, 2], eco$XY, int=100, smax=1000, nsim = 99)
plot(variog)

 
\# Computing NDVI with atmospheric correction 
\# over a time series and estimation of the temporal mean:

 data(eco.test)
 temp <- list()

##### we create 4 simulated rasters for the data included in the object tab:

 for(i in 1:4) {
 temp[[i]] <- runif(19800, 0, 254)
 temp[[i]] <- matrix(temp[[i]], 180, 110)
 temp[[i]] <- raster::raster(temp[[i]])
 extent(temp[[i]]) <- c(3770000, 3950000,6810000, 6920000)
 }

 
 writeRaster(temp[[1]], "20040719b4.tif", overwrite = T) 
 writeRaster(temp[[2]], "20040719b3.tif", overwrite = T)
 writeRaster(temp[[3]], "20091106b4.tif", overwrite = T)
 writeRaster(temp[[4]], "20091106b3.tif", overwrite = T)
 
 #### Computing NDVI images: 

 eco.NDVI(tab, "COST", "NDVI", "LT5")
 
 #### Mean NDVI image computed over the NDVI images that we created:
 
 eco.NDVI.post(tab, "COST", "NDVI", what = "mean")
 mean.ndvi <- raster("NDVI.COST.mean.tif")
 plot(mean.ndvi)
 
 #### Extraction of the mean NDVI for each point in the object eco and plot of the data:
 
 ndvi <- extract(mean.ndvi, eco$XY)
 ndvi<-eco.rescale(ndvi)
 plot(eco$XY[, 1], eco$XY[, 2], col=rgb(ndvi, 0, 0),
 pch=15, main = "Mean NDVI", xlab = "X", ylab  = "Y")
 
 
 \# Theil-sen estimation for a raster:
 
require("raster")
require("animation")

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

}
}
