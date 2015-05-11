# Generating atmospherically corrected NDVI and MSAVI2 images for 
# temporal series of Landsat 5 and 7

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015

setGeneric("eco.NDVI", 
           function(tab, 
                    correct = c("COST", "DOS"), 
                    method = c("NDVI", "MSAVI2"), 
                    landsat = c("LT5", "LT7.L", "LT7.H"),
                    datatype =  c("FLT4S", "FLT8S", "INT4U", "INT4S", 
                                  "INT2U", "INT2S", "INT1U", "INT1S", 
                                  "LOG1S")) {
             
             correct <- match.arg(correct)
             method <- match.arg(method)
             landsat <- match.arg(landsat)
             datatype <- match.arg(datatype)
             
             
             steps <- nrow(tab)
             
             if(landsat == "LT5") {
               bresb3 <- -2.213976
               gresb3 <- 1.043976
               bresb4 <- -2.386024
               gresb4 <- 0.8760236
               esunb4 <- 1031
               esunb3 <- 1536
               
             } else if(landsat == "LT7.L") {
               
               bresb3 <- -5.9425197
               gresb3 <- 0.9425197
               bresb4 <- -6.0692913
               gresb4 <- 0.9692913
               esunb4 <- 1039
               esunb3 <- 1533
               
             } else if(landsat == "LT7.H") {
               
               bresb3 <- -5.6216535
               gresb3 <- 0.6216535
               bresb4 <- -5.7397638
               gresb4 <- 0.6397638
               esunb4 <- 1039
               esunb3 <- 1533
               
             }
             
             radiancia <- function(banda, brescale, grescale) {
               banda[] <- grescale * banda[] + brescale
               banda
             }
             
             reflectancia <- function(banda, edist, Esun, coselev, TAUz) {
               banda[] <-  (pi * edist ^ 2 * banda[]) / (Esun * coselev * TAUz)
               banda
             }
             
             Lhaze <- function(Lhazeban, grescale, brescale, Esun, coselev,
                               edist, TAUz) {
               Lhazeban <-  (grescale * Lhazeban + brescale) - 
                 0.01 * (Esun * coselev * TAUz) / (pi * edist ^ 2)
               Lhazeban
             }
             
             
             
             esperar <- function(i) {
               cat ("\r", ceiling(100 * i / steps), "% ", "completed", sep = "")
             }
             
             
             ESdist <- function (adate) {
               edist <- julian(as.Date(adate),
                               origin = as.Date(paste(substring(adate, 1, 4),
                                                      "12", "31", sep = "-")))[[1]]
               edist <- 1 - 0.016729 * cos((2 * pi) * (0.9856 * (edist - 4) / 360))
               edist
             }
             
             cat("0% completed")
             
             j <- 0
             
             for(i in 1:steps) {
               sunelev <- tab[i, 2]
               edist <- ESdist(as.character(tab[i, 1]))
               suntheta <- (90 - sunelev) * pi / 180           
               suntheta <- cos(suntheta)
               if(correct == "DOS") {
                 TAUz <- 1} else {
                   TAUz <- suntheta
                 }
               
               b3 <- raster::raster(as.character(tab[i, 4]))
               
               
               message("band 3 loaded")
               
               Lhazeb3 <- tab[i, 6]
               
               b3 <- radiancia(b3, bresb3, gresb3)
               Lhazeb3 <- Lhaze(Lhazeb3, gresb3, bresb3, esunb3, suntheta, edist, TAUz)
               b3[] <- b3[] - Lhazeb3
               b3 <- reflectancia(b3, edist, esunb3, suntheta, TAUz)
               b3[b3 < 0] <- 0
               
               j <- j + steps / 4
               esperar(j)
               
               
               b4 <- raster::raster(as.character(tab[i, 3]))
               
               
               message("band 4 loaded")
               Lhazeb4 <- tab[i, 5]
               b4 <- radiancia(b4, bresb4, gresb4)
               Lhazeb4 <- Lhaze(Lhazeb4, gresb4, bresb4, esunb4, suntheta, edist, TAUz)
               b4[] <- b4[] - Lhazeb4
               b4 <- reflectancia(b4, edist, esunb4, suntheta, TAUz)
               b4[b4 < 0] <- 0
               
               j <- j + steps / 4
               esperar(j) 
               
               
               if(raster::extent(b4) != raster::extent(b3)) {  
                 dimension <- raster::intersect(b4, b3)
                 b4 <- crop(b4, dimension)
                 b3 <- crop(b3, dimension)
               }
               
               
               if (method  ==  "NDVI") {
                 NDVI = (b4 - b3) / (b4 + b3)
                 NDVI[NDVI <- 1] <- NA
                 #NDVI=(NDVI+1)*255/2    #OPTIONAL FOR RESCALING TO 8 BITS
                 
                 message("writing NDVI image on disk")
                 
                 raster::writeRaster(NDVI, as.character(paste("NDVI", correct, tab[i, 7],
                                                              sep = "")), format = "GTiff",
                                     datatype = datatype, overwrite = T)
                 
               } else if (method  ==  "MSAVI2") {
                 MSAVI2 = ((2 *b4 + 1) - sqrt((2 * b4 + 1) ^ 2 - 8 * (b4 - b3))) / 2
                 MSAVI2[MSAVI2 <- 1] <- NA
                 #MSAVI2 <- (MSAVI2+1)*255/2       #OPTIONAL FOR RESCALING TO 8 BITS
                 message("writing MSAVI2 image on disk")
                 
                 raster::writeRaster(MSAVI2, as.character(paste("MSAVI2", correct,
                                                                tab[i, 7], sep="")),
                                     format="GTiff", datatype = datatype, overwrite =T)
               }
               
             }
             cat("done!")
           })
