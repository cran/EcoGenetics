# Postprocessing for NDVI and MSAVI 2 temporal series of Landsat 5 and 7
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015


setGeneric("eco.NDVI.post", 
           function(tab,correct = c("COST", "DOS"), 
                    method = c("NDVI", "MSAVI2"), 
                    datatype =  c("FLT4S", "FLT8S", "INT4U", "INT4S", 
                                  "INT2U", "INT2S", "INT1U", "INT1S", 
                                  "LOG1S"), 
                    what = c("mean", "max", "min", "var", "none")) {
             
             correct <- match.arg(correct)              
             method <- match.arg(method)
             datatype <- match.arg(datatype)
             what <- match.arg(what, several.ok =  TRUE)
             
             y <- pmatch(what, c("mean", "max", "min", "var", "none"))
             
             
             esperar <- function(i){
               cat ("\r", 100 * i / steps, "% ", "complete",  sep = "")
             }
             
             
             steps <- nrow(tab)
             
             message("\n", "loading parameters", "\n")
             
             xmin <- rep(0, nrow(tab))
             xmax <- rep(0, nrow(tab))
             ymin <- rep(0, nrow(tab))
             ymax <- rep(0, nrow(tab))
             
             for(i in 1:nrow(tab)) {
               xmin[i] <- raster::extent(raster::raster(as.character(paste(method, correct,
                                                                           tab[i, 7], ".tif",
                                                                           sep = ""))))@xmin
               xmax[i] <- raster::extent(raster::raster(as.character(paste(method, correct,
                                                                           tab[i, 7], ".tif",
                                                                           sep = ""))))@xmax
               ymin[i] <- raster::extent(raster::raster(as.character(paste(method, correct,
                                                                           tab[i, 7], ".tif",
                                                                           sep = ""))))@ymin
               ymax[i] <- raster::extent(raster::raster(as.character(paste(method, correct,
                                                                           tab[i, 7], ".tif",
                                                                           sep = ""))))@ymax
             }
             
             dimension = c(max(xmin), min(xmax), max(ymin), min(ymax))
             
             message("\n", "stacking images, please wait...", "\n")
             
             raster::writeRaster(raster::crop(raster::raster(as.character(paste(method, correct,
                                                                                tab[1, 7], ".tif",
                                                                                sep = ""))), 
                                              raster::extent(dimension)), paste(method, correct, "time.tif",
                                                                                sep = ""), format = "GTiff",
                                 datatype = datatype, overwrite = T)
             
             
             for(i in 2:nrow(tab)) {
               raster::writeRaster(raster::addLayer(raster::brick(paste(method, correct, "time.tif",
                                                                        sep = "")), 
                                                    raster::crop(raster::raster(as.character(paste(method, correct,
                                                                                                   tab[i, 7], 
                                                                                                   ".tif", sep = ""))),
                                                                 raster::extent(dimension))), format = "GTiff",
                                   datatype = datatype,  "temporal.tif", 
                                   overwrite = T)    
               raster::writeRaster(raster::brick("temporal.tif"), paste(method, correct, "time.tif", 
                                                                        sep = ""), format = "GTiff",
                                   datatype = datatype,  overwrite = T)
               esperar(i)
               
             }
             file.remove("temporal.tif")
             
             
             message("computing...")
             
             time <- raster::brick(paste(method, correct, "time.tif",  sep = ""))
             
             if(any(y == 1)) {
               cat("\n", "computing mean", "\n")
               temp <- raster::calc(time, mean)
               raster::writeRaster(temp, paste(method, ".", correct, ".", "mean.tif",
                                               sep = ""), format="GTiff", 
                                   datatype = datatype,  overwrite = T)
             }
             
             if(any(y == 2)) {
               cat("\n", "computing max", "\n")
               temp <- raster::calc(time, max)
               raster::writeRaster(temp, paste(method, ".", correct, ".", "max.tif",
                                               sep = ""), format="GTiff",
                                   datatype = datatype, overwrite = T)
             }
             
             if(any(y == 3)) {
               cat("\n", "computing min", "\n") 
               temp <- raster::calc(time, min)
               raster::writeRaster(temp, paste(method, ".", correct, ".", "min.tif",
                                               sep = ""), format="GTiff", 
                                   datatype = datatype, overwrite = T)
             }
             
             
             if(any(y == 4)) {
               cat("\n", "computing variance", "\n")
               temp <- raster::calc(time, var)
               raster::writeRaster(temp, paste(method, ".", correct, ".", "var.tif",
                                               sep = ""), format = "GTiff", 
                                   datatype = datatype, overwrite = T)
             }
             
             cat("\n", "done", "\n")
             ReturnVal <- tcltk::tkmessageBox(title = "VI post process", 
                                              message = "process successful!",
                                              icon = "info", type = "ok")
           })
