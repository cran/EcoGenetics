# Leandro Roser leandroroser@ege.fcen.uba.ar
# June 17, 2015 


# Theil-sen regression for a raster time series

setGeneric("eco.theilsen", 
           function(stacked, date, 
                    adjust = "none") {
             
             
             adjust <- match.arg(adjust)
             
             esperar <- function(i)
             {
               cat ("\r", ceiling(100 * i / steps), "% ",
                    "completed", sep = "")
             }
             
             
             cat("starting...", "\n\n")
             
             fun <- function(date, data)
             {
               mod <- rkt::rkt(date, data)
               return(c(as.numeric(mod[3]), as.numeric(mod[1])))
             }
             
             cat("pre-processing data...", "\n\n")
             
             pendiente <- rep(NA, raster::ncell(stacked))
             pvalor <- rep(NA, raster::ncell(stacked))
             df <- as.data.frame(as.matrix(stacked))
             steps = nrow(df)
             
             for(i in 1:nrow(df)) {
               resultados <- fun(date,as.numeric(df[i, ]))
               pendiente[i] <- resultados[1]
               pvalor[i] <- resultados[2]
               esperar(i)
             }
             cat("\n\n")
             
             
             if(adjust != "none") {
               cat(paste("adjusting p values with", adjust, "method"), "\n")
               pvalor <- p.adjust(pvalor, "adjust")
             }
             
             cat("writing slope image to workspace...", "\n\n")
             
             pendiente <- matrix(pendiente, nrow = stacked@nrows,
                                 ncol = stacked@ncols, byrow = TRUE)
             pendiente <- raster::raster(pendiente, crs = stacked@crs,
                                         xmn = stacked@extent@xmin,
                                         ymn = stacked@extent@ymin, 
                                         xmx = stacked@extent@xmax,
                                         ymx = stacked@extent@ymax)
             raster::writeRaster(pendiente, "slope.tif", overwrite = T)
             
             cat("writing P-value image to workspace...", "\n")
             
             pvalor <- matrix(pvalor, nrow=stacked@nrows, ncol = stacked@ncols,
                              byrow = TRUE)
             pvalor <- raster::raster(pvalor, crs = stacked@crs,
                                      xmn = stacked@extent@xmin,
                                      ymn = stacked@extent@ymin, 
                                      xmx = stacked@extent@xmax,
                                      ymx = stacked@extent@ymax)
             
             cat("\n","done!","\n\n" )
             
             
             raster::writeRaster(pvalor, "pvalue.tif", overwrite = T)
             
             
           })
