#' Theil-sen regression for a raster time series.
#' @param stacked Stacked images ("RasterLayer"  or "RasterBrick").
#' @param date Vector with dates for each image.
#' @param adjust Adjustment method passed to \code{\link[stats]{p.adjust}}.
#' Default "none".
#' @description This function computes the theil-sen estimator and 
#' the p-value associated for each pixel over time in a stack of images,
#' writing the values in a raster (one for the estimators and one for 
#' the p-values). It is recommended to use a "RasterBrick", that
#' is more efficient in managing memory.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @seealso \code{\link[rkt]{rkt}}.
#' @examples
#' \dontrun{
#' 
#' require("raster")
#' require("animation")
#' temp <- list()
#'
#'for(i in 1:100) {
#' temp[[i]] <- runif(36,-1, 1)
#' temp[[i]] <- matrix(temp[[i]], 6, 6)
#' temp[[i]] <- raster(temp[[i]])
#'}
#'
#'temp <- brick(temp)
#'
#' oopt <- ani.options(interval = 0.01)
#' for (i in 1:ani.options("nmax")) {
#' plot(temp[[i]])
#' ani.pause()
#' }
#'ani.options(oopt)
#'
#'writeRaster(temp,"temporal.tif", overwrite=T)
#'rm(temp)
#'ndvisim <- brick("temporal.tif")
#'
#'date <- seq(from = 1990.1, length.out = 100, by = 0.2)
#'
#'eco.theilsen(ndvisim, date)
#'
#'pvalue <- raster("pvalue.tif")
#'slope <- raster("slope.tif")
#'par(mfrow = c(1, 2))
#'plot(pvalue, main = "p-value")
#'plot(slope, main = "slope")
#'
#'}
#' @export

setGeneric("eco.theilsen", 
					 function(stacked, date, 
					 				 adjust = c("none", "holm", "hochberg", 
					 										"hommel", "bonferroni", "BH",
					 										"BY", "fdr")) {

	
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
	
	ReturnVal <- tcltk::tkmessageBox(title = "Trend estimation",
																	 message = "process successful!",
																	 icon = "info", type = "ok")
	
	raster::writeRaster(pvalor, "pvalue.tif", overwrite = T)
	
	
})
