#'  Converting an ecogen genetic data frame into a gstudio object.
#' @param eco ecogen object
#' @param ... Further arguments passed to \code{\link[adegenet]{df2genind}}
#' @param type The type of data passed to \code{\link[gstudio]{locus}}. 
#' Default is "separated" (data as microsatellites or individual haplotypes);
#' "aflp" for presence - absence data. 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @examples
#' \dontrun{
#' 
#' data(eco.test)
#' gsteco <- eco.2gstudio(eco, "separated")
#' gsteco
#' 
#' }
#' @export

setGeneric("eco.2gstudio", 
					 function(eco, type = "separated", ...) {
						
						
						if(type == "separated") {
							dat <- adegenet::df2genind(eco$G)
							dat <- adegenet::genind2df(dat, sep = ":")
							for(i in 1:ncol(dat)) {	
								dat[, i] = gstudio::locus(dat[, i], type = "separated")
							}
						} else {
							dat<-eco$G
							for(i in 1:ncol(dat)) {
								dat[, i] = gstudio::locus(dat[, i], type = type)
							}
						}
						
						dat
					})
