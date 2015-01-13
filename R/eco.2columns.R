#' Exporting an ecogen genetic data frame with two locus per column
#' into a data frame with one locus per column format.
#' @param eco ecogen object.
#' @param ndig Number of digits coding each allele
#'  (e.g., 1: x, 2: xx, or 3: xxx).
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @include ecogen.definition.R
#' @examples
#' \dontrun{
#' 
#' data(eco.test)
#' eco.2al <- eco.2columns(eco, 1)
#' eco.2al
#' write.table(eco.2al, "one.per.column.txt", quote = FALSE)
#' 
#' }
#' @export

setGeneric("eco.2columns", 
					 function(eco, ndig) {
						
						
						geno <- as.matrix(eco@G)
						
						m1 <- substr(geno, 1, ndig)
						
						m1 <- matrix(m1, ncol = ncol(geno),
												 nrow = nrow(geno))
						
						m2 <- substr(geno, ndig + 1, 2 * ndig)
						
						m2 <- matrix(m2, ncol = ncol(geno),
												 nrow = nrow(geno))
						
						m3 <- matrix(, ncol = 2 * ncol(geno),
												 nrow = nrow(geno))
						
						colnames(m3) <- c(1:ncol(m3))
						
						j <- 1
						for(i in 1:ncol(m3)) {
							
							if((i %% 2) == 0) {
								
								m3[, i] <- m2[, j]
								j <- i/2
								colnames(m3)[i] <- paste(colnames(geno)[j], ".2", sep = "")
							} else {
								if(i != 1) {
									j <- j + 1
								}
								m3[, i] <- m1[, j]
								colnames(m3)[i] <- paste(colnames(geno)[j], ".1", sep = "")
							}
						}
						
						m3[m3 == "0"] <- "NA"
						m3[m3 == ""] <- "NA"
						m3 <- as.data.frame(m3)
						
						m3
						
					})
