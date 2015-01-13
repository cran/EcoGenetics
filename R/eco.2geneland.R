#' Creating input data for Geneland with an ecogen object.
#' @param eco ecogen object.
#' @param ndig Number of digits coding each allele (e.g. 1: x, 2: xx or 
#' 3: xxx) when there is more than one allele per individual locus.
#' @description This function creates four data frames in the workspace  
#' (XY.txt, NAMES.txt, P.txt, G.txt) which can be loaded in Geneland.
#' @return XY.txt matrix with coordinates.
#' @return NAMES.txt matrix with row names.
#' @return P.txt matrix with phenotypic data.
#' @return G.txt matrix with genotypic data.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @examples
#' \dontrun{
#' 
#' data(eco.test)
#' eco.2geneland(eco, 1)
#' 
#' }
#' @export


setGeneric("eco.2geneland", 
					 function(eco, ndig) {
	
 
  write.table(eco@XY, "XY.txt",
  						quote = FALSE,
  						row.names = FALSE,
              col.names = FALSE)
  
  write.table(rownames(eco@XY),
  						"NAMES.txt",
  						quote = FALSE, 
              row.names = FALSE,
  						col.names = FALSE)
  
  write.table(eco@P, 
  						"P.txt", 
  						quote = FALSE, 
  						row.names = FALSE,
              col.names = FALSE)
	
  a <- eco@GENIND$type
  b <- eco@GENIND$ploidy
  if(a == "codom" && b != 1) {
    
  write.table(eco.2columns(eco, ndig),
  						"G.txt",quote = FALSE,
              row.names = FALSE,
  						col.names = FALSE)
  
  } else if(a == "PA" | b == 1) {
    write.table(eco@G, 
    						"G.txt", 
    						quote = FALSE,
                row.names = FALSE, 
    						col.names = FALSE)
  }
  
})

