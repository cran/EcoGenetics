#' Sorting alleles of an ecogen genetic data frame.
#' @param eco ecogen object.
#' @param ndig Number of digits coding each allele (e.g. 2: xx, or 3: xxx).
#' @description This program returns a data fame with the alleles
#' of each individual \emph{i} and each loci \emph{j} in ascending order.
#' For example, a locus of type 51 is returned as 15.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @examples
#' \dontrun{
#' 
#' data(eco.test)
#' eco <- ecogen(XY = coordinates, P = phenotype, G = genotype, 
#' E = environment, S = as.data.frame(structure), missing = 0)
#' eco$G
#' eco <- eco.sortalleles(eco, 1)
#' eco$G
#' 
#' }
#' @export

setGeneric("eco.sortalleles",
					 function(eco, ndig)  {

geno <- as.matrix(eco$G)

m1 <- substr(geno, 1, ndig)
m2 <- substr(geno, ndig + 1, 2 * ndig)
m3 <- matrix(, ncol = ncol(geno), nrow = nrow(geno))
m4 <- matrix(, ncol = ncol(geno), nrow = nrow(geno))

m3[which(m1 < m2)] <- m1[which(m1 < m2)]
m3[which(m2 < m1)] <- m2[which(m2 < m1)]
m3[which(m2 == m1)] <- m2[which(m2 == m1)]

m4[which(m1 > m2)] <- m1[which(m1 > m2)]
m4[which(m2 > m1)] <- m2[which(m2 > m1)]
m4[which(m2 == m1)] <- m2[which(m2 == m1)]

m5 <- paste(m3, m4, sep = "")
m5 <- as.data.frame(matrix(m5, ncol = ncol(geno), nrow = nrow(geno)))
m6 <- eco
m6$G <- m5
m6

})
