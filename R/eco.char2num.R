#' Converting data coded with characters into data coded by numbers.
#' 
#' @param data Data frame or matrix for the conversion.
#' @param ndig Number of digits coding each allele when each cell has two 
#' alleles. Default takes each cell as having a character or string that 
#' indicates a single information (see examples). 
#' @details This function creates a data frame with data recoded in numeric
#' format.
#' The output possess a character formatting, for maintaining an uniform
#' number of digits (for example, data with the code "01" in a 
#' numeric format will be coded as "1", but as many programs requires inputs
#' with uniform number of digits, the program maintains a character format).
#' The program also creates a data frame with the information
#' of the numeric code that corresponds to each character/string. 
#' The program expects missing data coded as NA. 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @examples
#' \dontrun{
#' 
#' # Example with a single character:
#' ex <- c("A","T","G","C")
#' ex <- sample(ex, 100, rep= T)
#' ex <- matrix(ex, 10, 10)
#' colnames(ex) <- letters[1:10]
#' rownames(ex) <- LETTERS[1:10]
#' recoded <- eco.char2num(ex)
#' tab <- recoded$recoded_data
#' write.table(tab, "recoded.data.txt")
#' 
#' # Example with two strings per cell and missing values:
#' ex <- c("Ala", "Asx", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", 
#' "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", 
#' "Val", "Trp")
#' ex1 <- sample(ex, 100, rep= T)
#' ex2 <- sample(ex, 100, rep= T)
#' ex3 <- paste(ex1, ex2, sep="")
#' missing.ex3 <- sample(1:100, 20)
#' ex3[missing.ex3] <-NA
#' ex4 <- matrix(ex3, 10, 10)
#' colnames(ex4) <- letters[1:10]
#' rownames(ex4) <- LETTERS[1:10]
#' ex4
#' recoded <- eco.char2num(ex4, ndig = 3)
#' 
#' # Example with a vector, following the last example:
#' ex1 <- as.data.frame(ex1)
#' ex1
#' eco.char2num(ex1)
#' 
#' }
#' 
#' @export

setGeneric("eco.char2num", 
					 
					 function(data, ndig = 0) {
  
						
  singlechar <- function(x) {
    y <- as.vector(as.matrix(x))
    y <- as.factor(y)
    original.code <- levels(y)
    y <- as.numeric(y)
    max.char <- max(nchar(y))
    y <- formatC(y, width=max.char, flag="0")
    y <- as.factor(y)
    new.code <- levels(y)[levels(y) != "NA"]
    y <- as.character(y)
    y[y == "NA"] <- paste(rep(0,max.char), collapse = "")
    y.tab <- data.frame(original.code, new.code)
    res <- list(y, y.tab)
    res
 
  }
  
  if(ndig == 0) {
    res <- singlechar(data)
    res[[1]] <- data.frame(matrix(res[[1]], ncol  = ncol(data),
                                  nrow = nrow(data)))
    for(i in 1:ncol(res[[1]])) {
      res[[1]][,i] <- as.character(res[[1]][,i])
    }
    colnames(res[[1]]) <- colnames(data)
    rownames(res[[1]]) <- rownames(data)
    names(res) <- c("recoded_data", "conversion_table")
    res
      
  } else if(ndig != 0) {
    
    geno <- as.matrix(data)
    
    m1 <- substr(geno, 1, ndig)
    m2 <- substr(geno, ndig + 1, 2 * ndig)
    m2 <-c(m1, m2)
    m2[m2 == ""] <- NA
    recoding <- singlechar(m2)
    cadena <- recoding[[1]]
    m1 <- cadena[1:length(m1)]
    m2 <- cadena[(length(m1)+1): length(cadena)]
    m3 <- paste(m1, m2, sep ="")
    m3 <- matrix(m3, ncol  = ncol(data), nrow = nrow(data))
    m3 <-data.frame(m3)
    colnames(m3) <- colnames(data)
    rownames(m3) <- rownames(data)
    res <- list("recoded_data" = m3, "conversion_table" = recoding[[2]])
    for(i in 1:ncol(res[[1]])) {
      res[[1]][,i]<-as.character(res[[1]][,i])
    }
    res
  }
    
  
})
	
