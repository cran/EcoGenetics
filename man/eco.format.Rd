% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eco.format.R
\name{eco.format}
\alias{eco.format}
\title{Format tool for genetic data}
\usage{
eco.format(
  data,
  ncod = NULL,
  nout = 3,
  ploidy = 2,
  sep.in = "",
  sep.out = "",
  fill.mode = c("last", "first", "none"),
  recode = c("none", "all", "column", "paired"),
  replace_in = NULL,
  replace_out = NULL,
  show.codes = FALSE
)
}
\arguments{
\item{data}{Genetic data frame.}

\item{ncod}{Number of digits coding each allele in the input file.}

\item{nout}{Number of digits in the output.}

\item{ploidy}{Ploidy of the data.}

\item{sep.in}{Character separating alleles in the input data if present.}

\item{sep.out}{Character separating alleles in the output data. Default}

\item{fill.mode}{Add zeros at the beggining ("fist") or the end ("last")
of each allele. Default = "last".}

\item{recode}{Recode mode: "none" for no recoding (defalut), "all" for recoding
the data considering all the individuals values at once (e.g., protein data), 
"column" for recoding the values by column (e.g., microsatellite data), "paired" 
for passing the values of allelic states and corresponding replacement values, using 
the replace_in and replace_out arguments (e.g. replace_in = c("A", "T", "C", "G"),
replace_out = c(1,2,3,4)).}

\item{replace_in}{vector with states of the data matrix to be replaced, when recode = "paired".
This argument must be used in conjunction with the argument "replace_out".}

\item{replace_out}{vector with states of the data matrix used for replacement, when recode = "paired".
This argument must be used in conjunction with the argument "replace_in".}

\item{show.codes}{May we returned tables with the equivalence between the old 
and new codes when recode = "all" or recode = "column"?}
}
\description{
Format tool for genetic data
}
\details{
The function can format data with different ploidy levels. 
It allows to:
- add/remove zeros at the beginning/end of each allele
- separate alleles with a character
- divide alleles into columns - bind alleles from separate columns
- transform character data into numeric data

"NA" is considered special character (not available data).

Note that the function can also work with other type of data as well, 
where the "alleles" represent the different states of the variables.
}
\examples{
\dontrun{

data(eco.test)

# Adding zeros

example <- as.matrix(genotype[1:10,])
mode(example) <- "character"
# example data
example                
recoded <- eco.format(example, ncod = 1, ploidy = 2, nout = 3)
# recoded data
recoded


# Tetraploid data, separating alleles with a "/"
tetrap <- as.matrix(example)
# simulated tetraploid example data
tetrap <- matrix(paste(example,example, sep = ""), ncol = ncol(example)) 
recoded <- eco.format(tetrap, ncod = 1, ploidy = 4, sep.out = "/")
# recoded data
recoded

# Example with a single character
ex <- c("A","T","G","C")
ex <- sample(ex, 100, rep= T)
ex <- matrix(ex, 10, 10)
colnames(ex) <- letters[1:10]
rownames(ex) <- LETTERS[1:10]
# example data
ex  
recoded <- eco.format(ex, ploidy = 1, nout = 1,  recode = "all", show.codes = TRUE) 
# recoded data 
recoded

# Example using values-replacement pairs
samp1 <- sample(c("A","T","G","C"), 100, replace = TRUE)
samp2 <- sample(c("A","T","G","C"), 100, replace = TRUE)
paired <- matrix(paste0(samp1, samp2), 10, 10)
# Generate some NAs
paired[sample(1:100, 10)]<-NA
out <- eco.format(paired, recode = "paired", replace_in = c("A", "T", "G", "C"),
replace_out = c(1, 2, 3, 4))
out


# Example with two strings per cell and missing values:
ex <- c("Ala", "Asx", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", 
"Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", 
"Val", "Trp")
ex1 <- sample(ex, 100, rep= T)
ex2 <- paste(ex1, ex1, sep="")
missing.ex2 <- sample(1:100, 20)
ex2[missing.ex2] <-NA
ex2 <- matrix(ex2, 10, 10)
colnames(ex2) <- letters[1:10]
rownames(ex2) <- LETTERS[1:10]
# example data
ex2
recoded <- eco.format(ex2, ncod = 3, ploidy = 2, 
                      nout = 2, recode = "column")
# recoded data
recoded

# Example with a vector, following the latter example:
ex1 <- as.data.frame(ex1)
# example data
ex1
recoded <- eco.format(ex1, ploidy = 1, recode = "all")
# recoded data
recoded

}

}
\author{
Leandro Roser \email{learoser@gmail.com}
}
