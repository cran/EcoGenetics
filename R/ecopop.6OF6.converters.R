
#' Conversion form ecogen to ecopop 
#' 
#' @description This function creates an ecopop object from an ecogen object
#' @param from Object of class "ecogen"
#' @param to Object of class "ecopop"
#' @rdname ecogen2ecopop
#' @examples
#' 
#' \dontrun{
#'
#' data(eco.test)
#' ecogen2ecopop(eco, hier = "pop")
#'
#' }
#' 
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @export

setGeneric("ecogen2ecopop", function(from, hier, factor_to_dummy = TRUE,
                                     aggregator = function(x) mean(x, na.rm = TRUE)) { 
  
which_pop <- which(colnames(from@S)== hier)
if(length(which_pop) == 0) {
  stop("non matching S column name")
}

pop <- from@S[, hier]

to <- new("ecopop")

to@XY <-  aue.aggregated_df(from@XY, pop, aggregator, factor_to_dummy = FALSE)

to@P <-   aue.aggregated_df(from@P, pop, aggregator, factor_to_dummy = factor_to_dummy )

if(from@INT@type == "codominant") {
  to@AF <- as.matrix(apply(from@A, 2, tapply, structure, sum, na.rm=TRUE))
} else {
  to@AF <-  as.matrix(apply(from@G, 2, tapply, structure, sum, na.rm=TRUE))
}

to@E <-  aue.aggregated_df(from@P, pop, aggregator, factor_to_dummy = factor_to_dummy )

to@S <-  factor(levels(pop))

to@C <-   aue.aggregated_df(from@C, pop, aggregator, factor_to_dummy = factor_to_dummy )

popdat <- new("int.popdata")
popdat@ploidy <- from@INT@ploidy
popdat@type <- from@INT@type
popdat@aggregator <- aggregator
popdat@factor_to_dummy <- factor_to_dummy
popdat@loc.fac <- from@INT@loc.fac
popdat@all.names <- from@INT@all.names
to@INT <- popdat

to
})


#' Conversion form ecopop to genpop and genpop to ecopop
#' 
#' @description These functions export from ecopop to genpop and viceversa
#' @param from Object of class "ecopop" / "genpop"
#' @rdname ecopop2genpop
#' @examples
#' 
#' \dontrun{
#' data(eco.test)
#' my_ecopop <- ecogen2ecopop(eco, hier = "pop")
#' ecpop2genpop(my_ecopop)
#' }
#' 
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @export


setGeneric("ecopop2genpop", function(from) { 
  
  if(!require(adegenet)) stop("Please install the adegenet package first")
  
  to <- adegenet::genpop()
  
  if(!any(dim(from@XY) == 0)) {
    to@other$xy <- from@XY
  }
  
  if(!any(dim(from@AF) == 0)) {
    to@tab <- from@AF
    to@loc.fac <- from@INT@loc.fac
    nomloc <- levels(from@INT@loc.fac)
    temp<- tapply(rep(1, length(from@INT@loc.fac)), from@INT@loc.fac, sum)
    #array to list ("as.list()" do not works with the array)
    temp <- as(temp, "numeric")
    names(temp) <- nomloc
    to@loc.n.all <- temp
    temp <-  tapply(from@INT@all.names, names(from@INT@all.names), function(x) return(unname(x)),simplify = FALSE)
    #reorder temp an convert to list
    temp <- temp[pmatch(nomloc, names(temp))]
    nomloc <- names(temp)
    temp <- as(temp, "list")
    names(temp) <- nomloc
    to@all.names <- temp
    to@ploidy <- rep(from@INT@ploidy, nrow(from@AF))
    to@type <- ifelse(from@INT@type == "codominant", "codom", "PA")
  }

  to
})

#' genind2ecogen
#' @rdname ecogen2genind
#' @export


setGeneric("genpop2ecopop", function(from) { 
  
  if(!require(adegenet)) stop("Please install the adegenet package first")
  
  to <- ecopop(AF = from@tab, S = as.factor(rownames(from@tab)))
  to@INT@loc.fac <- from@loc.fac
  
  counts <- lapply(from@all.names, length)
  xnames <- names(xx@all.names)
  xnames <- rep(xnames, counts)
  to@INT@all.names <- unlist(from@all.names)
  names(to@INT@all.names) <- xnames
  
  to@INT@ploidy <- unique(from@ploidy)
  if(length(to@INT@ploidy) > 1) {
    stop("multiple ploidy levels are not supported by ecopop objects")
  }
  
  to@INT@type <- ifelse(from@type == "codom", "codominant", "dominant")
  
  if(!is.null(from@other$xy)) {
    to@XY <- from@other$xy
  }
  
  to@ATTR$whereIs <- parent.frame()
  to@ATTR$.call <- match.call()

  to
})