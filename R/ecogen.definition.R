##########################################
#ECOGEN CLASS: DEFINITION AND METHODS#####
##########################################

#' ecogen class
#' @name ecogen-class
#' @keywords internal
#' @import adegenet
#' @slot XY P data frame
#' @slot P P data frame
#' @slot G G data frame
#' @slot E E data frame
#' @slot S S data frame
#' @slot GENIND genind genetic data
#' @slot C C data frame
#' @slot OUT results
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @aliases ecogen-class


setClass("ecogen",
				 
				 representation(XY = "data.frame",
				 							 P = "data.frame",
				 							 G = "data.frame",
				 							 E = "data.frame",
				 							 S = "data.frame" ,
				 							 GENIND = "genind",
				 							 C = "data.frame",
				 							 OUT = "list"),
				 
				 prototype(XY = data.frame(), 
				 					P = data.frame(),
				 					G = data.frame(),
				 					E = data.frame(),
				 					S = data.frame(), 
				 					GENIND = new("genind"),
				 					C = data.frame(),
				 					OUT = list()
				 ),
				 
)

#----Constructor----------------------------------------------------------##

#' Creating a new ecogen object.
#' @param XY Data frame of 2 columns by n rows (individuals), in UTM format.
#' @param P Data frame of n rows (individuals) with phenotypic data
#' (traits in columns).
#' @param G Data frame of n rows (individuals) with genotypic data 
#' (loci in columns).
#' Data should be in a format accepted by the adegenet 
#' \code{\link[adegenet]{df2genind}} function.  
#' @param E Data frame of n rows (individuals) with environmental data 
#' in columns.
#' @param S Data frame of n rows (individuals) with group factors in columns.
#' Different columns indicate levels or grouping schemas (see details).
#' @param C data frame of n rows (individuals) with custom variables. 
#' @param missing Argument passed to \code{\link[adegenet]{df2genind}} when 
#' G input is a data.frame. Else, when data is of classes "DNAbin" or
#' "alignment", this function can operate over  missing with the same 
#' three values ("0", "NA", or "MEAN"). Missing elements are treated as
#'  zeros in the default option.
#' @param ploidy Ploidy of the G data frame. For an haploid data frame, 
#' ploidy =1 must be included in the function. Default ploidy = 2 for this
#' class of data.
#' @param ... further arguments passed to \code{\link[adegenet]{df2genind}} 
#' when G data is of class "data.frame", otherwise are passed to 
#' \code{\link[adegenet]{DNAbin2genind}} when G data is of class "DNAbin" or to
#' \code{\link[adegenet]{alignment2genind}} when data is of class "alignment".
#' @details This is a generic function for creating an ecogen object. The non 
#' genetic data must be of class "data.frame". For the genetic data, there are
#' three options: Data of "data.frame" class, (for example, microsatellite, 
#' AFLP, RAPD, etc.); data of "DNAbin" class (ape data). The program 
#' recognizes the classes and can discern a codominant data frame from a
#' presence - absence data frame. When data is composed entirely of two
#' values, the program takes these data is of presence - absence
#' type. 
#' Cells in data frames with missing data must be filled with 0.  
#' @seealso \code{\link[adegenet]{df2genind}}
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @examples
#' \dontrun{
#' 
#' #Example with G data of class "data.frame", corresponding to
#' #microsatellites of a diploid organism:
#' data(eco.test)
#' eco <- ecogen(XY = coordinates, P = phenotype, G = genotype,
#' E = environment, S = as.data.frame(structure))
#' eco <- eco.sortalleles(eco, 1)
#' 
#' #Example with G data of class "data.frame", corresponding to a
#' #presence - absence molecular marker:
#' dat <- sample(c(0,1),100,rep = TRUE)
#' dat <- data.frame(matrix(dat,10,10))
#' eco <- ecogen(G = dat)
#' 
#' #Example with G data of class "DNAbin":
#' require(ape)
#' data(woodmouse)
#' G <- woodmouse
#' eco <- ecogen( G = G)
#' 
#' #Example with G data of class "alignment":
#' require(seqinr)
#' data(mase)
#' G <- mase
#' eco <- ecogen( G = G)
#' 
#' }
#' @export ecogen


setGeneric("ecogen",				 
  function(XY = data.frame(),
				 P = data.frame(),
				 G = data.frame(), 
				 E = data.frame(),
				 S = data.frame(),
				 C = data.frame(),
				 missing = c("0", "NA", "MEAN"),
				 ploidy = 2,
				 ...) {				
	
  	
  missing <- match.arg(missing)

	Object <- new("ecogen")
	format <- class(G)
	pre<-c("data.frame",
				 "DNAbin",
				 "alignment")
	
	if(!any(pre %in% format)) {
		stop("data has not a valid format (<data.frame>,
				 <DNAbin> or <alignment>).
				 Check the class of your data.")
	}
	
	
	if(format == "data.frame") {
		

		if(dim(G)[1] != 0) {
			
			if(any(is.na(as.matrix(G))) | length(which(G == "NA")) != 0) {
	stop("NA cells detected in the data frame. Replace NA values with 0.")
			}
			type <- as.factor(as.vector(as.matrix(G)))
			if(length(levels(type)) == 2) {
				type <- "PA"
			} else {
				type <- "codom"
			}
			
			ndat <- nchar(as.matrix(G)[G != 0])
			if((type == "codom") && (any(ndat < 2)) && (ploidy == 2)) {
				stop("\n\n\r"," Some data is coded by 1 character (and are not = 0), 
						 but ploidy = 2. If data is haploid, set ploidy = 1 in the function.","\n\n")
			}
			
			
			obj <- adegenet::df2genind(G, type = type, 
																 missing = missing,
																 ploidy = ploidy, 
																 ...)
			Object@G <- G
			Object@GENIND <- obj
			}
		
	} else if(format == "DNAbin") {
		
		obj <-adegenet::DNAbin2genind(G, ...)
		a0 <- adegenet::genind2df(obj)
		Object@G <- as.data.frame(a0)
		Object@G[is.na(Object@G)] <- 0
		colnames(Object@G) <- colnames(a0)
		rownames(Object@G) <- rownames(a0)
		
		if (missing == "0") {
			obj@tab[is.na(obj@tab)] <- 0
		}
		
		if (missing == "MEAN") {
			moy <- apply(obj@tab, 2,
									 function(c) mean(c, na.rm = TRUE))
			
			for (j in 1:ncol(obj@tab)) {
				temp@tab[, j][is.na(obj@tab[, j])] <- moy[j]
			}
		}
		
		Object@GENIND <- obj
		
	} else if(format == "alignment") {
		
		obj <- adegenet::alignment2genind(G, ...)
		a0 <- adegenet::genind2df(obj)
		Object@G <- as.data.frame(a0)
		Object@G[is.na(Object@G)] <- 0
		colnames(Object@G) <- colnames(a0)
		rownames(Object@G) <- rownames(a0)
		
		if (missing == "0") {
			obj@tab[is.na(obj@tab)] <- 0
		}
		if (missing == "MEAN") {
			moy <- apply(obj@tab, 2,
									 function(c) mean(c, na.rm = TRUE))
			for (j in 1:ncol(obj@tab)) {
				obj@tab[, j][is.na(obj@tab[, j])] <- moy[j]
			}
		}
		
		Object@GENIND <- obj
	}
	
	Object@XY <- XY
	Object@P <- P
	Object@E <- E
	Object@S <- S
	Object@C <- C
	
	attr(Object, "format") <- format
	attr(Object, "type") <-  Object$GENIND$type
	attr(Object, "missing") <- missing
	attr(Object, "ploidy") <- ploidy
	
	return(Object)
	})

##-----Show------------------------------------------------------------##

#' show method
#' @keywords internal
#' @rdname ecogen-methods
#' @aliases show,ecogen-method


setMethod("show", "ecogen", function(object) {
  
  espacio<-function(x, ndig = 10) {
    ndi <- ndig - (3 + nchar(nrow(x)) + nchar(ncol(x)))
    fin <- rep(".", ndi)
    return(fin)
  }
  
  l1 <- paste(nrow(object$XY), "x", ncol(object$XY))
  l2 <- paste(nrow(object$P), "x", ncol(object$P))
  l3 <- paste(nrow(object$G), "x", ncol(object$G))
  l4 <- paste(nrow(object$E), "x", ncol(object$E))
  l5 <- paste(nrow(object$S), "x", ncol(object$S))
  if(ncol(object$S) != 0) {
    l51 <- paste( "-->", ncol(object$S), "structures found") 
  } else { 
    l51 <- ""
  }
  l6 <- paste(ifelse(length(object@GENIND$tab) == 0,"empty", "OK"))
  l7 <- paste(nrow(object$C), "x", ncol(object$C))
  l8 <- paste(length(object$OUT))
 
  out.len<-length(object$OUT)
  
  
  if(out.len != 0) {
    
    out.clas <- character()
    
    for(i in seq(along = object$OUT)) { 
      out.clas[i] <- class(object$OUT[[i]])[1]
    }
    
    out.names <- paste (names(object$OUT)," ", "(", out.clas,")",
                        c(rep(", ",(length(object$OUT) - 1)),""), sep="")
    
  } else {
    out.names <- ""
  }
  
  e <- function(x) rep("", x)
  
  
  cat("\n", e(7), "############################", "\n", e(9), 
      "| ecogen class object |")
  cat( "\n", e(7), "############################")
  cat( "\n\n", "Analysis of phenotypic, genotypic and environmental data",
       "\n")
  cat( "\n", "|", "$XY:", e(6), "|","---->", l1, e(1), e(13 - nchar(l1)),
       "variables")
  cat( "\n", "|", "$P:", e(7), "|", "---->", l2, e(1), e(13 - nchar(l2)),
       "variables")
  cat( "\n", "|", "$G:", e(7), "|", "---->", l3, e(1), e(13 - nchar(l3)),
       "variables")
  cat( "\n", "|", "$E:", e(7), "|", "---->", l4, e(1), e(13 - nchar(l4)), 
       "variables")
  cat("\n", "|", "$S:", e(7), "|", "---->", l5, e(1), e(13 - nchar(l5)),
      "structures", l51)
  cat("\n", "|", "$GENIND:", e(2), "|", "---->", l6)
  cat("\n", "|", "$C:", e(6), "", "|", "---->", l7, e(1),
      e(13 - nchar(l7)), "variables")
  cat("\n","___________________________________________________")
  cat("\n\n", "|", "$OUT:", e(5), "|", "---->", l8, e(1), e(13 - nchar(l8)),
      "analyses")
  if(length(object$OUT) != 0) { 
    cat("\n\n","Results stored:", out.names,"\n\n")
  }
  
})


##-----Summary------------------------------------------------------------##

	#' Summary for ecogen objects
	#' @param object Ecogen object.
	#' @param x The name of the S slot column with the grouping factor 
	#' when a summary taking in account groups is required.
	#' @return Redundancy analysis for the P data frame in function of G 
	#' (constraining data frame)  and E (unconstraining data frame).
	#' @return Mean and sd table for the P data frame.
	#' @return Genetic analysis for the G data frame via hierfstat. This analysis
	#' can be performed when the G data frame has an appropriate hierfstat
	#' input format (i.e., columns of class numeric or integer, neither factor nor character).
	#' @seealso \code{\link{rda}} \code{\link{basic.stats}}
  #' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
	#' @examples
  #' \dontrun{
  #' 
	#' data(eco.test)
	#' summary(eco, "structure")
  #' 
  #' }
  #' @rdname ecogen-summary
  #' @aliases summary,ecogen,ANY-method
  #' @exportMethod summary
 
	

	setMethod("summary", "ecogen", function(object, x = NULL) {
		
		eco <- object
		
		if(is.null(x)){
			eco$S <- data.frame(rep(1, nrow(eco$P)))  
			colnames(eco$S) <- "x"
			x <- "x"
		}
		grupo <- eco$S
		fact<-match(x, colnames(eco$S), nomatch = 0)
		fact <- fact[fact != 0]
		if(length(fact) == 0)
			stop("incorrect factor name")
		
		eco.meansd <- function(eco) {
			
			mi <- function(z, fac, d) {
				x <- tapply(z, fac, mean)
				y <- tapply(z, fac, sd)
				x <- round(x, d)
				y <- round(y, d)
				medias <- rep(0, length(x))
				
				for(i in 1:length(x)) {
					temp <- paste(x[i], " ", "(", y[i], ")", sep = "")
					medias[i] <- temp
				}
				medias <- as.data.frame(medias)
				rownames(medias) <- names(x)
				return(medias)
			}
			
			tab <- apply(eco$P, 2, mi, eco$S[, fact], 3)
			tab <- do.call(cbind, tab)
			colnames(tab) <- colnames(eco$P)
			tab
		}
		cat("\n\n")
		
		cat(paste("====================================",
							"=======================", sep=""),"\n")
		cat(paste("RDA for P~ G and",
							"E (conditioning matrix)"), "\n")
		cat(paste("====================================", 
							"=======================", sep=""), "\n\n")
		test<-vegan::rda(eco$P ~ eco@GENIND$tab, eco$E)
		print(test)
		cat("\n")
		cat(paste("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"), "\n\n")
		print(anova(test))
		cat(paste("-------------------------------------------"), "\n\n")
		
		cat(paste("==================================="), "\n")
		cat(paste("P mean and sd per factor",x),"\n")
		cat(paste("==================================="), "\n\n")
		print(eco.meansd(eco))
		cat(paste("-------------------------------------------"), "\n\n")
		
		cat(paste("========================================="), "\n")
		cat(paste("Basic genetic stats per factor", x), "\n")
		cat(paste("========================================="), "\n\n")
		if((class(eco$G[,1]) != "numeric") & (class(eco$G[,1]) != "integer")) {
			warning("hierfstat does not allow G data frames with 
							factor or character variables")
		print(hierfstat::basic.stats(eco.2hierfstat(eco, x)))
		}
		cat("\n")
		
	})

##-----"$"------------------------------------------------------------##
#' $
#' @rdname ecogen-methods
#' @aliases $,ecogen,character-method
#' @exportMethod $

setMethod("$","ecogen",
					function(x,name) {				
						return(slot(x,name))
					}		
)

##-----"$<-"------------------------------------------------------------##

#' $<-
#' @rdname ecogen-methods
#' @aliases $<-,ecogen,character,ANY-method
#' @exportMethod $<-

setMethod("$<-","ecogen",
					function(x,name,value) {
						
						slot(x,name,check=TRUE) <- value
						return(x)
					}
)

##-----"Names"------------------------------------------------------------##

#' names 
#' @rdname ecogen-methods
#' @aliases names,ecogen-method
#' @exportMethod names

setMethod("names", "ecogen",
					function(x){
						
	return(c("XY", "P", "G",
					 "E","GENIND",
					 "S", "C", "OUT"))
})



##----setAs----------------------------------------------------------------##

# #setAs 
# #@name setAs
# #@keywords internal

#as.ecogen <- ecogen

#setAs ("ecogen" , "genind",
#			 function ( from , to ){
#			 	to <- from@GENIND
#			 	
#			 })

#setAs ("genind" , "ecogen",
#			 function ( from , to ){
#			 	to <- ecogen()
#			 	to@GENIND <-from
#			 	to@G <- genind2df(from)
#			 	to@G[is.na(to@G)] <- 0
#			 	to
#			 })


#' Conversion between genind and ecogen, and vice versa
#' @name genind2ecogen
#' @param gen genind object
#' @rdname ecogen2genind-genind2ecogen
#' @exportMethod genind2ecogen

setGeneric("genind2ecogen", function(gen) {

	z <- ecogen()
	z@GENIND <- gen
	return(z)
	
})

#' 
#' @name ecogen2genind
#' @param eco ecogen object
#' @rdname ecogen2genind-genind2ecogen
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @exportMethod ecogen2genind

setGeneric("ecogen2genind", function(eco) {
	
	return(eco@GENIND)
	
})

##-----is -----------------------------------------------------------------##
#' is.ecogen
#' @rdname ecogen-methods
#' @aliases is,ecogen-method
#' @export is.ecogen

is.ecogen <- function(x) {
	value <- is(x, "ecogen")
  value
}

##----"["-----------------------------------------------------------------##
#' [
#' @rdname ecogen-methods
#' @aliases [,ecogen,integer,missing-method
#' @exportMethod [

	setMethod("[", c("ecogen", "integer", "missing", "ANY"),
						
						function(x, i, j,..., drop=FALSE) {
							
							initialize(x, 
												 XY = x@XY[i, , drop = FALSE],
												 P = x@P[i, , drop = FALSE],
												 G = x@G[i, , drop = FALSE], 
												 GENIND = df2genind(x@G[i, , drop = FALSE], 
												 									 missing = attr(x, "missing"),
												 									 ploidy = attr(x, "ploidy"),
												 									 ...),
												 E = x@E[i, , drop = FALSE], 
												 S = x@S[i, , drop = FALSE],
												 C = x@C[i, , drop = FALSE], 
												 OUT = list())
						})
						

##----------------------------------------------------------------------------##
