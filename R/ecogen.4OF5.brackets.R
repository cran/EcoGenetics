
################################################
#### $, $<-, "[", "[<-", "[[, AND "[[<-" METHODS
################################################


## "$"-------------------------------------------------------------------------#
#' $
#' @rdname ecogen-methods
#' @aliases $,ecogen,character-method 
#' @exportMethod $
#' 

# DEPRECTED

setMethod("$","ecogen",
          function(x,name) {	
mess <- message(
paste("NOTE: This method has been deprecated in EcoGenetics 1.2.0-2.\n",
"     Use instead the accessor", aue.access(name, deparse(substitute(x))), "or double square brackets,\n",
"    ", paste(deparse(substitute(x)), "[[", deparse(substitute(name)),"]].", sep = ""), 
"See help('EcoGenetics accessors')."))
return(mess)
          })


## "$<-"-----------------------------------------------------------------------#
#' $<-
#' @rdname ecogen-methods 
#' @aliases $<-,ecogen,character,ANY-method 
#' @exportMethod $<-

# DEPRECTED

setMethod("$<-","ecogen",
          function(x,name,value) {
            tmp <- "name_of_this_object"
            mess <- message(
              paste("NOTE: This method has been deprecated in EcoGenetics 1.2.0-2.\n",
                    "     Use instead the accessor, ", 
                    paste(aue.access(name, tmp), "<-,\n", sep = ""),
                    "     or double square brackets,",
                    "", paste(tmp, "[[", deparse(substitute(name)),"]]<-.", sep = ""), 
                    "\n      See help('EcoGenetics accessors')."))
            print(mess)
            return(x)
          })


## "["-------------------------------------------------------------------------#
#' [ 
#' @rdname ecogen-methods 
#' @aliases [,ecogen,integer,missing-method 
#' @exportMethod [

setMethod("[", c("ecogen", "numeric", "missing", "ANY"), 
          
          function(x, i, j, ..., drop = FALSE) {
            
            
            # create an int.genind object if nrow(G) != 0
            if(all(dim(x@G) != 0)) {
              tempo <- int.df2genind(x@G[i, , drop = FALSE], 
                                     missing = x@INT@missing,
                                     ploidy = x@INT@ploidy,
                                     type =  x@INT@type,
                                     ...)
            } else {
              tempo <- new("int.genind")
            }
            
            z <- new("ecogen")
            
            # if(all...) condition required because subsetting over matrices of 
            # dimension 0 generates a matrix of dimension i x 0 (undesired)
            if(all(dim(x@XY) != 0)) {
            z@XY <- x@XY[i, , drop =FALSE]
            }
            if(all(dim(x@P) != 0)) {
            z@P <- x@P[i, , drop =FALSE]
            }
            
            if(all(dim(x@G) != 0)) {
            z@G <- as.data.frame(int.genind2df(tempo), 
                                 stringsAsFactors = FALSE)
            }
            
            # IF TYPE == DOMINANT -> A is empty
            if(x@INT@type == "codominant") {
            if(all(dim(x@A) != 0)) {
            z@A <- data.frame(tempo@tab)
            }
            }
            
            if(all(dim(x@E) != 0)) {
            z@E <- x@E[i, , drop =FALSE]
            }
            
            #all S columns as factors
            if(all(dim(x@S) != 0)) {
              Sout <- x@S[i, , drop = FALSE]
              for(n in 1:(ncol(Sout))) {
                Sout[, n] <- factor(Sout[, n])
              }
            } else {
              Sout <- data.frame()
            }
            z@S <- Sout
            
            if(all(dim(x@C) != 0)) {
            z@C  <- x@C[i, , drop =FALSE]
            }
            z@OUT  <- list()
            z@INT <- int.genind2gendata(tempo)
            
           z
          
          })

## [<--------------------------------------------------------------------------#
#[<-
# @rdname ecogen-methods
# @aliases [,ecogen-method
# @exportMethod [<-

#setReplaceMethod("[","ecogen",  function (x,i,j,value) {
#message(paste("undefined operation for ecogen objects (would you mean",
#               paste("[[", x, "]]", sep = ""), "?"))
#})

## [[ internal-----------------------------------------------------------------#
#' [[
#' @keywords internal

# THIS IS AN INTERNAL METHOD. IT INCLUDES THE SLOT INT.

setMethod("[[", c("ecogen","numeric", "missing"), function(x, i, j) {
  if(i == 1) return(x@XY)
  if(i == 2) return(x@P)
  if(i == 3) return(x@G)
  
  # if x is dominant, A return G
  if(i == 4) {
    if(x@INT@type == "codominant") {
      return(x@A)
    } else {
      return(NULL)
    }
  }
  
  if(i == 5) return(x@E)
  if(i == 6) return(x@S)
  if(i == 7) return(x@C)
  if(i == 8) return(x@OUT)
  if(i == 9) return(x@INT) #INVISIBLE FOR PUBLIC CHARACTER METHOD 
  if(!any(c(1:9) %in% i)) {
    return(NULL)
  }
})

## [[--------------------------------------------------------------------------#
#' [[
#' @rdname ecogen-methods 
#' @aliases [[,ecogen,character,missing-method
#' @exportMethod [


setMethod("[[", c("ecogen","character", "missing"), function(x, i, j) {
  if(toupper(i) ==  "XY") return(ecoslot.XY(x))
  if(toupper(i) == "P") return(ecoslot.P(x))
  if(toupper(i) == "G") return(ecoslot.G(x))
  if(toupper(i) == "A") return(ecoslot.A(x)) # DOMINANCE / CODOMINANCE HANDLED WITH THE ACCESSOR
  if(toupper(i) == "E") return(ecoslot.E(x))
  if(toupper(i) == "S") return(ecoslot.S(x))
  if(toupper(i) == "C") return(ecoslot.C(x))
  if(toupper(i) == "OUT") {
    if(length(x@OUT) != 0) {
      return(x@OUT)
    } else {
      return("OUT is empty")
    }
  }
  if(!toupper(i) %in% c("XY", "P", "G", "A","E", "S", "C", "OUT")) {
    message(paste(paste("<", i, ">", sep = ""), "is an undefined ecogen slot"))
  }
})


## [[<- internal---------------------------------------------------------------#
#' [[<-
#' @keywords internal

# THIS IS AN INTERNAL METHOD. IT INCLUDES THE SLOT INT.

setReplaceMethod("[[", c("ecogen","numeric", "missing"), function (x, i, j, ..., value) {
  if(i == 1) ecoslot.XY(x) <- value
  if(i == 2) ecoslot.P(x) <- value
  if(i == 3) ecoslot.G(x, ...) <- value  # if x is dominant, A return G
   
  if(i == 4) { # CODOMINANT / DOMINANT DEPENDENT
      if(x@INT@type == "codominant") {
        ecoslot.A(x) <- value
      } else {
        return(x)
      }
    } 
  
  if(i == 5) ecoslot.E(x) <- value
  if(i == 6) ecoslot.S(x) <- value
  if(i == 7) ecoslot.C(x) <- value
  if(i == 8) ecoslot.OUT(x) <- value
  if(i == 9) int.ecoslot.INT(x) <- value
  if(!any(c(1:9) %in% i)) {
    message("invalid slot")
  }
  return (x)
})


#' [[<- -----------------------------------------------------------------------#
#' @rdname ecogen-methods 
#' @aliases [,ecogen,character,missing,ANY-method
#' @exportMethod [[<-

setReplaceMethod("[[",c("ecogen","character", "missing"),  function (x, i, j,..., value) {
  
  if(toupper(i) == "XY") ecoslot.XY(x) <- value
  if(toupper(i) == "P") ecoslot.P(x) <- value
  if(toupper(i) == "G") ecoslot.G(x, ...) <- value
  if(toupper(i) == "A") ecoslot.A(x) <- value # DOMINANCE / CODOMINANCE HANDLED WITH THE ACCESSOR
  if(toupper(i) == "E") ecoslot.E(x) <- value
  if(toupper(i) == "S") ecoslot.S(x) <- value
  if(toupper(i) == "C") ecoslot.C(x) <- value
  
  # SLOT OUT - deparse-substitute arguments are complex to 
  # handle with ecoslot.OUT(x) <- value. The code below is a copy
  # of ecoslot.OUT
  #---------------------------------------------------------------------
  if(toupper(i) == "OUT") {
   
    object <- x
    
    #empty list -> new slot OUT empty
    if(is.list(value) && length(value) == 0) {
      object@OUT <- list()
      return(object)
    } # return object
    
    # obtention of names of the argument <value>
    ## split arguments if several present
    if(is.list(value)) {
      res.names <- substitute(value)
      res.names <- lapply(res.names, deparse)[-1]
      res.names <- unlist(res.names)
    } else {
      ## one argument
      res.names <- deparse(substitute(value))
    }
    
    # if vector -> conversion into list
    if(!is.list(value)) {
      value <- list(value)
    }
    
    # remotion of quotation marks and spaces
    res.names <- gsub("[\"]| ", "", res.names)
    
    #abbreviate names if required 
    # if(!is.null(abbr)) {
    # res.names <- as.vector(abbreviate(res.names, abbr))
    # }
    
    #-------------------------------#
    Z <- object
    
    # fill OUT slot
    
    # original names 
    tmp <- names(Z@OUT)
    # add elements to Z
    Z@OUT <- c(Z@OUT, value)
    # add names to Z. Names must be unique
    names(Z@OUT)<- make.unique(c(tmp, eval(res.names)))
    # order names
    orden <- order(names(Z@OUT))
    Z@OUT <- Z@OUT[orden]
    
    return(Z)
    
  }
  #---------------------------------------------------------------------
  
  if(!toupper(i) %in% c("XY", "P", "G", "A","E", "S", "C", "OUT")) {
    message(paste(paste("<",i,">", sep = ""), "is an undefined ecogen slot"))
  }
  return(x)
})
