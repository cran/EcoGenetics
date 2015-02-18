##########################################
#ECOGEN CLASS: DEFINITION AND METHODS#####
##########################################
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015

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

##-----"$"----------------------------------------------------------------##


setMethod("$","ecogen",
     function(x,name) {    
      return(slot(x,name))
     }  
)

##-----"$<-"--------------------------------------------------------------##


setMethod("$<-","ecogen",
     function(x,name,value) {
      
      slot(x,name,check=TRUE) <- value
      return(x)
     }
)

##-----"Names"------------------------------------------------------------##


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
#    function ( from , to ){
#     to <- from@GENIND
#     
#    })

#setAs ("genind" , "ecogen",
#    function ( from , to ){
#     to <- ecogen()
#     to@GENIND <-from
#     to@G <- genind2df(from)
#     to@G[is.na(to@G)] <- 0
#     to
#    })


#' Conversion between genind and ecogen, and vice versa


setGeneric("genind2ecogen", function(gen) {

 z <- ecogen()
 z@GENIND <- gen
 return(z)
 
})


setGeneric("ecogen2genind", function(eco) {
 
 return(eco@GENIND)
 
})

##-----is -----------------------------------------------------------------##


is.ecogen <- function(x) {
 value <- is(x, "ecogen")
  value
}

##----"["-----------------------------------------------------------------##


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
