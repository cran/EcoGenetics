# Exporting an ecogen genetic data frame into SPAGeDI format

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015

eco.2spagedi <- function(eco, 
                         pop = NULL, 
                         ndig, 
                         name = "infile.spagedi.txt", 
                         smin = 0,
                         smax= NULL,
                         int = NULL, 
                         nclass = NULL,
                         seqvec = NULL,
                         size = NULL,
                         bin = c("sturges", "FD"),
                         distmat = NULL,
                         latlon = FALSE) {
  
  bin <- match.arg(bin)
  if(is.null(pop)) {
    pop <- rep(1, nrow(eco$XY))
  } else {
    pop <- as.numeric(eco$S[, which(colnames(eco$S)== pop)])
  }
  
  if(sum(pop == 0)) {
    stop("non matching S column name")
  }
  
  matriz <- data.frame(rownames(eco$P), pop, eco$XY, eco$G)
  matriz <- as.matrix(matriz)
  colnames(matriz) <- c("Individual", "Population",
                        colnames(eco$XY), colnames(eco$G))
  
  arriba <- rep("", ncol(matriz))
  arriba[1] <- nrow(matriz)
  arriba[2] <- max(pop)
  arriba[3] <- ncol(eco$XY)
  arriba[4] <- ncol(eco$G)
  arriba[5] <- ndig
  arriba[6] <- eco$GENIND$ploidy
  
  if(!is.null(smax) | !is.null(seqvec)) {
    
    xy <- eco$XY[,1:2]
    if(latlon) {
      dist(SoDA::geoXY(xy[,2], xy[,1], unit=1))
    }
    
    input <- int.break(XY = xy, 
                       int = int, 
                       smin =smin,
                       smax = smax,
                       nclass = nclass,
                       seqvec = seqvec, 
                       size = size,
                       bin = bin)
    
    breakpoints <- input$breakpoints
    
  } else  {
    breakpoints <- NULL
  }
  
  final <- rep("", ncol(matriz))
  final[1] <- "END"
  
  
  sink(name)
  cat(arriba, sep = "\t")
  cat("\n")
  if(!is.null(breakpoints)) {
    cat(length(breakpoints), breakpoints, sep = "\t")
    cat("\n")
  } else {
    cat(0)
  }
  cat("\n")
  cat(colnames(matriz), sep = "\t")
  cat("\n")
  write.table(matriz, sep = "\t", quote = FALSE, row.names = FALSE, 
              col.names = FALSE)
  cat(final, sep = "\t")
  
  if(!is.null(distmat)) {
    if(class(distmat) == "dist") {
      distmat <- as.matrix(distmat)
    } 
    if(class(distmat) != "matrix" & class(distmat) != "data.frame") {
      stop("invalid distance matrix format (It should be of class: dist, matrix or data.frame")
    }
    distnames <- rownames(distmat)
    distmat<-data.frame(distnames, distmat)
    colnames(distmat)[1]<-paste("M", nrow(eco$XY), sep="")
    cat("\n")
    write.table(distmat, sep = "\t", quote = FALSE, row.names = FALSE, 
                col.names = TRUE)
    cat("END")
  }
  
  sink()
  
}

