
#' Conversion of a non symmetric binary matrix into symmetric.
#' @param mat Matrix.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal


misc.2symmetric <- function(mat) {
  mat[row(mat) > col(mat)] <- mat[row(mat) > col(mat)] + mat[row(mat) < col(mat)]
  mat[row(mat) > col(mat)][mat[row(mat) > col(mat)] > 0] <- 1
  mat[row(mat) < col(mat)] <- mat[row(mat) > col(mat)]
  mat
}


#' Creates a matrix without diagonal, in row order
#' @param mat Matrix.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

misc.undimmattg <- function(mat) {
  ncolp <- ncol(mat) -1
  mat2 <- as.vector(t(mat))
  mat2 <-mat2[-which(col(mat) == row(mat))]
  mat2<-matrix(mat2, ncol = ncolp, byrow = T)
  mat2
}


#' Computing a distance matrix in meters among points in decimal degrees
#' under a spherical Earth model
#' @param XY data frame or matrix with latitude-longitude coordinates 
#' in decimal degrees format.
#' This program computes a distance matrix for Earth points in decimal degrees.
#' It assumes a spherical model with an Earth radius of 6371 km. 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal


misc.dlatlon2distm <- function(XY) {
  out <- matrix(,nrow(XY), nrow(XY)) 
  for(i in 1:nrow(XY)) {
    for(j in 1:nrow(XY)) {
      lat1 <- XY[i, 1] 
      lon1 <- XY[i, 2] 
      lat2 <- XY[j, 1] 
      lon2 <- XY[j, 2] 
      R <- 6371                                
      dLat <- (lat2 - lat1) * pi / 180
      dLon <- (lon2 - lon1) * pi / 180
      a <- sin((dLat/2)) ^ 2 + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * (sin(dLon / 2)) ^2
      c <- 2 * atan2(sqrt(a), sqrt(1-a))
      d <- R * c  
      out[i, j] <- d
    }
  }
  rownames(out) <- rownames(XY)
  colnames(out) <- rownames(XY)
  as.dist(out)
}


#' Filter a raster using a conditional expression and values in a conditional vector
#' @param data raster data
#' @param expr logical expression to apply a a character string. Must have the
#' argument "filter" and the form "filter < 3", "filter == 3" || filter< 4"
#' @param filter vector with the values to use for filtering the data (included
#' or excluded of the matrix). If a column or row of the data, filter is 
#' the value data[row, ], data[, col].
#' @param byrow filter the rows? default FALSE (filter the columns)
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

misc.parse.filter <- function(data, expr, filter, byrow = TRUE) {
  
  
  #conditional vector checkpoint
  l.cond <- length(filter)
  if(byrow) {
    if(l.cond != nrow(data)) {
      stop("the conditional vector do not match the number of rows of the data")
    } 
  } else {
    if(l.cond != ncol(data)) {
      stop("the conditional vector do not match the number of columns of the data")
    }
  }
  
  # conditional expression checkpoint
  ## obtain expression. Logical o character allowed
  cond <- deparse(substitute(expr))
  cond <- gsub("\\\"", "", cond)
  cond <- gsub(" ", "", cond)
  
  # evaluate expression
  cond <- eval(parse(text = cond))
  
  # all FALSE -> stop
  if(all(cond == FALSE)) {
    stop("the condition is not satisfied by any value")
  }
  
  # if incorrect length -> stop
  if(length(cond) != length(filter)) {
    stop("bad condition syntax")
  }
  
  # now filter the matrix
  if(byrow) {
    out <- data[cond, , drop = FALSE]
  } else {
    out <- data[, cond , drop = FALSE]
  }
  
  out
}

