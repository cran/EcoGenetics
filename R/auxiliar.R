###############################
####$  AUXILIAR FUNCTIONS #####
###############################

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015 

# Union of classes


setClassUnion("dataframeORmatrix",
              c("data.frame", "matrix"))
setClassUnion("characterORnull",
              c("character", "NULL"))


# Conversion of a non symmetric binary matrix into symmetric

int.2symmetric <- function(mat) {
  mat[row(mat) > col(mat)] <- mat[row(mat) > col(mat)] + mat[row(mat) < col(mat)]
  mat[row(mat) > col(mat)][mat[row(mat) > col(mat)] > 0] <- 1
  mat[row(mat) < col(mat)] <- mat[row(mat) > col(mat)]
  mat
}


# Creates a matrix without diagonal, in row order

int.undimmattg <- function(mat) {
  ncolp <- ncol(mat) -1
  mat2 <- as.vector(t(mat))
  mat2 <-mat2[-which(col(mat) == row(mat))]
  mat2<-matrix(mat2, ncol = ncolp, byrow = T)
  mat2
}


# Checks a connection network. If "con" is of class "listw", 
# it returns a matrix

int.check.con <- function(con) {
  ccon <- class(con)[1]
  if(ccon == "listw") {
    listwg <- sapply(con$neighbours, c, simplify = FALSE)
    weig <- sapply(con$weights, c, simplify = FALSE)
    wg <- outer(z, z)
    wg[] <- 0
    for(i in 1:nrow(wg)) {
      wg[i, ][listwg[[i]]] <- weig[[i]]
    }
  } else if(ccon == "matrix"){
    wg <- con
  } else {
    stop("weight object provided is not of class listw or matrix")
  }
  wg
}


# Computing a distance matrix in meters among points in decimal degrees
# under a spherical Earth model


int.dlatlon2distm <- function(XY) {
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


# Scaling a data frame or matrix to 0 - 1 range

aue.rescale  <- function(dfm) {
  dfm <- as.data.frame(dfm)
  col <- apply(dfm, 2, function(X) { 
    (X - min(X, na.rm = TRUE)) / diff(range(X, na.rm = TRUE))
  })
  return(col)
}


# Ordering cells in a data frame or matrix

"aue.sort" <- function(x, ndig)  {
  
  geno <- as.matrix(x)
  
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
  rownames(m5) <- rownames(x)
  
  m5
  
}


# Converting data coded with characters into data coded by numbers

aue.char2num <- function(data, ndig = 0) {
  
  
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
}


#########################################################################
############ graphical functions ########################################
#########################################################################


### Graphical elements

# Circle perimeter

aue.circle <- function(mat, r, x0, y0, smooth = 100) { 
  
  for(theta in seq(0, 2*pi, pi/smooth)) {
    xd <- round(r*cos(theta)) 
    yd <- round(r*sin(theta)) 
    mat[x0 + xd, y0 + yd] <- 1
  }
  mat
}


# Solid circle

aue.point <- function(mat, r, x0, y0) { 
  xmat <- col(mat)
  ymat <- row(mat)
  mat2 <- mat - mat
  mat2[sqrt((xmat - x0)^2 + (ymat -y0)^2) <= r] <- 1
  mat2
}


# Solid diamond

aue.diamond<- function(mat, r, x0, y0) {
  xmat <- col(mat)
  ymat <- row(mat)
  mat2 <- mat - mat
  mat2[abs(xmat-x0) + abs(ymat -y0) <= r] <- 1
  mat2
}


# Solid square

aue.square <- function(mat, d, x0, y0) {
  xmat <- col(mat)
  ymat <- row(mat)
  mat2 <- mat - mat
  mat2[pmax(abs(xmat-x0), abs(ymat -y0)) <= d] <- 1
  mat2
}

### Other image manipulation functions

# Local filter

aue.filter <- function(mat, r, fun) {
  mresamp <- mat - mat
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      area <- aue.point(mresamp, r, i, j)
      tot <- sum(area != 0)
      mresamp[i, j] <- fun(mat * area)
    }
  }
  t(mresamp)
}


# Radial distance to a point

aue.circle.w <- function(mat, x0, y0) { 
  xmat <- col(mat)
  ymat <- row(mat)
  mat2 <- mat - mat
  mat2 <- sqrt((xmat - x0)^2 + (ymat -y0)^2)
  mat2
}


# Transforming a raster into a data frame 

aue.image2df <- function(mat) {
  xdim <- nrow(mat)
  ydim <- ncol(mat)
  mat.twocol <- data.frame(matrix(0, xdim * ydim, 3))
  colnames(mat.twocol) <- c("x", "y", "z")
  count <- 1
  for(j in 1:ydim) {
    for(i in 1:xdim) {
      mat.twocol[count, ] <- c(i, j, mat[i, j])
      count <- count + 1
    }
  }
  mat.twocol
}
