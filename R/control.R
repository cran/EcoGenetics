# checking weights

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015 

int.check.con <- function(con) {
	
	ccon <- class(con)[1]
	
	if(ccon == "listw") {
		listwg <- sapply(con$neighbours, c, simplify = FALSE)
		weig <- sapply(con$weights, c, simplify = FALSE)
		Z<- 1:length(con$weights)
		wg <- outer(Z, Z)
		wg[] <- 0
		for(i in 1:nrow(wg)) {
			wg[i, ][listwg[[i]]] <- weig[[i]]
		}
	} else if(ccon == "matrix"){
		wg <- con
	} else if(ccon == "eco.weight"){ 
	  wg <- con@W
	} else {
		stop("weight object provided is not of class listw, matrix or eco.weight")
	}
	wg
}


# checking if columns are in numeric format in matrix/ data frame objects

int.check.numeric <- function(mat) {
  
  x <- mat
  clases <- character()
  for(i in 1:ncol(x)) {
    clases[i] <- class(x[, i])
  }
  
if(any(clases != "numeric" | clases != "integer")) {
  x <- as.matrix(x)
  colhier <- ncol(x)
  rowhier <- nrow(x)
  x <- matrix(as.numeric(x), ncol = colhier, nrow= rowhier)
  if(class(mat) == "data.frame") {
  x <- as.data.frame(x)
  }
  colnames(x) <- colnames(mat)
  rownames(x) <- rownames(mat)
  
  x
}
}
