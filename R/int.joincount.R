# Join-count statistic, internal

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015 

int.joincount <- function(Z, con, ncod, nsim,
                          alternative, test = "permutation", 
                          adjust.n = FALSE, adjust) {
  
  con <- int.check.con(con)
  con <- as.vector(con)
  
  if(test == "permutation") {
    replace <- FALSE
  } else {
    replace <- TRUE
  }
  
  
  #internal function 
  
  jcfun <- function(input) {
    
    outmat <- outer(input, input, FUN = "paste", sep = "")
    if(is.null(ncod)) {
      stop("a ncod parameter was not given") 
    } else {
      outmat <- (as.matrix(aue.sort(outmat, ncod)))   #symmetric matrix
    }
    
    
    outmat <- as.factor(outmat)
    
    out <- numeric()
    temp <- list()
    
    for(i in seq(along = levels(outmat))) {
      temp[[i]] <- as.numeric(outmat == levels(outmat)[i])
    }
    
    out <- sapply(temp, function(x) sum(x * con) / 2) 
    
    out
  }
  
  obs <- jcfun(Z)
  
  #simulated datasets. permuting rows and columns, mantaining structure
  monte.c <- matrix(0, nrow = nsim, ncol = length(obs))
  for(i in 1:nsim) {
    samp <- sample(length(Z), replace = replace)
    outsamp <- Z[samp]
    monte.c[i, ] <- jcfun(outsamp)
  }
  
  ran <- int.random.test(repsim = monte.c, obs = obs, 
                         nsim = nsim, test = test,
                         alternative = alternative,
                         adjust = adjust)
  
  
  #labeling rows
  outmat <- outer(Z, Z, FUN = "paste", sep = "")
  outmat <- (as.matrix(aue.sort(outmat, ncod))) 
  outmat <- as.factor(outmat)
  rownames(ran) <- levels(outmat)
  
  res <- list("analysis" = "Join-count", 
              "nsim" = nsim,
              "results" = ran)
  
  res
  
}
