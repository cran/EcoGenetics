#' Moran's I, Geary's C and Join-Count  correlograms for an ecogen object.
#' @param eco ecogen object.
#' @param int Interval distance for the analysis.
#' @param smax Maximum distance for the analysis.
#' @param nsim Number of Monte - Carlo simulations.
#' @param select The analysis type ("moran", "geary", "join-count").
#' @param d.class The class of data to analyze ("alleles", "locus",
#' "phenotype").
#' @param style Weights matrix style passed to \code{\link[spdep]{nb2listw}}.
#' Can take  the values "W", "B", "C", "U", "minmax" and "S".
#' @param fact The name of the S slot column when analysis per-group is required. 
#' Default analysis takes in account all sites.
#' @param zero.policy Default TRUE assign zero to the lagged value
#' of zones without neighbors, if FALSE assign NA, if NULL use global option 
#' value.
#' @param grp The group (contained in x) to analyze if analysis
#' per-group is required. Default analysis takes in account all sites.
#' @description This program can manipulate different class of data. 
#' For data of class "allele" the program takes the data frame allocated
#' in the GENIND slot. For data of class "locus", the program takes
#' the data frame of the G slot and transforms the genotypes 
#' of the individual loci into factor levels (e.g., for a biallelic loci
#' with alleles A and B, three levels could be present in the data:
#' "AA","AB","BB"). Before the "locus" analysis, it is necessary
#' to run \code{\link{eco.sortalleles}} in codominant diploid
#'  data if the alleles are not ordered (e.g., 15 and 51 will 
#'  be interpreted as two levels).  For data of class "phenotype"
#'   the data frame used is the allocated in the P slot.
#' For "allele" and "phenotype" data, "moran", "geary" and
#' "join-count" analysis are available (the user must match
#' their data with the  statistic, according to the class of analysis,
#' data is taken as factor with join-count). 
#' For "locus" data, only join-count is available. 
#' @return An object of class "RelDist" with data frames for the statistic
#' and the 95\% upper and lower confidence interval bounds estimated by
#' bootstrap. 
#' @seealso  \code{\link[spdep]{moran}} \code{\link[spdep]{geary}}
#' \code{\link[spdep]{joincount.test}}
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @examples 
#' 
#' \dontrun{
#' data(eco.test)
#' eco.ac <- eco.autocor(eco, int = 50, smax = 1000, d.class = "phenotype")
#' plot(eco.ac, var = 6)
#' plot(eco.ac, var = 3)
#' 
#' }
#' @export

setGeneric("eco.autocor",  function(eco, fact = NA, 
																		int, 	smax,
																		nsim = 99, 	select = c("moran", "geary", "join-count"),
																		d.class = c("alleles", "locus", "phenotype"),	grp = NA,
																		style = c("B", "W", "C", "U", "minmax", "S"),	
																		zero.policy = TRUE) {
				
	
	
  cat("\n")
  cat(" analysis started at", as.character(Sys.time()),"\n")
  cat("\n")
  
  select <- match.arg(select)
  d.class <- match.arg(d.class)
  style <- match.arg(style)
  if(pmatch(as.character(zero.policy), c("TRUE", "FALSE"), nomatch= -1) == -1) {
  	stop("zero.policy must be TRUE or FALSE")
  }
  
  if(!is.na(fact)) {
    grupo <- eco@S
    fact <- match(fact, colnames(eco@S), nomatch = 0)
    fact <- fact[fact != 0]
    
    if(length(fact) == 0) {
      stop("incorrect factor name")
    }
    
    if(is.na(grp)) {
    	stop("grp argument not provided")
    }
    
    where <- which(eco@S[, fact] == grp)
    hmuch <- sum(dist(eco@XY[where, ]) < int)
    if(hmuch < 5) {
      stop("Scale not apropiated.Increase distance interval")
    }
    
    hlast <- sum(dist(eco@XY[where, ]) > smax - int)
    if(hlast <5 ) {
      stop("Range not apropiated. Decrease smax value")
    }
    
    M <- as.matrix(eco@P[where, ])
    G <- as.matrix(eco@G[where, ])
    xy <- as.matrix(eco@XY[where, 1:2])
    if(ncol(xy)>2) {
    	cat("The XY data frame contains more that 2 columns.
(maybe altitude data, but it is ok). The program takes the 
first two columns as latitude -longitude information.", "\n\n")
    }
  } else {
    hmuch <- sum(dist(eco@XY) < int)
    if(hmuch < 5) {
      stop("Scale not apropiated.Increase distance interval")
    }
    hlast <- sum(dist(eco@XY) > smax - int)
    if(hlast < 5) {
      stop("Range not apropiated. Decrease smax value")
    }
    
    M <- as.matrix(eco@P)
    G <- as.matrix(eco@G)
    xy <- as.matrix(eco@XY[,1:2])
    if(ncol(xy)>2) {
    	cat("The XY data frame contains more that 2 columns.
(maybe altitude data, but it is ok). The program takes the 
first two columns as latitude -longitude information.", "\n\n")
    }
  }
  
  salida <- new("RelDist")
  
  
  
  analysis = c("moran", "geary", "join-count")
  test = which(analysis %in% select == TRUE)
  if(length(test) == 0) {
    stop("analysis type do not match with moran, geary or join-count")
  }
  
  if(d.class == "alleles") {
    x <- eco@GENIND$tab
  } else if(d.class == "phenotype") {
    x <- M
  } else if(d.class == "locus") {
    if(!is.factor(G[, 1])) {
      if(select != "join-count") {
        stop("only join-count method is available with factor data")
      }
      
      for(i in 1:ncol(G)) {
        x <- G
        x[, i] <- as.factor(x[, i])
      }
    } 
  } else {
    stop("please enter a correct variable type (alleles, locus, phenotype)")
  }
  
  
  nloc <- ncol(x)
  d.max<- seq(int, smax, int)
  d.min <- d.max - int
  d.min[1] <- 1e-10
  j <- 0
  
  
  
 ####-Computing distance intervals-###
  
  
  medint <- function(coordenadas, int, smax) {
    
    dist.range <- function(coordenadas, int, smax)
    {
      distancia <- dist(coordenadas)
      lista <- list()
      j <- 1
      for (i in seq(int, smax, int)) {
        temp <- which((distancia <= i) & (distancia > i - int))
        lista[[j]] <- temp
        names(lista)[j] <- i
        j <- j + 1
      }
      lista
    }
    rangos <- dist.range(coordenadas, int, smax)
    listamedias <- vector()
    dis <- dist(coordenadas)
    for(i in 1:length(rangos)) {
      listamedias[i] <- mean(dis[rangos[[i]]])
    }
    listamedias
  }
  
  
 ####-Funtion to estimate the stat in each iteration-###
  
  choose_stat <- function(...) {
    
    terms <- list(...)
    
    var.t <- terms$var
    test.t <- terms$test
    n.t <- terms$n
    n1.t <- terms$n1
    S0.t <- terms$S0
    zero.policy.t <- terms$zero.policy
    listw.t <- terms$listw
    
    if(test.t == 1) {
      out <- spdep::moran(x = var.t, n = n.t, S0 = S0.t,
                   zero.policy = zero.policy.t, listw = listw.t)$I
    } else if(test.t == 2) {
      out <- spdep::geary(x = var.t, n = n.t, n1 = n1.t, S0 = S0.t,
                   zero.policy = zero.policy.t, listw = listw.t)$C
    } else if(test.t == 3) {
      out <- spdep::joincount.test(fx = as.factor(var.t), 
                            zero.policy = zero.policy.t,
                            listw = listw.t)[[1]]$statistic 
    }
    out
  }
  
  
 ####-Function that iterate the latter-###
  
  estimate <- function(dat, xy, d.min, d.max, nsim) { 
    
    lista<-list()
    niter<-ncol(dat)*length(d.min)
    dist.dat<-paste("d=", d.min, "-", d.max)
    dist.dat[1]<-paste("d=","0", "-", d.max[1])
    
    for(k in 1:ncol(dat)) {
      lista[[k]] <- matrix(0, nrow = length(d.min), ncol = (nsim + 1))
      rownames(lista[[k]])<-paste("d=", d.min, "-", d.max)
    }
    
    for(i in 1:length(d.min))
    {
      con <- spdep::dnearneigh(sp::coordinates(xy), d1 = d.min[i], d2 = d.max[i])
      con <- spdep::nb2listw(con, zero.policy = TRUE, style = style)
      wc <- spdep::spweights.constants(con, zero.policy = TRUE, adjust.n = TRUE)
      n <- wc$n
      n1 <- wc$n1
      S0 <- wc$S0
      
      
      for(j in 1:ncol(dat)) {
        
        var<-dat[,j]
        
        lista[[j]] [i, 1:nsim] <- replicate(nsim, 
        																	choose_stat(var = sample(var), 
                                                     test = test,
                                                     n = n, n1 = n1, 
        																						 S0 = S0,
                                                     zero.policy = zero.policy,
                                                     listw = con))
        lista[[j]] [i, nsim+1] <- choose_stat(var = var, test = test,
                                            n = n, n1 = n1, S0 = S0,
                                            zero.policy = zero.policy,
        																			listw = con)
        
        cat("\r", "Computing", "class distance",dist.dat[i]," ",
            ceiling(j*100/ncol(dat)), "%")
      }
      cat("\n")
    }
    lista
  }
  
  ###- Calling the latter-###
  
  res <-estimate(x, xy, d.min, d.max, nsim)
  
  
  ls.test <- list()
  tabla <- matrix(, length(d.min), 5)
  tabla[,1] <- medint(xy, int, smax)
  colnames(tabla) <- c("d.mean", "est", "sd", "lwr", "upr")
  rownames(tabla) <- paste("d=", d.min, "-", d.max, sep = "")
  
  for(i in 1:ncol(x))
  {
    for(j in 1:length(d.min))
    {
      sta <- c(res[[i]][j, nsim+1], sd(res[[i]][j, ]), quantile(res[[i]][j, ],
                                                                probs = c(0.05, 0.95),  na.rm = TRUE))
      tabla[j, 2:5] <- sta
    }
    ls.test[[i]] <- tabla
  }
  names(ls.test) <- colnames(x)
  cat("\n")
  
  
  salida@OUT <- ls.test
  salida@NAMES <- colnames(x)
  salida@INTERVAL <- int
  salida@MAX <- smax
  salida@TYPE <- d.class
  salida@METHOD <- select
  assign(paste(substitute(eco),"@OUT$RELGEN",sep=""), salida)
  
  cat("\n")
  cat("done!")
  cat("\n\n")
  
  ReturnVal <- tcltk::tkmessageBox(title = "Correlogram", 
                                   message = "process successful!",
                                   icon = "info", type = "ok")
  
  salida
 
})
