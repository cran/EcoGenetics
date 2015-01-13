#' Heterozygosity estimates with bootstrap confidence intervals over 
#' individuals and loci, and jackkninfe estimates over populations.
#' @param eco ecogen object.
#' @param grp Column of the slot S with the populations data.
#' @param nrep Number of repetitions for bootstrap resampling.
#' @param boot.indiv Logical. Should be performed bootstrap over individuals? 
#' Default TRUE.
#' @param boot.loci Logical. Should be performed bootstrap over loci? 
#' Default TRUE.
#' @param jack.pop Logical. Should be performed jackknife over populations? 
#' Default TRUE.
#' @param pres Logical. Should the output be in presentation mode? 
#' (see details). 
#' Default FALSE.
#' @details This program estimates  bootstrap confidence intervals for the 
#' observed and expected heterozygosity. The intervals could be computed for
#'  individual loci (bootstrap over individuals),
#' for mean values over loci (bootstrap over loci), 
#' and for mean values across populations (by jackknife). 
#' The program also could provide tables with the results formatted
#' in a presentation style.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @examples
#' \dontrun{
#' 
#' data(eco.test)
#' eco.boothet(eco,"structure",nrep = 20)
#' eco.boothet(eco,"structure",nrep = 20, pres = T)
#' 
#'}
#'
#' @export


setGeneric("eco.boothet", 
					 
					 function(eco, grp, nrep = 99, 
					 				 boot.indiv = TRUE, 
					 				 boot.loci = TRUE,
					 				 jack.pop = TRUE, 
					 				 pres = FALSE)  {                               #MISSING DATA AS NA
   
	
	inicio <- as.character(Sys.time())
	cat("\n\n", "analysis started at",inicio,"\n\n" )	
	
	
	if(pmatch(as.character(boot.indiv), c("TRUE", "FALSE"), nomatch= -1) == -1) {
		stop("boot.indiv must be TRUE or FALSE")
	}
	
	if(pmatch(as.character(boot.loci), c("TRUE", "FALSE"), nomatch= -1) == -1) {
		stop("boot.loci must be TRUE or FALSE")
	}
	
	if(pmatch(as.character(jack.pop), c("TRUE", "FALSE"), nomatch= -1) == -1) {
		stop("boot.indiv must be TRUE or FALSE")
	}
	
	if(pmatch(as.character(pres), c("TRUE", "FALSE"), nomatch= -1) == -1) {
		stop("boot.indiv must be TRUE or FALSE")
	}
	
	
  grupo <- eco$S
  fact <- match(grp,colnames(eco$S),nomatch = 0)
  fact <- fact[!fact == 0]
  
  if(length(fact) == 0)
    stop("incorrect group name")
  
  x <- eco.2hierfstat(eco,grp)
  colnames(x)[1]<-"pop"
  
  
  orden <- order(x[, 1])
  x <- x[orden,]
  nmar <- ncol(x[, -1])
  if(nmar < 2)
  {
    stop("too few markers!")
  }
  npop <- max(x[, 1])
  
  basicas <- function (data, diploid = TRUE, locus = FALSE) #VERSION SIMPLIFICADA DE BASIC.STATS
  {
    loc.names <- names(data)[-1]
    if (length(table(data[, 1])) < 2) 
      data[dim(data)[1] + 1, 1] <- 2
    p <- hierfstat::pop.freq(data, diploid)
    n <- t(hierfstat::ind.count(data))
    if (diploid) {
      dum <- hierfstat::getal.b(data[, -1])
      Ho <- dum[, , 1] == dum[, , 2]
      sHo <- (1 - t(apply(Ho, 2, fun <- function(x) tapply(x, data[, 1], mean, na.rm = TRUE))))
      mHo <- apply(sHo, 1, mean, na.rm = TRUE)
    } else {
      sHo <- NA
      mHo <- NA
    }
    sp2 <- lapply(p, fun <- function(x) apply(x, 2, fun2 <- function(x) sum(x^2)))
    sp2 <- matrix(unlist(sp2), nrow = dim(data[, -1])[2], byrow = TRUE)
    if (diploid) {
      Hs <- (1 - sp2 - sHo/2/n)
      Hs <- n/(n - 1) * Hs
    } else {
      Hs <- n/(n - 1) * (1 - sp2)
    }
    
    if(locus == TRUE)
    {
      np <- apply(n, 1, fun <- function(x) sum(!is.na(x)))
      mn <- apply(n, 1, fun <- function(x) {
        np <- sum(!is.na(x))
        np/sum(1/x[!is.na(x)])
      })
      msp2 <- apply(sp2, 1, mean, na.rm = TRUE)
      mp <- lapply(p, fun <- function(x) apply(x, 1, mean, na.rm = TRUE))
      mp2 <- unlist(lapply(mp, fun1 <- function(x) sum(x^2)))
      if (diploid) {
        mHs <- mn/(mn - 1) * (1 - msp2 - mHo/2/mn)
      } else {
        mHs <- mn/(mn - 1) * (1 - msp2)
      }
      
      res <- cbind(mHo, mHs)
      names(res) <- c("Ho", "Hs")
      if (diploid) {
        rownames(sHo) <- loc.names
      }
    }
    
    if(locus == TRUE)
    {
      return(res)} else {
        return(list(Ho = sHo, Hs = Hs))
      }
  }
  
  
  if(boot.indiv == T)
  {
    vec <- rep(0, npop)
    
    for(i in 1:npop)
    {
      vec[i] <- sum(x[, 1] == i)
    }
    
    lista <- list()
    lista.loci <- list(Ho = matrix(nrow = nmar, ncol = npop),
    									 He = matrix(nrow = nmar, ncol = npop))
    Hop <- matrix(nrow = nrep,ncol = nmar)
    Hep <- matrix(nrow = nrep,ncol = nmar)
    
    
    fun <- function(x, vec, i)
    {
      samp <- x[sample(filas <- which(x[, 1] == i) ,vec[i] , T), ]
      samp$pop <- rep(1, nrow(samp))
      return(samp)
    }
    
    repp <- function(nrep, x, vec, i){
    	replicate(nrep, eval(fun(x, vec, i), new.env()), simplify = F)
    }
    
    for(i in 1:npop)
    { 
      cat ("\r","bootstrap over individuals"," ","pop"," ",i," ",as.character(Sys.time()),"\n",sep="")
      muestras <- lapply(repp(nrep,x,vec,i),basicas)
      Hop <- sapply(muestras,function(x){return(x$Ho[,1])})
      Hep <- sapply(muestras,function(x){return(x$Hs[,1])})
      
      Hofinal <- Hefinal <- matrix(ncol = 2,nrow = nmar)
      Hofinal <- apply(Hop,1,quantile,probs = c(0.05, 0.95), na.rm = T)
      Hefinal <- apply(Hep,1,quantile,probs = c(0.05, 0.95), na.rm = T)
      
      x2 <- x
      x2[, 1] <-rep(1, nrow(x2))
      fac <- x[,1]
      temp <- basicas(x2[fac == i, ])     
      Ho <- temp$Ho[,1]
      He <- temp$Hs[,1]
      lista.loci[[1]][, i] <- Ho            
      lista.loci[[2]][, i] <- He           
      Hofinal <- t(rbind(Ho,Hofinal))
      Hefinal <- t(rbind(He,Hefinal))
      colnames(Hofinal) <-colnames(Hefinal) <- c("est", "lwr", "upr")
      rownames(Hofinal) <-rownames(Hefinal) <- colnames(x[, -1])
      lista[[i]]<-list("Ho" = Hofinal,"He" = Hefinal)
    }
    
  } else {
    lista <- NULL
  }
  
  if(boot.loci == T) {
    
    if(nmar < 4) {
      warning("warning, bootstrap over less than 4 
      				loci may be meaningless!", "\n")
    } 

    
    temp.Holoci <- temp.Heloci <- data.frame(matrix(nrow = nrep, 
    																								ncol = npop))
    for(i in 1:nrep)  
    {
      Holoci <- lista.loci[[1]][sample(1:nmar, nmar, T), ]
      Heloci <- lista.loci[[2]][sample(1:nmar, nmar, T), ]
      Holoci.p <- as.vector(apply(Holoci, 2, mean))
      Heloci.p <- as.vector(apply(Heloci, 2, mean))
      temp.Holoci[i,] <- Holoci.p
      temp.Heloci[i,] <- Heloci.p
      
      cat ("\r","bootstrap over loci "," ", 
      		 format(round(100*(i/nrep),1), trim = T,nsmall = 1), "% ", sep="")
    }
    
    cat("\n")
    qHo.media <- apply(temp.Holoci,2,quantile,probs = c(0.05, 0.95))
    qHe.media <- apply(temp.Heloci,2,quantile,probs = c(0.05, 0.95))
    Ho.media <- rbind(apply(basicas(x)$Ho, 2, mean)[1:npop], qHo.media)
    He.media <- rbind(apply(basicas(x)$Hs, 2, mean)[1:npop], qHe.media)
    colnames(Ho.media) <- colnames(He.media) <- c(1:npop)
    rownames(Ho.media) <- rownames(He.media) <- c("est", "lwr", "upr")
    
  } else {
    Ho.media <- NULL
    He.media <- NULL
  }
  
  if(jack.pop == T)
    
  {
    if(npop == 1)
    {
      warning("1 population for jackknife is not allowed!")
      Ho.pop <- NULL
      He.pop <- NULL
    } else {
      if(npop < 4) {
        warning("jackknife over less than 4 population
        				may be meaningless!", "\n")
      }

      
      
      ho.pop <- (basicas(x, locus = TRUE))[, 1]
      he.pop <- (basicas(x, locus = TRUE))[, 2]
      jac.ho <- data.frame(matrix(nrow = nmar, ncol = npop))
      jac.he <- data.frame(matrix(nrow = nmar, ncol = npop))
      
      for(i in 1:npop)
      {
      	jac<- basicas(x[x[, 1] != i, ], locus = TRUE)
        jac.ho[, i] <- jac[,1]
        jac.ho[, i] <- (npop * ho.pop) - ((npop - 1) * jac.ho[, i])
        jac.he[, i] <- jac[,2]
        jac.he [, i] <- (npop * he.pop) - ((npop - 1) * jac.he[, i]) 
        cat ("\r", "jackknife over populations "," ",
        		 format(round(100 * (i / nmar), 1), trim = T , nsmall = 1),
        		 "% ", sep = "")
      }
      
      cat("\n")
      
      mean.ho <- apply(jac.ho, 1, function(x) sum(x) / npop)
      mean.he <- apply(jac.he, 1, function(x) sum(x) / npop)
      
    
      desv.ho <- apply(jac.ho, 2, function(x) (x - mean.ho) ^ 2)
      desv.he <- apply(jac.he, 2, function(x) (x - mean.he) ^ 2)
      desv.ho <- apply(jac.ho, 1, sum)
      desv.he <- apply(jac.he, 1, sum)
      
      desv.ho <- (1 / (npop - 1)) * desv.ho
      desv.he<- (1 / (npop - 1)) * desv.he
      temp <- 1.96 * sqrt(desv.ho)
      Ho.pop <- cbind(ho.pop,mean.ho - temp,mean.ho + temp)
      temp <- 1.96 * sqrt(desv.he)
      He.pop <- cbind(he.pop,mean.he - temp,mean.he + temp)
      colnames(Ho.pop) <- colnames(He.pop) <- c("est", "lwr", "upr")
      
      
    }
  }

  resultados <- list()
  resultados$individual.bootstrap <- lista
  resultados$loci.bootstrap <- list(Ho = Ho.media, He = He.media)
  resultados$jackknife.pop <- list(Ho = Ho.pop, He = He.pop)
  
  if(pres==T)
   {
    
    bootsummary <- function(boot){
      
      
      simple <- function(tab){
        
        h <- format(round(tab, 3),nsmall = 3)
        apply(h, 1, function(y){paste(y[1]," ","(", y[2], "-", y[3], ")",
        															 sep = "")})
      }
      
      funtabla <- function(u,k) {
        if(k == 1) {
          a <- simple(u$Ho)
        } else {
          a <- simple(u$He)
        }
        return(a)
      } 
      
      dato <- boot$individual.bootstrap
      index <- 0
      for(i in 1:length(dato))
      {
        if(is.null(dato[[i]]$Ho)||is.null(dato[[i]]$He))
        {
          boot.ind <- NULL
          index <- 1
        }
      } 
      if(index == 0){
        boot.ind <- list(Ho = as.data.frame(sapply(dato, 
        																					 funtabla,1, 
        																					 simplify = T)), 
        								 He = as.data.frame(sapply(dato,funtabla, 
        								 													2,simplify = T)))
      }
      
      
      dato <- (boot$loci.bootstrap)
      if(is.null(dato$Ho)||is.null(dato$He))
      {
        boot.loci <- NULL
      } else {
        boot.loci <- list(Ho = as.data.frame(funtabla(dato, 1)), 
        									He = as.data.frame(funtabla(dato, 2)))
        colnames(boot.loci$Ho) <- "Ho"
        colnames(boot.loci$He) <- "He"
      }
      
      
      dato <- boot$jackknife.pop
      if(is.null(dato$Ho)||is.null(dato$He))
      {
        jack.pop <- NULL
      } else{
        jack.pop <- list(Ho = as.data.frame(funtabla(dato, 1)),
        								 He = as.data.frame(funtabla(dato, 2)))
        colnames(jack.pop$Ho)<-"Ho"
        colnames(jack.pop$He)<-"He"
      }
      
      return(list("boot.ind" = boot.ind, "boot.loci" = boot.loci,
      						"jack.pop" = jack.pop))
    
    }
    resultados <- bootsummary(resultados)
    cat("\n\n", "done!", "\n\n")
    resultados
  } else {
  	
  cat("\n\n", "done!", "\n\n")
  resultados
  }
  
})
