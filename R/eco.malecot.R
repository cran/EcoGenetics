# Leandro Roser leandroroser@ege.fcen.uba.ar
# June 17, 2015 


# Global and local kinship analyses (beta version)

eco.malecot <- function(eco, 
                         method = c("global", "local"), 
                         kinmatrix =NULL,
                         int = NULL,
                         smin = 0,
                         smax = NULL,
                         nclass = NULL,
                         kmax = NULL,
                         seqvec = NULL,
                         size = NULL,
                         type = c("knearest", "radialdist"),
                         cubic = TRUE,
                         testclass.b = TRUE,
                         testmantel.b = TRUE,
                         jackknife = TRUE,
                         cummulative = FALSE,
                         nsim = 99, 
                         test = c("permutation", "bootstrap"), 
                         alternative = c("auto","two.sided", 
                                         "greater", "less"),
                         sequential = TRUE, 
                         conditional = c("AUTO", "TRUE", "FALSE"),
                         bin = c("sturges", "FD"),
                         row.sd = FALSE,
                         adjust = "holm",
                         latlon = FALSE) {
  
  XY <- eco$XY
  
  if(ncol(XY)>2) {
    message("XY with > 2 columns. The first two are taken as X-Y coordinates")
    XY <-XY[,1:2]
  } 
  
  if(latlon == TRUE) {
    XY <- SoDA::geoXY(XY[,2], XY[,1], unit=1)
  } 
  
  distancia <- dist(XY)
  logdistancia <- log(distancia)
  
  control <- c(!is.null(smax), !is.null(kmax), !is.null(seqvec))
  if(sum(control) == 0) {
    smax <- max(distancia)
  } 
  if(sum(control) > 1) {
    stop("multiple selection of smax, kmax, seqvec; please enter only one parameter")
  }
  conditional <- match.arg(conditional)
  method <- match.arg(method)
  type <- match.arg(type)
  test <- match.arg(test)
  alternative <- match.arg(alternative)
  bin <- match.arg(bin)
  geno <- eco$GENIND$tab
  loc.fac <- as.numeric(eco$GENIND$loc.fac)
  nloc <- max(loc.fac)
  nind <- nrow(eco$GENIND$tab)
  crit <- abs(qt(0.975, nloc - 1)) #critical value for CI
  
  #kinship matrix
  if(is.null(kinmatrix)) {
  kin <- int.kin.loiselle(geno, loc.fac)
  } else {
    kin <- as.matrix(kinmatrix)
    diag(kin) <- 0
  }
  
  #tests
  
  if(test == "permutation") {
    replace <- FALSE
  } else {
    replace <- TRUE
  }
  
  #method selection
  if(method == "global") {
    type <- "correlogram"
    
    if(is.null(smax) & is.null(seqvec)) {
      stop("please provide a smax argument or a vector of classes for dist method")
    }
    conditional <- FALSE
    modelmatrix.comp <- eco.lagweight(XY, 
                                      int = int, 
                                      smin = smin,
                                      smax = smax, 
                                      nclass = nclass,
                                      size = size,
                                      seqvec = seqvec,
                                      row.sd = row.sd,
                                      bin = bin,
                                      cummulative = FALSE)
    
    modelmatrix <- modelmatrix.comp@W
    nbins <- length(modelmatrix)
      
  } else if(method == "local") {
    
    sequential <- FALSE
    
    if(conditional == "AUTO") {
      conditional <- FALSE
    }
    
    if(type == "knearest" & is.null(kmax)) {
      stop("kmax argument must be non-null with local knearest analysis")
    }
    if(type == "radialdist" & is.null(smax)) {
      stop("smax argument must be non-null with local radialdist analysis")
    }
    if(type == "knearest") {
      modelmatrix <- eco.weight(XY = XY, 
                                method = "knearest", 
                                k = kmax,
                                row.sd = row.sd)
      modelmatrix <- modelmatrix@W
    } else if (type == "radialdist") {
      modelmatrix <- eco.weight(XY = XY, 
                                method = "circle", 
                                d1 = 0, 
                                d2 = smax,
                                row.sd = row.sd)
      modelmatrix <- modelmatrix@W
    }
    nbins <- nind
  }
  
  
  ### function for computing each case
  
  select_method <- function(kinmat) {
    out <- numeric()
    
    if(method != "local") {
      
      for(i in 1:nbins) {
        if(sum(modelmatrix[[i]]) == 0) {
          stop("no individuals in class", i)
        }
        out[i] <- sum(kinmat * modelmatrix[[i]]) / sum(modelmatrix[[i]])
      }
      
    } else {
      out <- apply((kinmat * modelmatrix), 1, sum) / apply(modelmatrix, 1, sum)
      out[is.infinite(out)] <- NA             #when there are no connections
      out[is.na(out)] <- NA 
    }
    out
  }
  
  obs <- select_method(kin)
  
  # Here stats the permutation test. The kinship matrix is randomized, 
  #then select_method is called each time
  
  random.kin <- matrix(0, nbins, nsim)
  samp <- 1:nind
  kin.perm <- kin
  cat("\n\n")
  
  counter <- 1
  
  
  for(i in 1:nsim) {
    
    #permuted matrix for each repetition
    #conditional case
    if(conditional) {
      for(k in samp) {
        order.z <- sample(samp[-k], replace = replace)
        kin.perm[k,-k] <- kin.perm[k, ][order.z]
      }
      
      #free case
    } else {
      shuffle.kin <- sample(samp, replace = replace)
      kin.perm <- kin[shuffle.kin, shuffle.kin]
    }
    
    
    random.kin[ , i] <- select_method(kin.perm)
    cat(paste("\r", "simulations...computed",
              100 * round(counter / nsim, 1), "%"))
    counter <- counter + 1
  }
  
  
  if(method == "global") {
    meandistance <- modelmatrix.comp@MEAN
    logdist <- modelmatrix.comp@LOGMEAN
    breaks.kin <- modelmatrix.comp@BREAKS
    cardinal <- modelmatrix.comp@CARDINAL
    if(!cummulative) {
      d.max <- round(breaks.kin[-1], 3)
      d.min <- round(breaks.kin[-length(breaks.kin)], 3)
      dist.dat<-paste("d=", d.min, "-", d.max, sep = "")
    } else {
      dist.dat <- paste("d=0-d=", breaks.kin[-1], sep="")
    }
    
  } else if(method == "local") {
    sequential <- FALSE
    distancia <- as.matrix(distancia)
    meandistance <- apply(modelmatrix * distancia, 1, sum) / apply(modelmatrix, 1, sum)
    cardinal  <- apply(modelmatrix, 1, sum)
    logdist <- numeric()
    breaks.kin <- 0
    d.max <- 1
    d.min <- rep(1, nrow(modelmatrix))
    dist.dat<- rownames(XY)
  }
  
  
  if(test == "permutation") {
    tab <- data.frame(matrix(0, nrow = length(dist.dat), ncol= 13))
    rownames(tab) <- dist.dat
    colnames(tab) <- c("d.mean","d.log","obs", "exp","alter", 
                       "p.val", "mean.jack", "sd.jack", "Jack.CI.inf",
                       "Jack.CI.sup", "null.lwr", "null.uppr", 
                       "cardinal")
    tab[, 1] <- meandistance
    tab[, 2] <- logdist
    for(i in 1:nrow(tab)) {
      ran <-  int.random.test(random.kin[i, ], obs[i],
                              test = "permutation", 
                              alternative = alternative, 
                              nsim = nsim) 
      
      if(!is.na(ran[[1]])) {
        
        tab[i, 3] <- round(ran[[1]], 4) #obs
        tab[i, 4] <- round(ran[[2]], 4) #exp
        tab[i, 5] <- ran[[3]]           #alter
        tab[i, 6] <- round(ran[[4]], 4) #p.val
        tab[i, 11] <- round(ran[[5]][1], 4) #null.lwr
        tab[i, 12] <- round(ran[[5]][2], 4) #null.uppr
        
        tab[, 3:4] <- round(tab[, 3:4], 4) # rounding obs, exp
        tab[, 13] <- cardinal       #cardinal
        
        if(sequential) {
          for(i in 1:nrow(tab)) {
            tab[i, 6] <- (p.adjust(tab[1:i, 6], method= adjust))[i]
          }
        } else {
          tab[ , 6] <- p.adjust(tab[ , 6], method = adjust)
        }
        tab[, 6] <- round(tab[, 6], 5)
        
      } else  {
        tab[i, 3] <- NA
        tab[i, 4] <- NA
        tab[i, 5] <- NA
        tab[i, 6] <- NA
        tab[i, 9] <- NA
        tab[i, 10] <- NA
      }
      
    } 
  } else if(test == "bootstrap") {
    tab <- data.frame(matrix(nrow = length(d.min), ncol=10))
    rownames(tab) <- dist.dat
    colnames(tab) <- c("d.mean","d.log", "obs", "mean.jack", "sd.jack", 
                       "Jack.CI.inf", "Jack.CI.sup", "null.lwr", "null.uppr",  
                       "cardinal")
    tab[, 1] <- round(meandistance, 3)
    tab[, 2] <- round(logdist, 3)
    for(i in 1:nrow(tab)) {
      rand <-  int.random.test(random.kin[i, ],
                               obs[i], 
                               test = "bootstrap", 
                               nsim = nsim)
      if(!is.na(obs[i])) {
        tab[i, 3] <- round(rand[[1]], 4) # obs
        tab[i, 8:9] <- round(rand[[2]], 4) #null lwr, uppr
      }
    }
    tab[, 10] <- cardinal     #cardinal
  }
  
  
  
  #SD AND SP ANALYIS FOR "global"
  
  if(method == "global") {  
    ##############calculo de sp 
    dist.sp <- as.matrix(distancia)
    logdistmat <- as.matrix(logdistancia)
    Fij.sp <- kin
    dmin <- min(breaks.kin)
    dmax <- max(breaks.kin)
    if(dmin == 0) {
    dmin <- exp(-100)
    }
    restricted <- (dist.sp >= dmin & dist.sp <= dmax)
    logdist.sp <- logdistmat[restricted]
    Fij.sp <- Fij.sp[restricted]
    
    modelo <- lm(Fij.sp ~ logdist.sp)
    bhat <- as.numeric(coef(modelo)[2])
    sp.stat <- - bhat /(1 - (tab$obs)[1])
    x.intercept <- exp(1)^(-coef(modelo)[1]/coef(modelo)[2])
    
    ###cubic interpolation
    if(cubic) {
    Fij.detrended <- modelo$residuals
    cubica <- lm(Fij.detrended ~ poly(logdist.sp, degree = 3))
    capture.output(cubica.final <- step(cubica, scope = list(lower = Fij.detrended ~ 1,
                                                             upper = cubica)))
    }
    
    ####permutation test for b. Permutation over locations. 
    #if(testclass.b) {
    #cat("\n")
    #bhat.test <- rep(0, nsim)
    #for(k in 1:nsim) {
    #  samp <- sample(nrow(kin))
    #  kin2 <- kin[samp, samp]
    #  obs.sim <- numeric()
    #  for(i in 1:nbins) {
    #    if(sum(modelmatrix[[i]]) == 0) {
    #      stop("no individuals in class", i)
    #    }
    #    obs.sim[i] <- sum(kin2 * modelmatrix[[i]]) / sum(modelmatrix[[i]])
    #  }
    #  modelo.sim <- lm(kin2[restricted] ~ logdist.sp)
    #  bhat.sim <- as.numeric(coef(modelo.sim)[2])
    #  bhat.test[[k]] <- bhat.sim
    #  cat(paste("\r", "testing slope...computed",
    #            100 * round(k /nsim, 1), "%"))
    #}
    
    #btest <- int.random.test(bhat.test, bhat, nsim)
    #}
    
    #--------------------MANTEL TEST FOR SP-------------------------------------#
    
    if(testmantel.b) {
    mantelcells <- (row(logdistmat)>col(logdistmat)) & restricted
    logdist.mantel <- logdistmat[mantelcells]
    Fij.mantel <- kin[mantelcells]
    obs <- cor(logdist.mantel, Fij.mantel)
    
    repsim <- numeric()
    N <- nrow(logdistmat)
    for(i in 1:nsim){
      samp <- sample(N)
      temp <- (kin[samp, samp])[mantelcells]
      repsim[i] <- cor(logdist.mantel, temp)
    }
    
    resmantel <- int.random.test(repsim = repsim, 
                                 obs = obs,
                                 nsim = nsim,
                                 test = "permutation",
                                 alternative = "auto")
    mantel.obs <- resmantel$obs
    mantel.pval <- resmantel$p.val
    }
    
    ##################################################################
    if(jackknife) {
    #jackknife over loci
    obs.jack <- matrix(0, nrow = nloc, ncol = nbins)
    
    
    #jackknifing kinship per class-----------#
    cat("\n")
    kin.loci <- list()
    nrepet <- nloc * nbins
    counter <- 1
    for(i in 1: nloc) {
      col.delete <- which(loc.fac == i)
      loc.fac.jack <- as.numeric(as.factor(loc.fac[-col.delete]))  # renumber the factor from 1 to nloc-1
      geno.jack <- geno[, -col.delete]
      kin.jack <- int.kin.loiselle(geno.jack, loc.fac.jack)
      kin.loci[[i]] <- kin.jack
      for(j in 1:nbins) {
        obs.jack[i, j] <- sum(kin.jack * modelmatrix[[j]]) / sum(modelmatrix[[j]])
        cat("\r", paste("jackknife over loci...", 100 * round(counter / (nrepet), 1), "%"))
        counter <- counter + 1
      }
    }
    cat("\n")
    #---------#
    
    ### error bars
    obs.real <- t(replicate(nloc, tab$obs))                                                     #value without jackknife                                                        #mean of jackknife values
    pseudo <- (nloc * obs.real) - ((nloc - 1) * obs.jack)                   #pseudo-value
    theta <- apply(pseudo, 2, mean)
    pseudo.variance <- apply(pseudo, 2, var)
    sd.jack <- sqrt(pseudo.variance / nloc) 
    ####
    #F1 CONFIDENCE INTERVALS
    temp.F1 <- crit * sd.jack[1]
    F1.confint <- drop(cbind(theta[1] - temp.F1,theta[1] + temp.F1))
    ###
    #ALL CONFIDENCE INTERVALS
    temp.FCI <- crit * sd.jack
    FCI <- drop(cbind(theta - temp.FCI, theta + temp.FCI))
    colnames(FCI) <- c("Jack.lwr", "Jack.uppr")
    if(test == "permutation") {
      tab[, 7] <- round(theta, 4)
      tab[, 8] <- round(sd.jack, 4)
      tab[, 9:10] <- round(FCI, 4)
    } else {
      tab[, 4] <- round(theta, 4)
      tab[, 5] <- round(sd.jack, 4)
      tab[, 6:7] <- round(FCI, 4)
    }
    
    #---------------------------------------------------------------------------#
    #JACKKNIFE OF THE SLOPE OVER LOCI, USING INDIVIDUALS 
    bhat.jack <- numeric()
    for(i in 1:nloc) {
      mod.jack <- lm(kin.loci[[i]][restricted] ~ logdist.sp)
      bhat.jack[i] <- as.numeric(coef(mod.jack)[2])
    }
    pseudo.bhat <- (nloc * rep(bhat, nloc)) - ((nloc - 1) * bhat.jack)                   #pseudo-value
    theta.bhat <- mean(pseudo.bhat)
    pseudovar.bhat <- var(pseudo.bhat)
    sd.bhat <- sqrt(pseudovar.bhat / nloc)
    temp <- crit * sd.bhat
    bhat.confint <- drop(cbind(theta.bhat - temp,theta.bhat + temp))
    sp.confint <- - bhat.confint /(1 - (tab$obs)[1])
    sp.confint <- drop(c(sp.confint[2], sp.confint[1]))
    names(sp.confint) <- c("5%", "95%")
    
    bhat <- round(c(bhat, sd.bhat, theta.bhat, bhat.confint[1],
                     bhat.confint[2], x.intercept, 
                     (tab$obs)[1], F1.confint[1], F1.confint[2]), 4)
    
    names(bhat) <- c("bhat", "SD","theta", "CI.5%", "CI.95%", 
                      "X-intercept", "F1", "F1.CI.5%", "F1.CI.95%")
    
    sp.stat <- c(sp.stat, sp.confint[1], sp.confint[2])
    
    names(sp.stat) <- c("sp", "CI.5%", "CI.95%")
    
    }
    
    distance <- data.frame(meandistance, logdist)
    colnames(distance) <- c("dist", "log.dist")
    rownames(distance) <- 1:nrow(distance)
    distance <- t(distance) 
    
  
    out.sp <- list(model = modelo,
                   distance = distance,
                   restricted = round(c(dmin, dmax), 4),
                   bhat = bhat,
                   mantel.obs.b = round(mantel.obs, 4),
                   mantel.pval.b = round(mantel.pval, 4),
                   sp = round(sp.stat, 5), 
                   cubic_model = ifelse(cubic, cubica.final, NA)
    )
    

  #######################
    
  }
  
  if(method == "global") {
  
  salida <- new("eco.IBD")
  salida@OUT <- list(tab)
  salida@IN <- list(XY = XY, ECOGEN = eco$GENIND)
  salida@BREAKS <- breaks.kin
  salida@CARDINAL <- cardinal
  salida@NAMES <- colnames(eco$G)
  salida@METHOD <- c("Kinship", type)
  salida@DISTMETHOD <- method
  salida@TEST <- c(test, conditional) 
  salida@NSIM <- nsim
  salida@PADJUST <- paste(adjust, "-sequential:", sequential)
  salida@SP <- out.sp
  
  } else if (method == "local") {
    
    salida <- new("eco.lsa")
    if(test == "permutation") {
      tab <- tab[, c(1,3,4,5,6,11,12,13)]
    } else {
      tab <- tab[, c(1,3,8, 9,10)]
    }
    salida@OUT <- tab
    salida@METHOD <- "local kinship"
    salida@TEST <- test
    salida@NSIM <- nsim
    salida@COND <- ifelse(conditional ==  "TRUE", TRUE, FALSE)
    salida@PADJ <- adjust
    salida@XY <- data.frame(XY)
  }
  
  cat("\n")
  cat("done!")
  cat("\n\n")
   
  salida
}
