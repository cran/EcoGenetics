# Local Moran's I statistic with Monte-Carlo test
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015


setGeneric("eco.localmoran",
      function(z, con, zerocon = NA, nsim = 99, 
           test = c("permutation", "bootstrap"), 
           alternative = c("auto", "two.sided", 
                   "greater", "less"),
           adjust = c("fdr", "holm", "hochberg", 
                 "hommel", "bonferroni", "BH",
                 "BY", "none")) {
       
       
       test <- match.arg(test)
       alternative.i <- match.arg(alternative)
       adjust <- match.arg(adjust)
       
       zerocon2 <- match(zerocon, c(0, NA))
       if(is.na(zerocon2)) {
        stop("zerocon argument must be 0 or NA")
       }
       

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
       
       
       n <- length(z)
       
       anselin <- function(z) {
       z2 <- z - mean(z)
       coef.sup <- t(replicate(length(z), return(z2))) 
       diag(coef.sup) <- 0
       coef.sup2 <- apply(coef.sup * wg, 1, sum)
       coef.inf <-  sum(z2 ^ 2) / n 
       out <- z2 * coef.sup2 / coef.inf
       out
       }

      obs <- anselin(z)
       
       
       if(test == "permutation") {
        
        monte.c <- replicate(nsim, anselin(sample(z)))
        monte.c <- t(monte.c)
        p.val <- numeric()
        expected <- numeric()
        altern <- rep("", ncol(monte.c))
        
        for (i in 1:ncol(monte.c)) {
         
         alternative <- alternative.i
         
         repet <- monte.c[, i]
         expected[i] <- median(repet)
         if(alternative == "auto") {
          alter <- expected[i] - obs[i]
          if(alter > 0) {
           alternative <- "less"
          } else if (alter <= 0) {
           alternative <- "greater"
          }
         }
         
         
         altern[i] <- alternative
         
         if(alternative == "greater") {
          howmany <- sum(repet >= obs[i])
          p.val[i] <- (howmany + 1) / (nsim + 1)
          
         } else if(alternative == "less") {
          howmany <- sum(repet <= obs[i])
          p.val[i] <- (howmany + 1) / (nsim + 1)
          
         } else if(alternative == "two.sided") {
          howmany <- sum(abs(repet) >= abs(obs[i]))
          p.val[i] <- (howmany + 1) / (nsim + 1)
         }
        }
       
        p.val <- p.adjust(p.val, method = adjust)
        
        tab <- data.frame(outer(obs, c(1:4)))
        tab[, 1] <- round(obs, 4)
        tab[, 2] <- round(p.val, 4)
        tab[, 3] <- round(expected, 4)
        tab[, 4] <- altern
        
        conect <- apply(wg, 1, sum)
        if(is.na(zerocon)) {
         tab[which(conect == 0), ] <- NA
        } else {
         tab[which(conect == 0), ] <- 0
        }
        
        colnames(tab) <- c("obs", "P", "exp", "alter")
        if(!is.null(rownames(z))) {
         rownames(tab) <- rownames(z)
        }
        
        
       } else if(test == "bootstrap") {
        
        monte.c <- replicate(nsim, anselin(sample(z, replace = TRUE)))
        
        if(!is.null(dim(monte.c))) {
        qq <- apply(monte.c, 1, quantile, probs = c(0.05, 0.95),  
              na.rm = TRUE)
       } else {
        qq <- quantile (monte.c, probs = c(0.05, 0.95),  na.rm = TRUE)
       }
      
      tab <- outer(obs, c(1:3))
      tab[, 1] <- round(obs, 4)
      tab[, 2:3] <- round(t(qq), 4)
      
      conect <- apply(wg, 1, sum)
      if(is.na(zerocon)) {
       tab[which(conect == 0), ] <- rep(NA, 3)
      } else {
       tab[which(conect == 0), ] <- rep(0, 3)
      }
      
      colnames(tab) <- c("obs", "lwr", "uppr")
      if(!is.null(rownames(z))) {
       rownames(tab) <- rownames(z)
      }
      
      tab <- as.data.frame(tab)
      
      }
      
      res <- list()
      
       global.moran <- round(sum(obs) / n, 4)
      
      if(nsim != 0) {
       if(test == "bootstrap") {
        res <- list("analysis" = "local-Moran", 
              "test" = test,
              "nsim" = nsim,
              "results" = tab,
              "global-Moran" = global.moran)
       } else {
        res <- list("analysis" = "local-Moran", 
              "test" = test,
              "nsim" = nsim,
              "p.adjust.method" = adjust,
              "results" = tab,
              "global-Moran" = global.moran)
       }
       
      } else {
       res <- list("analysis" = "local-Moran", 
             "results" = obs,
             "global-Moran" = global.moran)
      }
      class(res) <- "eco.gm"
      res
      
      
      })
      
