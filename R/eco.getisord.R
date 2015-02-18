# Local Getis-Ord's G and G* for analysis of hot spots
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015



setGeneric("eco.getisord",
       function(z, con, zerocon = NA, 
                classG = c("G", "G*"),
                nsim = 99, 
                test = c("permutation", "bootstrap"),
                alternative = c("auto", "two.sided", 
                                "greater", "less"), 
               adjust = c("fdr", "holm", "hochberg", 
                          "hommel", "bonferroni", "BH",
                          "BY", "none")) {
        
        classG <- match.arg(classG)
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
        
        Gm <- t(replicate(n, return(z)))
        
        getis <- function(Gm) {
        
        if (classG == "G") {
         
         diag(wg) <- 0
          diag(Gm) <- 0
        
        W <- apply(wg, 1, sum)
        S1 <- apply(wg ^ 2, 1, sum)
        
        G <- wg * Gm
        G <- apply(G, 1, sum)
        Gm2 <- Gm ^ 2
        
        meanG <- apply(Gm, 1, sum) / (n -1)
        desvsqG <- apply(Gm2, 1, sum) / (n -1)
        desvsqG <- desvsqG - meanG ^ 2
        W <- apply(wg, 1, sum)
              S1 <- apply(wg ^ 2, 1, sum)
        denom <- (n - 1) * S1 - W ^ 2
        denom <- sqrt(desvsqG * denom / (n - 2))
        numer <- (G - W * meanG)
        G <- numer / denom

        } else if(classG == "G*") {
         
         W <- apply(wg, 1, sum)
         S1 <- apply(wg ^ 2, 1, sum)
         
         G <- wg * Gm
         G <- apply(G, 1, sum)
         Gm2 <- Gm ^ 2
         
         meanG <- apply(Gm, 1, sum) / n
         desvsqG <- apply(Gm2, 1, sum) / n
         desvsqG <- desvsqG - meanG ^ 2
       
         denom <- n * S1 - W ^ 2
         denom <- sqrt(desvsqG * denom / (n - 1))
         numer <- (G - W * meanG)
         G <- numer / denom
        }
        G
        
        } 
        
        obs <- getis(Gm)
        
        if(test == "permutation") {
         
          
         monte.c <- replicate(nsim, 
                    getis(t(replicate(n, 
                             return(sample(z))))))
         
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
           } else if (alter < 0) {
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
          tab[which(conect == 0), ] <- rep(NA, 4)
         } else {
          tab[which(conect == 0), ] <- rep(0, 4)
         }
         
         colnames(tab) <- c("obs", "P", "exp", "alter")
         if(!is.null(rownames(z))) {
          rownames(tab) <- rownames(z)
         }
         
        } else if(test == "bootstrap") {
         
         monte.c <- replicate(nsim, 
                    getis(t(replicate(n, 
                             return(sample(z, replace = TRUE))))))
         
         
         if(!is.null(dim(monte.c))) {
          qq <- apply(monte.c, 1,
                quantile, probs = c(0.05, 0.95),  
                na.rm = TRUE)
         } else {
          qq <- quantile (monte.c, 
                  probs = c(0.05, 0.95),  na.rm = TRUE)
         }
         
         tab <- data.frame(outer(obs, c(1:3)))
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
         
        }
        
        res <- list()
        if(nsim != 0) {
        if(test == "bootstrap") {
        res <- list("analysis" = "Getis-Ord", 
              "statistic" = classG,
              "test" = test,
              "nsim" = nsim,
              "results" = tab)
        } else {
         res <- list("analysis" = "Getis-Ord", 
               "statistic" = classG,
               "test" = test,
               "nsim" = nsim,
               "p.adjust.method" = adjust,
               "results" = tab)
        }
         
        } else {
         res <- list("analysis" = "Getis-Ord", 
               "statistic" = classG,
               "results" = obs)
        }
        
        class(res) <- "eco.gm"
        res
        
       })
