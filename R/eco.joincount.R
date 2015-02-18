# Join-Count statistic with Monte-Carlo test
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015


setGeneric("eco.joincount",
           function(z, con, nsim = 99, 
                    test = c("permutation", "bootstrap"),
                    alternative = c("auto", "two.sided", 
                                    "greater", "less"), 
                    adjust = c("fdr", "holm", "hochberg", 
                      "hommel", "bonferroni", "BH",
                      "BY", "none"), 
                    ncod = NULL) {
             
             alternative.i <- match.arg(alternative)
             adjust <- match.arg(adjust)
             test <- match.arg(test)
             
             
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
             
             
             outmat <- outer(z, z, FUN = "paste", sep = "")
             if(is.null(ncod)) {
               stop("a ncod parameter was not given") 
             } else {
               outmat <- (as.matrix(eco.sort(outmat, ncod)))
             }
             outmat <- as.vector(outmat)
             outmat <- as.factor(outmat)
             
             obs <- numeric()
             temp <- list()
             
             for(i in seq(along = levels(outmat))) {
               temp[[i]] <- as.numeric(outmat == levels(outmat)[i])
             }
             
             obs <- sapply(temp, function(x) sum(x*wg) / 2) 
             
             monte.c <- list()
             if(test == "permutation") {
               monte.c <- (replicate(nsim, 
                                     sapply(temp,
                                            function(x) sum(sample(x) * wg) / 2)))
               if(!is.null(dim(monte.c))) {
                 monte.c <- t(monte.c)
               } else {
                 monte.c <- as.matrix(monte.c)
               }
               
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
               tab[, 1] <- round(obs, 3)
               tab[, 2] <- round(p.val, 4)
               tab[, 3] <- round(expected, 3)
               tab[, 4] <- altern
               
               
               colnames(tab) <- c("jc.stat", "P", "exp", "alter")
               rownames(tab) <- levels(outmat)
               
               
             }  else {
              monte.c <- replicate(nsim, 
                                   sapply(temp, 
                                   function(x) sum(sample(x, 
                                                         replace = TRUE) * wg) / 2))
               
               if(!is.null(dim(monte.c))) {
                 qq <- apply(monte.c, 1, quantile, probs = c(0.05, 0.95),  
                             na.rm = TRUE)
               } else {
                 qq <- quantile (monte.c, probs = c(0.05, 0.95),  na.rm = TRUE)
               }
               
               tab <- outer(obs, c(1:3))
               tab[, 1] <- round(obs, 4)
               tab[, 2:3] <- t(qq)
               colnames(tab) <- c("joincount.stat", "lwr", "uppr")
               rownames(tab) <- levels(outmat)
               
             }
             
             if(nsim != 0) {
             if(test == "bootstrap") {
             res <- list("analysis" = "join count test", 
                   "test" = test,
                         "nsim" = nsim,
                   "p.adjust.method" = adjust,
                         "results" = tab)
             } else {
              res <- list("analysis" = "join count test", 
                    "test" = test,
                    "nsim" = nsim,
                    "p.adjust.method" = adjust,
                    "results" = tab)
             }
             } else {
              res <- list("analysis" = "join count test", 
                    "results" = obs)
             }
             
             res
             
           })


