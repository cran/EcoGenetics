# Moran's I statistic with Monte-Carlo test
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015


setGeneric("eco.moran",
           function(z, con, nsim = 99, 
                    alternative = c("auto", "two.sided", "greater", "less"),
                    test = c("permutation", "bootstrap"), adjust.n = TRUE, 
           				  plotit =TRUE) {
             
             
             alternative <- match.arg(alternative)
             test <- match.arg(test)
             
             n <- length(z)
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
             
             if(adjust.n == TRUE) {
               colTRUE <- apply(wg, 1, sum)
               colTRUE[colTRUE != 0] <- 1
               n <- sum(colTRUE)
             }
             
             wg <- as.numeric(wg)
             
             z2 <- z - mean(z)
             coef.sup <- outer(z2, z2)
             diag(coef.sup) <- 0
             coef.sup <- as.numeric(coef.sup)
             
             coef.sup2 <- n * sum(wg * coef.sup) / sum (wg)
             coef.inf <- (length(z) - 1)
             coef.inf1 <- var(z) * coef.inf
             
             obs <- coef.sup2 /coef.inf1
             
             
             coefsup.mc <- n  / (coef.inf * sum (wg))
             
             repsim <- numeric()
             
             if(test == "permutation") {
               for(i in 1:nsim){
                 
                 z.boot <- sample(z2)
                 coef.sup <- outer(z.boot, z.boot)
                 diag(coef.sup) <-0
                 coef.sup <- as.numeric(coef.sup)
                 repsim[i] <- coefsup.mc * sum(wg * coef.sup) /  var(z)
                 
                 expected <- median(repsim)
                 if(alternative == "auto") {
                   alter <- expected - obs
                   if(alter > 0) {
                     alternative <- "less"
                   } else if (alter < 0) {
                     alternative <- "greater"
                   }
                 }
                 
                 if(alternative == "greater") {
                   howmany <- sum(repsim >= obs)
                   p.val <- (howmany + 1) / (nsim + 1)
                   
                 } else if(alternative == "less") {
                   howmany <- sum(repsim <= obs)
                   p.val <- (howmany + 1) / (nsim + 1)
                   
                 } else if(alternative == "two.sided") {
                   howmany <- sum(abs(repsim) >= abs(obs))
                   p.val <- (howmany + 1) / (nsim + 1)
                 }
                 
                 res <- list("analysis" = "Moran's I", 
                             "alternative"= paste("alternative hypothesis:", 
                                                  alternative), 
                             "observation" = round(obs, 4),
                             "expectation" = round(expected, 3),
                             "nsim" = nsim, "p.value" = round(p.val, 5))
                 
               }
               
             } else {
               for(i in 1:nsim){
                 z.boot  <- sample(z, replace = TRUE)
                 z2 <- z.boot -mean(z.boot)
                 coef.sup <- outer(z2, z2)
                 diag(coef.sup) <-0
                 coef.sup <- as.numeric(coef.sup)
                 repsim[i] <- coefsup.mc * sum(wg * coef.sup) / var(z.boot)
                 qq <- quantile(repsim, probs = c(0.05, 0.95),  na.rm = TRUE)
                 res <- list("analysis" = "Moran's I",
                             "observation" = round(obs, 4),
                             "nsim" = nsim,
                             "quantile" = qq)
               }
             }
             
             
             if(plotit == TRUE) {
               hist(c(repsim, obs), xlab = "Moran's I", main = "Monte Carlo test")
               abline(v = obs, col = "red")
               points(obs, 0, col = "green", pch = 15, cex = 3.6)
               text(obs, 0, "obs")
             }
             
             res
             
           })
