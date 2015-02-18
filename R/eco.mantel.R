# Mantel test
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015


setGeneric("eco.mantel", 
           function(distm1, distm2, nsim = 99, 
                    alternative = c("auto", "two.sided", "less", 
                                    "greater"), plotit = TRUE) {
             
             if(class(distm1) != "dist" | class(distm1) != "dist") {
               stop("non dist class matrix in the arguments. 
                    The imput data must be of class dist")
             }
             
             alternative <- match.arg(alternative)
             
             m1 <- as.numeric(distm1)
             m1 <-  (m1 - mean(m1)) / sd(m1)
             m2 <- as.numeric(distm2)
             m2 <-  (m2 - mean(m2)) / sd(m2)
             
             d <- length(m1) - 1
             
             obs <- sum(m1 * m2) /d
             
             repsim <- replicate(nsim, sum(sample(m1) * m2)) / d
             
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
             
             if(plotit == TRUE) {
               hist(c(repsim, obs), xlab = "Mantel statistic", 
                    main = "Monte Carlo test")
               abline(v = obs, col = "red")
               points(obs, 0, col = "green", pch = 15, cex = 3.6)
               text(obs, 0, "obs")
             }
             
             list("analysis" = "Mantel test",
                  "alternative"= paste("alternative hypothesis:", alternative),
                  "observation" = round(cor(m1, m2), 4), 
                  "expectation" = round(expected, 4), 
                  "nsim" = nsim,  "p.value" = p.val)
             
             })
