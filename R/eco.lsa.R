# Local spatial analysis

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015

setGeneric("eco.lsa",
           function(Z, con, method = c("G*","G", "I", "C"),
                    zerocon = NA, nsim = 99, 
                    conditional = c("auto", TRUE, FALSE),
                    test = c("permutation", "bootstrap"),
                    alternative = c("auto", "two.sided", 
                                    "greater", "less"), 
                    adjust = "fdr") {
             
             method <- match.arg(method)
             conditional <- match.arg(conditional)
             
             test <- match.arg(test)
             alternative.i <- match.arg(alternative)
             
             
             if(conditional == "auto") {
               if(method == "G*") {
                 conditional <- FALSE
               } else {
                 conditional <- TRUE
               }
             }
             
             
             
             #weight configuration
             
             zerocon2 <- match(zerocon, c(0, NA))
             if(is.null(zerocon2)) {
               stop("zerocon argument must be 0 or NA")
             }
             
             
             if(class(con) == "eco.weight") {
               XY <- con@XY
               con <- con@W
             } else {
               con <- int.check.con(con)
               if(attr(con, "xy") == NULL) {
                 stop("The weight matrix requires an attribute <xy>")
               }
               XY <- attr(con, "xy")
             }
             
             
             #general configuration. method selection
             n <- length(Z)
             counter <- 0
             cat("\n")
             
             ############## getis ####################
             
             if(method == "G*"| method == "G") {
               classG <- method
               
               Gm <- t(replicate(n, return(Z)))
               
               
               if (classG == "G") {
                 if(any(diag(con) != 0)) {
                   diag(con) <- 0
                   warning(paste("Non zero elements in the diagonal of the weight matrix", 
                                 "self-included individuals). These values are set to 0 for G"))
                 }
                 diag(Gm) <- 0
                 n2 <- n - 1
               } else if (classG == "G*"){
                 if(any(diag(con) == 0)) {
                   stop(paste("Individuals non self-included in the weight matrix",
                              "(zeros present in the diagonal). Self-included individuals",
                              "are required for G*"))
                 }
                 n2 <- n
               }
               
               #specific function for G* or G
               select.method <- function(Gm) {
                 
                 G <- con * Gm
                 G <- apply(G, 1, sum)
                 Gm2 <- Gm ^ 2
                 
                 meanG <- apply(Gm, 1, sum) / n2
                 desvsqG <- apply(Gm2, 1, sum) / n2
                 desvsqG <- desvsqG - meanG ^ 2
                 W <- apply(con, 1, sum)
                 S1 <- apply(con ^ 2, 1, sum)
                 denom <- n2 * S1 - W ^ 2
                 denom <- sqrt(desvsqG * denom / (n2-1))
                 numer <- (G - W * meanG)
                 G <- numer / denom
                 G
               }
               
               
               ############## moran ####################
               
             } else if(method == "I") {
               
               Z2 <- Z - mean(Z)
               m2 <-  sum(Z2 ^ 2) / n 
               Gm <- t(replicate(length(Z), return(Z2)))
               select.method <- function(Gm) {
                 
                 coef.sup <- apply(Gm * con, 1, sum)
                 out <- Z2 * coef.sup / m2
                 out
               }
               
               ############## geary ####################
               
             } else if(method == "C") {
               
               Z2 <- Z - mean(Z)
               m2 <-  sum(Z2 ^ 2) / n 
               Gm <- as.matrix(dist(Z, upper = T))
               
               select.method <- function(Gm) {
                 
                 num <- Gm ^ 2
                 num <- con * Gm 
                 num <- apply(num, 1, sum)
                 out <- num / m2
                 out
                 
               }
               
             }
             
             obs <- select.method(Gm)
             
             
             ################ TESTING #################################
             if(test == "permutation") {
               replace <- FALSE
             } else {
               replace <- TRUE
             }
             
             monte.c <- matrix(0, nrow(Gm), nsim)
             #fixed pivot
             if(conditional) {
               for(k in 1:nsim) {
                 Gm.test <- Gm
                 samp <- 1:nrow(Gm)
                 for(i in samp) {
                   order.Z <- sample(samp[-i], replace = replace)
                   Gm.test[i, -i]<- Gm.test[i, order.Z]
                 }
                 monte.c[, k] <- select.method(Gm.test)
                 
                 cat(paste("\r", "simulations...computed",
                           ceiling(100 * counter / nsim), "%"))
                 counter <- counter + 1
                 
               }
               #free sampling
             } else {
               for(i in 1:nsim) {
                 monte.c[, i] <- select.method(Gm[,sample(ncol(Gm), 
                                                          replace = replace)])
                 
                 cat(paste("\r", "simulations...computed",
                           ceiling(100 * counter / nsim), "%"))
                 counter <- counter + 1
               }
             }
             
             
             monte.c <- t(monte.c)
             
             tab <- int.random.test(repsim = monte.c, 
                                    obs = obs,
                                    nsim = nsim,
                                    test = test, 
                                    alternative = alternative,
                                    adjust = adjust)
             
             
             ################ END OF  TESTING #################################
             
             connect <- apply(con, 1, sum)
             if(is.na(zerocon)) {
               tab[which(connect == 0), ] <- rep(NA, 3)
             } else {
               tab[which(connect == 0), ] <- rep(0, 3)
             }
             
             
             if(!is.null(rownames(Z))) {
               rownames(tab) <- rownames(Z)
             }
             
             
             sel <- match(method,  c("G", "G*", "I", "C"))
             name <- c("Getis Ord's G", "Getis Ord's G*", 
                       "local Moran's I", "local Geary's C")
             method <- name[sel]
             
             
             res <- new("eco.lsa")
             
             
             if(test == "bootstrap") {
               
               res@METHOD <- method
               res@TEST <- test
               res@NSIM <- nsim
               res@COND <- conditional
               res@OUT <- tab
               res@XY <- data.frame(XY)
               
             } else {
               
               res@OUT <- tab
               res@METHOD <- method
               res@TEST <- test
               res@NSIM <- nsim
               res@COND <- conditional
               res@PADJ <- adjust
               res@XY <- data.frame(XY)
               
             }
             
             
             cat("\n", "done!", "\n\n")
             res
             
           })
