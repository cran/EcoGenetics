# Leandro Roser leandroroser@ege.fcen.uba.ar
# June 17, 2015 


# Global spatial analysis

setGeneric("eco.gsa",
           function(Z, Y = NULL, con, 
                    method =c("I", "C", "CC", "JC"),
                    ncod = NULL,
                    nsim = 99, 
                    alternative = c("auto", "two.sided", 
                                    "greater", "less"),
                    adjust = "fdr",
                    row.sd = FALSE,
                    plotit =TRUE) {
             
             
             alternative <- match.arg(alternative)
             method <- match.arg(method)
             
             if(any(class(con) == "ecoweight")) {
               XY <- attr(con, "XY")
             } else {
               XY <- attr(con, "xy")
               con <- int.check.con(con)
             }
             
             Z.class <- class(Z)
             if(method != "JC") {
             if(Z.class == "matrix" | Z.class == "data.frame") {
               c.Z <- apply(Z, 2, class)
               if(any(c.Z != "integer") & any(c.Z != "numeric")) {
                 stop(paste("Non numeric data.", method, "requires numeric data"))
               }
             } else if(Z.class != "numeric" & Z.class != "integer") {
               stop(paste("Non numeric data.", method, "requires numeric data"))
             }
             } else {
               if(Z.class == "matrix" | Z.class == "data.frame") {
                 c.Z <- apply(Z, 2, class)
                 if(any(c.Z != "factor")) {
                   warning(paste("Non-factor data present. Elements in variables will be 
                                 treated as factor levels"))
                 }
               } else if(Z.class != "numeric" & Z.class != "integer" & Z.class != "factor" & Z.class != "character") {
                 stop("unknown data class")
               }
             }
               
              
             if(!is.null(Y)) {
               Y.class <- class(Y)
               if(Y.class != "numeric" & Y.class != "integer" & Y.class != "vector") {
                 stop("Y must me a numeric vector")
               }
             }
             
             
             nvar <- ncol(Z)
             
             
             ##### END OF CHECKS #####
             
             
             
             multiple <- NULL
             if(ncol(as.matrix(Z)) == 1) {
               multiple <- FALSE
             } else {
               multiple <- TRUE
               plotit <- FALSE
             }
             
             ###test#####
             
             # Funtion to estimate the stat in each iteration
             
             select_method <- function(u, ...) {
               
               if(method == "I") {
                 out <- int.moran(Z = u, plotit = plotit,  ...)
               } else if(method == "C") {
                 out <- int.geary(Z = u, plotit = plotit,  ...)
               } else if(method == "CC") {
                 out <- int.crosscor(Z = u, Y = Y, plotit = plotit, ...)
               } else if(method == "JC") {
                 out <- int.joincount(Z = u, ncod = ncod, 
                                      adjust = adjust, ...)
               }
               out
             }
             
             sel <- match(method,  c("I", "C", "CC", "JC"))
             name <- c("Moran' I", "Geary's C", 
                       "Bivariate Moran's Ixy", "Join-count")
             name <- name[sel]
             salida <- new("eco.gsa")
             
             #one test
             if(!multiple) {
               
               res <- select_method(Z, con = con, nsim = nsim,
                                    alternative = alternative) 
             } 
             
             #rearranging no JC
             if(!multiple & method != "JC") {
               
               salida@METHOD <- name
               salida@OBS <- res$obs
               salida@EXP <- res$exp
               salida@PVAL <- res$p.val
               salida@ALTER <- res$alter
               salida@NSIM <- res$nsim
             }

               
             if(multiple) {
               #multiple tests
            
               res <- list()
               
               #test
               for(i in 1:nvar) {
                 
                 res[[i]] <- select_method(Z[, i], con = con, nsim = nsim,
                                           alternative = alternative) 
                 
               }
               if(method != "JC") {
               res <- int.multitable(res)
               salida@METHOD <- name
               salida@NSIM <- nsim
               salida@ADJUST <- adjust
               salida@MULTI <- res$results
               }
             }
                 
                 #rearranging JC
                 if(method == "JC") {
                   salida@METHOD <- name
                   salida@NSIM <- nsim
                   salida@ADJUST <- adjust
                   salida@MULTI <- res$results[,1:4]
                 }
             
             salida
             })
            


