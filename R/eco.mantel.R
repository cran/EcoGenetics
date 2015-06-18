# Leandro Roser leandroroser@ege.fcen.uba.ar
# June 17, 2015 


# Mantel and partial Mantel tests

setGeneric("eco.mantel", 
					 function(d1, d2, dc = NULL, 
					          method = c("pearson", "spearman", "kendall"),
					          nsim = 99,  
					 				 alternative = c("auto", "two.sided", "less", 
					 				 								"greater"), 
					 				 ...) {
					 	
					 	alternative <- match.arg(alternative)
					 	method <- match.arg(method)
					 
					 	control <- c(class(d1), 
					 										class(d2), 
					 										class(dc))%in% "dist"
					 	sumcontrol<- sum(control)
					 	
					 	if(sumcontrol != 3 & !is.null(dc)) { 
					 		nondist <- which(!(control))
					 		dcont <- c("d1", "d2", "dc")
					 		dcont <- dcont[nondist]
					 		dcont <- paste(dcont, collapse=", ")
					 		stop(paste("non dist object/s found in the arguments. 
					 				 The imput data must be of class dist:", dcont))
					 	}
					 
					 	res <- int.mantel(d1 = d1, d2 = d2, dc = dc,
					 	                  method = method, nsim = nsim,
					 	                  test = "permutation", 
					 	                  alternative = alternative, 
					 	                  plotit = TRUE)
					 	
					 	
					 	if(is.null(dc)) { 
					 	  method.mantel <- "Mantel test"
					 	} else {
					 	  method.mantel <- "Partial Mantel test"
					 	}
					 	
					 	salida <- new("eco.gsa")
					 	salida@METHOD <- c(method.mantel, method)
					 	salida@OBS <- round(res$obs, 4)
					 	salida@EXP <- round(res$exp, 4)
					 	salida@PVAL <- res$p.val
					 	salida@ALTER <- res$alter
					 	salida@NSIM <- nsim

					 	salida
					 	
					 	})
