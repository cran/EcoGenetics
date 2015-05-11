# Detrending spatial data with polynomial interpolation

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015

setGeneric("eco.detrend", 
           function(Z, XY, degree, center = TRUE, 
                    scale = FALSE, raw = FALSE, 
                    latlon = FALSE) {
             
             if(latlon) {
               XY <- SoDA::geoXY(XY[,2], XY[,1], unit=1)
             }
             
             XY.scale <- scale(XY, center = center, scale = scale)
             
             polinomio <- data.frame(poly(as.matrix(XY.scale), degree = degree, raw = raw))
             
             
             capture.output(surf.ls <- eco.lmtree(Z, polinomio, analysis = "mlm"))
             
             
             residuos <- surf.ls@RESIDUALS
             formulas  <- lapply(1:length(surf.ls@MLM), function(i) surf.ls@MLM[[i]]$call)
             names(formulas) <- colnames(Z)
             
             res <- new("eco.detrend")
             res@POLY.DEG <- degree
             res@RES <- residuos
             res@XY <- XY
             res@MODEL <- formulas
             res@ANALYSIS <- surf.ls
             
             res
             
           })
