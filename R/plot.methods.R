######################
#### PLOT METHODS ####
######################

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015 

# Plot method for correlograms and variograms

setMethod("plot", "eco.correlog", function(x, var = NULL, 
                                           xlabel = NULL, 
                                           ylabel = NULL, 
                                           title = NULL, 
                                           legend = TRUE,
                                           background = c("grey", "white"),
                                           errorbar = FALSE,
                                           intervals = TRUE,
                                           significant = FALSE,
                                           xlim = NULL,
                                           ylim = NULL,
                                           axis.size = 14,
                                           title.size = 16) {
  
  
  if(length(x@OUT) == 1) {
    var2 <- 1
  } else if(is.null(var)) {
    stop("var argument not found")
  } else {
    var2 <- which(names(x@OUT) %in% var)
  }
  
  randtest <- x@TEST
  if(length(randtest) == 0) { 
    randtest <- "notest"
  }
  
  method <- (x@METHOD)[1]
  
  method2 <- pmatch(method, c("Mantel test",
                              "Partial Mantel test", 
                              "Moran's I", 
                              "Geary's, C",
                              "bivariate Moran's Ixy",
                              "empirical variogram"
  ))
  if(length(method2) == 0) {
    stop("invalid input to plot")
  }
  
  
  datos <- as.data.frame(x@OUT[[var2]])
  
  
  ########## title,  x and y labels ############
  
  if(is.null(title)) {
    title <- names(x@OUT[var2])
  }
  
  if(is.null(xlabel)) {
    xlabel <- "Great circle distance"
  }
  
  if(is.null(ylabel)) {
    ylabel <- method
  }
  
  ########## x and y axes ############
  
  if(!is.null(xlim)) {
    xlim <- ggplot2::scale_x_continuous(limits = xlim)
  }
  
  if(!is.null(ylim)) {
    ylim <- ggplot2::scale_y_continuous(limits = ylim)
  }
  
  localenv <- environment()
  
  
  ########## basic plot ############
  
  z <- ggplot2::ggplot(datos, environment = localenv) + 
    ggplot2::geom_line(ggplot2::aes(x = d.mean, y = obs)) + 
    ggplot2::xlab(xlabel) + 
    ggplot2::ylab(ylabel) + 
    ggplot2::labs(title = title) 
  
  
  ######## background ########
  
  theme <- match.arg(background)
  if(theme == "grey") {
    
    themecol <-  ggplot2::theme_grey()
  } else {
    themecol <- ggplot2::theme_bw()
  }
  
  z <- z + themecol
  
  ######## permutation and bootstrap cases ########
  
  if(randtest == "permutation") {
    #labeling S and NS points
    pval2 <- as.numeric(datos$pval < 0.05)
    if(sum(pval2) == 0) {
      labelp <- "NS"
    } else if (sum(pval2) == length(pval2)) {
      labelp <- "P<0.05"
    } else {
      labelp <- c("P < 0.05", "NS")
    }
    pval2[pval2 == 0] <- "#F8766D"
    pval2[pval2 == 1] <- "#00B0F6"
    datos$pval2 <- pval2
    puntos <- ggplot2::geom_point(ggplot2::aes(x = d.mean, y = obs, colour = pval2), 
                                  size = 5)
    z <- z + puntos
  }
  
  if(randtest == "bootstrap")  {
    intervalo <- ggplot2::geom_ribbon(ggplot2::aes(x = d.mean, ymax = uppr, 
                                                   ymin = lwr), 
                                      fill = "blue",
                                      alpha = 0.2)
    
    puntos <- ggplot2::geom_point(ggplot2::aes(x = d.mean, y = obs, colour = "blue"), 
                                  size = 5)
    z <- z + puntos + intervalo
    legend <- FALSE
  }
  
  if(randtest == "notest" & method == "empirical variogram") {
    z <- z + ggplot2::ylab("Semivariance") + 
      ggplot2::geom_point(ggplot2::aes(x = d.mean,y = obs), 
                          colour = "#F8766D", size =5)
    
    legend <- FALSE
  }
  
  ########## error bar ############
  
  if(errorbar) {
    datos$sd <- 1.96 * datos$sd.jack
    datos$min.sd <- datos$mean.jack - datos$sd
    datos$max.sd <- datos$mean.jack + datos$sd
    z <- z + ggplot2::geom_errorbar(data =datos,
                                    ggplot2::aes(x = d.mean, 
                                                 ymax = max.sd , 
                                                 ymin= min.sd))
  }
  
  ########## legend (permutation case) ############
  
  if(legend) {
    z <- z +  ggplot2::scale_colour_discrete(name  ="P value",
                                             labels= labelp) + 
      ggplot2::theme(axis.text = ggplot2::element_text(size = axis.size), 
                     axis.title = ggplot2::element_text(size = title.size, face = "bold"), 
                     legend.position = "right",  
                     plot.title = ggplot2::element_text(size = title.size, face = "bold"))
  }
  
  if(!legend) {
    z <- z + theme(legend.position="none", axis.text = ggplot2::element_text(size = axis.size), 
                   axis.title = ggplot2::element_text(size = title.size, face = "bold"))
  }
  
  ########## xlim, ylim ############
  
  if(!is.null(xlim)) {
    z <- z + xlim
  }
  if(!is.null(ylim)) {
    z <- z + ylim
  }
  
  ######################
  
  return(z)
  
})

#############################################################

# plot eco.lsa

setMethod("plot", "eco.lsa", function(x, ...) {
  if(x@TEST == "permutation") {
    eco.rankplot(x, ...)
  } else if(x@TEST == "bootstrap") {
    eco.forestplot(x, ...)
  }
})
