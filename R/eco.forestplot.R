##### Forestplot graphs #####

# Leandro Roser leandroroser@ege.fcen.uba.ar
# May 11, 2015

# generic function

setGeneric("eco.forestplot", 
           function(input, 
                    xlabel = NULL,
                    ylabel = NULL,
                    titlelabel = NULL,
                    legendlabel = NULL) {
             standardGeneric("eco.forestplot")
           })


# eco.forestplot,eco.lsa-method

setMethod("eco.forestplot", 
          "eco.lsa",
          function(input,
                   xlabel,
                   ylabel,
                   titlelabel,
                   legendlabel) {
            
            
            if(input@TEST != "bootstrap") {
              stop("this method is available for eco.lsa with bootstrap test")
            }
            
            datos <- input@OUT
            
            ind <- rownames(datos)
            obs <- datos$obs
            lwr <- datos$lwr
            uppr <- datos$uppr
            method <- input@METHOD
            
            data.select <- data.frame(ind, obs, lwr, uppr)
            
            
            
            if(is.null(xlabel)) {
              xlabel <- method
            } 
            
            if(is.null(ylabel)) {
              ylabel <- "Individual"
            } 
            
            if(is.null(titlelabel)) {
              titlelabel <- " "
            }
            
            if(is.null(legendlabel)) {
              legendlabel <- paste("  ", method)
            }
            
            p <- ggplot2::ggplot(data.select, ggplot2::aes(x = c(1:length(ind)), 
                                                           y = obs, 
                                                           ymin = lwr,
                                                           ymax = uppr)) + 
              ggplot2::geom_pointrange(ggplot2::aes(colour = obs), size=0.9) +
              ggplot2::geom_point(size=2.5,  shape=1) +
              ggplot2::geom_point(size=2.7,  shape=1) +
              
              #ggplot2::geom_line(ggplot2::aes(x= XVALUE, y = YVALUE), lwd = 1, alpha = 0.4) + 
              ggplot2::coord_flip() + 
              ggplot2::scale_colour_gradient(legendlabel, 
                                             low = "green", 
                                             high = "red") +
              ggplot2::ylab(xlabel) +
              ggplot2::xlab(ylabel) +
              ggplot2::labs(title = titlelabel)+
              ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
                             axis.title = ggplot2::element_text(size = 14, 
                                                                face = "bold"), 
                             plot.title = ggplot2::element_text(size = 16, 
                                                                face = "bold"))
            attr(p, "data") <- data.select
            p
            
          })


# eco.forestplot,dataframeORmatrix-method

setMethod("eco.forestplot", 
          "dataframeORmatrix",
          function(input,
                   xlabel,
                   ylabel,
                   titlelabel,
                   legendlabel) {
            
            
            datos <- input
            
            ind <- rownames(datos)
            obs <- datos[, 1]
            lwr <- datos[, 2]
            uppr <- datos[, 3]
            
            data.select <- data.frame(ind, obs, lwr, uppr)
            
            if(is.null(xlabel)) {
              xlabel <- "value"
            } 
            
            if(is.null(ylabel)) {
              ylabel <- "Individual"
            } 
            
            if(is.null(titlelabel)) {
              titlelabel <- ""
            }
            
            if(is.null(legendlabel)) {
              legendlabel <- ""
            }
            
            p <- ggplot2::ggplot(data.select, ggplot2::aes(x = c(1:length(ind)), 
                                                           y = obs, 
                                                           ymin = lwr,
                                                           ymax = uppr)) + 
              ggplot2::geom_pointrange(ggplot2::aes(colour = obs), size=0.9) +
              ggplot2::geom_point(size=2.5,  shape=1) +
              ggplot2::geom_point(size=2.7,  shape=1) +
              ggplot2::coord_flip() +  
              ggplot2::ylab(xlabel) +
              ggplot2::xlab(ylabel) +
              ggplot2::labs(title = titlelabel)+
              ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
                             axis.title = ggplot2::element_text(size = 14, 
                                                                face = "bold"), 
                             plot.title = ggplot2::element_text(size = 16, 
                                                                face = "bold"))
            
            attr(p, "data") <- data.select
            p
            
          })
