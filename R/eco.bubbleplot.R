#' Bubbleplot graphs
#' @name eco.bubbleplot
#' @description This function generates a bubble plot for spatial variables.
#' @param data Data frame or matrix with X-Y coordinates and variable Z.
#' @param xlabel Optional label for x axis.
#' @param ylabel Optional label for y axis.
#' @param title Optional title label.
#' @param background color of the background ("grey" or "white")-
#' @param legendlabel Optional legend label.
#' @param significant should be colored only the individuals with significant 
#' result?. This argument can be used with \code{\link{eco.lsa}} results. 
#' Default TRUE
#' @param ns color for non significant individuals, when significant = TRUE.
#' This argument can be used with \code{\link{eco.lsa}} results.
#' @param ... Additional elements to the generic.
#' 
#' @examples
#' \dontrun{
#' data(eco3)
#' 
#' # The data set eco3 has  50 points in two sites, 
#' # but they are not visible in a usual X-Y plot 
#' due to the small distance among them
#' 
#' var <- eco3[["P"]][,1]
#' plot(eco3[["XY"]], col = var)
#' x <- sample(1:100, 30)
#' y <- sample(1:100, 30)
#' 
#' # in a rankplot graph, the inter-individual distances are
#' # reduced to a single scale
#' rankeco3 <- eco.rankplot(var, eco3[["XY"]])
#' rankeco3
#' 
#' # the rankplot method support the use of ggplot2 syntax
#' rankeco3 <- rankeco3 + theme_bw() + theme(legend.position="none")
#' rankeco3
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @rdname rankplot-methods
#' 
#' @aliases eco.rankplot,genetic-method
#' 
#' @exportMethod eco.rankplot


setGeneric("eco.bubbleplot",
           function(data,
                    xlabel =  NULL,
                    ylabel=  NULL,
                    legendlabel=  NULL,
                    title =  NULL,
                    col.ns = "white",
                    significant = TRUE,
                    alpha = 0.05,
                    point.size = 1,
                    bg.size = 1,
                    ...) {
                
   
   if(is.null(xlabel)) {
     xlabel <- "Longitude"
   } 
   
   if(is.null(ylabel)) {
     ylabel <- "Latitude"
   } 
   
   if(is.null(title)) {
       title <- " "
     }
   
   if(is.null(legendlabel)) {
     legendlabel <- ""
   }
   
   if(significant) {
     colnames(data)[4] <- "p"
     data2 <- data
     data2[which(data[, 4] < alpha), ] <- NA 
   } else {
     data2 <- data
   }
   
   data[which(is.na(data[, 3])), ] <- NA
   data[, 3] <- point.size * data[, 3]
   
   bg.size <- 9*bg.size
   bg.size <- rep(bg.size, nrow(data))
   data[, 3]  <- data[, 3] * point.size
   data <- cbind(data, bg.size)
   colnames(data)[1:3] <- c("X", "Y", "Z")   
   
   bubble <-  ggplot(data = data, aes(x = X, y = Y)) +
     ggplot2::geom_point(data = data, colour = "grey", size = bg.size, 
                         alpha =0.4) +
     ggplot2::geom_point(data = data2, ggplot2::aes(size = Z),
                                   colour ="red", alpha = 0.5) +
     ggplot2::scale_size_continuous(range= c(2,15)) +
     theme_bw()+
     ggplot2::xlab(xlabel) +
     ggplot2::ylab(ylabel) +
     ggplot2::labs(title = title) +
     ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
                    axis.title = ggplot2::element_text(size = 14, 
                                                       face = "bold"), 
                    plot.title = ggplot2::element_text(size = 16,
                                                       face = "bold"))
   
   
  
  xy.out <- data[, 1:2, drop = FALSE]

  attr(bubble, "data") <- xy.out

  message(paste("plot option: significant =", significant))
  
  bubble
  
})

setMethod("eco.bubbleplot", "eco.lsa", 
          function(data, 
                   xlabel =  NULL,
                   ylabel=  NULL,
                   legendlabel=  NULL,
                   title =  NULL,
                   col.ns = "white",
                   significant = TRUE,
                   alpha = 0.05,
                   point.size = 1,
                   bg.size = 1,
                   ...) {
            
          if(data@TEST == "permutation") {
          XYZ <- cbind(data@XY, data@OUT$obs, data@OUT$p.val)
          } else {
          XYZ <- cbind(data@XY, data@OUT$obs) 
          }

          callGeneric(data = XYZ, 
                      title = data@NAMES,
                      legendlabel =  paste("  ",data@METHOD),
                      xlabel =  xlabel,
                      ylabel =   ylabel,
                      col.ns = col.ns,
                      significant =  significant,
                      alpha = alpha,
                      point.size = point.size,
                      bg.size = bg.size,
                      ...)          

          })

