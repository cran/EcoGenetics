# Plot for Getis- Ord's G or local-Moran's I analysis
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015

setGeneric("eco.rankplot", function(input, XY, var = NULL) {
 standardGeneric("eco.rankplot")
})

#' @rdname rankplot-methods
#' @aliases eco.rankplot,eco.multiboot-method
#' @exportMethod eco.rankplot

setMethod("eco.rankplot", 
     c("eco.multiboot",
      "dataframeORmatrix",
      "character"),
     function(input, XY, var) {

 
 if(length(input@OUT) == 1) {
  var2 <- 1
 } else if(is.null(var)) {
  stop("a variable (var) argument must be selected")
 } else {
 var2 <- which(names(input@OUT) %in% var)
 }
 
 x <- rank(XY[, 1])
 y <- rank(XY[, 2])
 
  seqq <- seq(input@INTERVAL, input@MAX, input@INTERVAL)
 
 for(u in seq(along = seqq)) {
 
z <- input@OUT[[var2]]$observed[, u]
data <- data.frame(x, y, z)
 
rankplot <- ggplot2::ggplot(data, ggplot2::aes(x = x,y = y)) +
 ggplot2::geom_point(colour="grey50", size = 4.5)+
 ggplot2::geom_point(ggplot2::aes(colour = z), size=3.5)+
 ggplot2::scale_color_gradient2("Getis-Ord's", 
                low = "blue",
                high = "red")+
 ggplot2::xlab("X rank") +
 ggplot2::ylab("Y rank") +
 ggplot2::labs(title = paste(var, 
               "(mean lag distance =", 
               colnames(input@OUT[[var2]]$observed)[u], "m)")) +
ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
         axis.title = ggplot2::element_text(size = 14, 
                           face = "bold"), 
         plot.title = ggplot2::element_text(size = 16,
                           face = "bold")) 

print(rankplot)

}

xy.out <- data.frame(x, y)
rownames(xy.out) <- rownames(XY)
colnames(xy.out) <- c("X rank", "X rank")

return(xy.out)

})


#' @rdname rankplot-methods
#' @aliases eco.rankplot,eco.gm-method
#' @exportMethod eco.rankplot

setMethod("eco.rankplot", 
     c("eco.gm", 
      "dataframeORmatrix", 
      "missing"),
     function(input, XY, var) {
 
 
 x <- rank(XY[, 1])
 y <- rank(XY[, 2])
 analysis <- input$analysis
  
  z <- input$results$obs
  data <- data.frame(x, y, z)
  
  rankplot <- ggplot2::ggplot(data, ggplot2::aes(x = x,y = y)) +
    ggplot2::geom_point(colour="grey50", size = 4.5)+
   ggplot2::geom_point(ggplot2::aes(colour = z), size=3.5)+ 
   ggplot2::scale_color_gradient2(analysis, 
                   low = "blue",
                   high = "red")+
   ggplot2::xlab("X rank") +
   ggplot2::ylab("Y rank") +
   ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
           plot.title = ggplot2::element_text(size = 16,
                             face = "bold")) 
  
  print(rankplot)
  
 
 xy.out <- data.frame(x, y)
 rownames(xy.out) <- rownames(XY)
 colnames(xy.out) <- c("X rank", "X rank")
 
 return(xy.out)
 
})

