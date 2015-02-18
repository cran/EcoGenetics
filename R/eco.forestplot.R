# Plot for correlograms of Getis- Ord's G or local-Moran's I 
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015

setGeneric("eco.forestplot", function(input, var) {
 standardGeneric("eco.forestplot")
})


setMethod("eco.forestplot", 
     c("eco.multiboot", "character"), 
     function(input, var) {
 
 if(length(input@OUT) == 1) {
  var2 <- 1
 } else if(is.null(var)) {
  stop("a variable (var) argument must be selected")
 } else {
  var2 <- which(names(input@OUT) %in% var)
 }
 
 data <- input@OUT[[var2]]
 
 cat("the following interval distances are available: ", "\n")
 mean_distance <- colnames(data$observed)
 select <- 1:length(mean_distance)
 print(data.frame(select, mean_distance))
 
 sel = -9
 while(sum(select %in% sel) == 0) {
 cat("\n", "please select one interval of the following: ", "\n",  select)
 cat("\n\n", "your choice: ")
 sel <- scan(n = 1)
 if(sum(select %in% sel) == 0) {
  cat("incorrect choice", "\n")
 }
 }
 cat(sel)
 
 
 ind <- rownames(data$observed)
 obs <- data$observed[, sel]
 lwr <- data$lwr[, sel]
 uppr <- data$uppr[, sel]
 
 data.select <- data.frame(ind, obs, lwr, uppr)
 
 p <- ggplot2::ggplot(data.select, ggplot2::aes(x = c(1:length(ind)), 
                y = obs, 
                ymin = lwr,
                ymax = uppr)) + 
  ggplot2::geom_pointrange(ggplot2::aes(colour = obs)) +
  coord_flip() + 
  ggplot2::geom_hline(ggplot2::aes(x = 0), lty = 2)+ 
  ggplot2::scale_colour_gradient("Getis-Ord's", 
                  low = "green", 
                  high = "red")+
  ggplot2::xlab("Individual") +
  ggplot2::ylab("Getis-Ord's statistic") +
  ggplot2::labs(title = paste(var, 
                "(mean lag distance =", 
                colnames(data$observed)[sel], "m)"))+
 ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
         axis.title = ggplot2::element_text(size = 14, 
                           face = "bold"), 
         plot.title = ggplot2::element_text(size = 16, 
                           face = "bold"))
 
 print(p)
 
 return(data.select)
 
 
})



setMethod("eco.forestplot", 
     c("eco.gm"),
     function(input) {
 
 if(input$test != "bootstrap") {
  stop("the input has not a bootstrap intervals (test = permutation)")
 }
 
 datos <- input$results
 
 ind <- rownames(datos)
 obs <- datos$obs
 lwr <- datos$lwr
 uppr <- datos$uppr
 method <- input$analysis
 
 data.select <- data.frame(ind, obs, lwr, uppr)
 
 p <- ggplot2::ggplot(data.select, ggplot2::aes(x = c(1:length(ind)), 
                         y = obs, 
                         ymin = lwr,
                         ymax = uppr)) + 
  ggplot2::geom_pointrange(ggplot2::aes(colour = obs)) +
  coord_flip() + 
  ggplot2::geom_hline(ggplot2::aes(x = 0), lty = 2)+ 
  ggplot2::scale_colour_gradient(method, 
                  low = "green", 
                  high = "red")+
  ggplot2::xlab("Individual") +
  ggplot2::ylab(method) +
  ggplot2::labs(title = paste("mean lag distance (m)"))+
  ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
          axis.title = ggplot2::element_text(size = 14, 
                            face = "bold"), 
          plot.title = ggplot2::element_text(size = 16, 
                            face = "bold"))
 
 print(p)
 
 return(data.select)
 
 
})
