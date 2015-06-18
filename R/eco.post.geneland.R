# Leandro Roser leandroroser@ege.fcen.uba.ar
# June 17, 2015 


# Log posterior probability plot for Geneland repetitions with fixed K

setGeneric("eco.post.geneland", 
           function(niter, burnin) {
             
             
             posterior <- data.frame(c(1:niter), rep(0, niter))
             colnames(posterior) <- c("repeticion", "posterior")
             
             for(i in 1:niter) {
               path <- getwd()
               path <- paste(path, "/", i, "/", 
                             "log.posterior.density.txt",
                             sep = "")
               logmedio <- read.table(path)
               temporal <- mean(logmedio[-c(1:burnin), 1])
               posterior[i, 2] <- temporal
             }
             
             plotfun <- ggplot2::ggplot(data = posterior) +
               ggplot2::geom_line(ggplot2::aes(x = repeticion,
                                               y = posterior),
                                  directions = "hv",
                                  linetype = 2, colour = "red") + 
               ggplot2::geom_point(ggplot2::aes(x = repeticion, y = posterior)) +
               ggplot2::scale_x_discrete() + 
               ggplot2::xlab("Repetition") + 
               ggplot2:: ylab("log posterior probability")
             
             plotfun 
             
           })
