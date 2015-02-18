# Allelic frequency histogram for an ecogen genetic data frame
# Leandro Roser leandroroser@ege.fcen.uba.ar
# February 18, 2015 

setGeneric("eco.alfreq", function(eco, fact = NULL) {
  
  if(!is.null(fact)) {
    cual <- which(colnames(eco$S) == fact)
    fact.num <- as.numeric(levels(eco$S[,cual]))[eco$S[,cual]] 
    nfact <- max(fact.num)
  } else {
    eco$S$dummy <- rep(1, nrow(eco$G))
    cual <- which(colnames(eco$S) =="dummy")
    fact.num <- as.numeric(levels(as.factor(eco$S[,cual])))[eco$S[,cual]] 
    nfact <- 1
  }
  
  for(i in 1:nfact) {
    
    eco2 <- eco[which(eco$S[,cual] == i)]
    clases<- as.numeric(eco2$GENIND$loc.fac)
    tabla <- eco2@GENIND$tab
    tabla <- 2 * tabla
    frecuencia <- apply(tabla, 2, sum)
    alelos.locus <- tapply(frecuencia, clases, sum)
    for( j in 1:length(clases)) {
      temp <- clases[j]
      clases[j] <- alelos.locus[temp]
    }
    frecuencia <- frecuencia / clases
    
    #tapply(frecuencia, clases, sum) #verification
    
    frecuencia <- as.data.frame(frecuencia)
    if(nfact >1) {
      tit <-paste("Pop",levels(eco$S[, cual])[i])
    } else {
      tit <-""
    }
    
    grafico<- ggplot2::ggplot(frecuencia, ggplot2::aes(frecuencia),
                              fill = "black") + 
      ggplot2::geom_histogram(ggplot2::aes(y = ..density..)) + 
      ggplot2::geom_density(alpha = 0.5, fill = "red") +
      ggplot2::labs(title = tit) +
      ggplot2::xlab("Frequency") +
      ggplot2::ylab("Density")
    suppressMessages(print(grafico))
  }
  
})
