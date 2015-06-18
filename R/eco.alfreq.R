# Leandro Roser leandroroser@ege.fcen.uba.ar
# June 17, 2015 


# Allelic frequency histograms for an ecogen genetic data frame

setGeneric("eco.alfreq", function(eco, grp = NULL) {
	
	if(!is.null(grp)) {
		cual <- which(colnames(eco$S) == grp)
		grp.num <- as.numeric(levels(eco$S[,cual]))[eco$S[,cual]] 
		nfact <- max(grp.num)
	} else {
		eco$S$dummy <- rep(1, nrow(eco$G))
		cual <- which(colnames(eco$S) =="dummy")
		grp.num <- as.numeric(levels(as.factor(eco$S[,cual])))[eco$S[,cual]] 
		nfact <- 1
	}
	
  out <- list()
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
			ggplot2::xlab("Frequency class") +
			ggplot2::ylab("Density")
		out[[i]] <- grafico
	}
	out
})

