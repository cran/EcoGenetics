# Leandro Roser leandroroser@ege.fcen.uba.ar
# June 17, 2015 


# Rankplot graphs

setGeneric("eco.rankplot", function(input, 
																		XY, 
																		xlabel = NULL,
																		ylabel = NULL,
																		title = NULL,
																		legendlabel = NULL,
																		background = c("grey", "white"),
																		...) {
	standardGeneric("eco.rankplot")
})


# eco.rankplot,eco.lsa-method

setMethod("eco.rankplot", 
					c("eco.lsa", 
						"missing", 
						"missing"),
					function(input, 
									 XY, 
									 xlabel,
									 ylabel,
									 title,
									 legendlabel,
									 background = c("grey", "white"),
									 significant = TRUE,
									 ns = NULL) {
	
					  
 if((input@TEST)[1] != "permutation") {
	 stop("this method is available for eco.lsa with permutation test")
 }
					  
	theme <- match.arg(background)
	if(!is.null(ns)) {
		col.ns <- ns
	}
	if(theme == "grey") {
		themecol <-  ggplot2::theme_grey()
		if(is.null(ns)) {
		col.ns <- "moccasin"
		}
	} else {
		themecol <- ggplot2::theme_bw()
		if(is.null(ns)) {
		col.ns <- "white"
		}
	}
	
	XY <- input@XY

	x <- rank(XY[, 1])
	y <- rank(XY[, 2])
	method <- input@METHOD
		
	z <- input@OUT$obs
	p <- input@OUT$p.val
	
	if(significant == TRUE) {
		alpha <- as.numeric(p < 0.05)
		z <- alpha * z
	}
	
	data <- data.frame(x, y, z)
	data[which(is.na(z)), ] <- NA
	
	if(significant == TRUE) {
		data2 <- data[-which(p < 0.05), ]
		p.points <- ggplot2::geom_point(data = data2, 
																		colour =col.ns, 
																		size = 3) 
	}
	
	
	if(is.null(xlabel)) {
		xlabel <- "X rank"
	} 
	
	if(is.null(ylabel)) {
		ylabel <- "Y rank"
	} 
	
	if(is.null(title)) {
		title <- " "
	}
	
	if(is.null(legendlabel)) {
		legendlabel <- paste(" ", method)
	}
		
		rankplot <- ggplot2::ggplot(data, ggplot2::aes(x = x,y = y)) +
		  ggplot2::geom_point(colour="grey50", size = 4.5)+
			ggplot2::geom_point(ggplot2::aes(colour = z), 
													size=3)+ 
		  ggplot2::scale_color_gradient2(legendlabel, 
		  															 high= scales::muted("red"),
		                                 low = scales::muted("blue"))+
			themecol +
      ggplot2::xlab(xlabel) +
			ggplot2::ylab(ylabel) +
		  ggplot2::labs(title = title) +
	  	ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
									 axis.title = ggplot2::element_text(size = 14, 
									 																	 face = "bold"), 
									 plot.title = ggplot2::element_text(size = 16,
									 																	 face = "bold")) 
		

	xy.out <- data.frame(x, y)
	rownames(xy.out) <- rownames(XY)
	colnames(xy.out) <- c("X rank", "X rank")
	
	attr(rankplot, "data") <- xy.out
	if(significant == TRUE) {
		rankplot <- rankplot + p.points
	}
	message(paste("plot option: significant =", significant))
	rankplot
})


# eco.rankplot,numeric-method

#plot with a numeric variable for colours vs XY

setMethod("eco.rankplot", 
					c("numeric",
						"dataframeORmatrix", 
						"missing"),
					function(input, 
									 XY, 
									 xlabel,
									 ylabel,
									 title,
									 legendlabel,
									 background = c("grey", "white")) {
						
						theme <- match.arg(background)
						if(theme == "grey") {
							themecol <-  ggplot2::theme_grey()
						} else {
							themecol <- ggplot2::theme_bw()
						}
						
						x <- rank(XY[, 1])
						y <- rank(XY[, 2])
						
						z <- input
						
						data <- data.frame(x, y, z)

						
						if(is.null(xlabel)) {
							xlabel <- "X rank"
						} 
						
						if(is.null(ylabel)) {
							ylabel <- "Y rank"
						} 
						
						if(is.null(title)) {
							title <- " "
						}
						
						if(is.null(legendlabel)) {
							legendlabel <- "Z"
						}
						
						rankplot <- ggplot2::ggplot(data, ggplot2::aes(x = x,y = y)) +
							ggplot2::geom_point(colour="grey50", size = 4.5)+
							ggplot2::geom_point(ggplot2::aes(colour = z), 
																	size=3)+ 
							ggplot2::scale_color_gradient2(legendlabel, 
																						 high= scales::muted("red"),
																						 low = scales::muted("blue"))+
							themecol + 
							ggplot2::xlab(xlabel) +
							ggplot2::ylab(ylabel) +
							ggplot2::labs(title = title) +
							ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
														 axis.title = ggplot2::element_text(size = 14, 
														 																	 face = "bold"), 
														 plot.title = ggplot2::element_text(size = 16,
														 																	 face = "bold")) 
						
						
						
						xy.out <- data.frame(x, y)
						rownames(xy.out) <- rownames(XY)
						colnames(xy.out) <- c("X rank", "X rank")
						
						attr(rankplot, "data") <- xy.out
						rankplot
						
					})


# rankplot-methods

#plot with a factor for colours vs XY

setMethod("eco.rankplot", 
					c("factor",
						"dataframeORmatrix", 
						"missing"),
					function(input, 
									 XY, 
									 xlabel,
									 ylabel,
									 title,
									 legendlabel) {
						
						theme <- match.arg(background)
						if(theme == "grey") {
							themecol <-  ggplot2::theme_grey()
						} else {
							themecol <- ggplot2::theme_bw()
						}
						
						x <- rank(XY[, 1])
						y <- rank(XY[, 2])
						
						z <- input
						
						data <- data.frame(x, y, z)
						
						
						if(is.null(xlabel)) {
							xlabel <- "X rank"
						} 
						
						if(is.null(ylabel)) {
							ylabel <- "Y rank"
						} 
						
						if(is.null(title)) {
							title <- " "
						}
						
						if(is.null(legendlabel)) {
							legendlabel <- "Z"
						}
						
						rankplot <- ggplot2::ggplot(data, ggplot2::aes(x = x,y = y)) +
							ggplot2::geom_point(colour="grey50", size = 4.5)+
							ggplot2::geom_point(ggplot2::aes(colour = z), 
																	size=3)+ 
							ggplot2::scale_color_discrete(legendlabel)+
							themecol +
							ggplot2::xlab(xlabel) +
							ggplot2::ylab(ylabel) +
							ggplot2::labs(title = title) +
							ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
														 axis.title = ggplot2::element_text(size = 14, 
														 																	 face = "bold"), 
														 plot.title = ggplot2::element_text(size = 16,
														 																	 face = "bold")) 
						
						
						
						
						xy.out <- data.frame(x, y)
						rownames(xy.out) <- rownames(XY)
						colnames(xy.out) <- c("X rank", "X rank")
						
						attr(rankplot, "data") <- xy.out
						rankplot
						
					})


