# Leandro Roser leandroroser@ege.fcen.uba.ar
# June 17, 2015 


# Converting a diploid ecogen genetic data frame into a gstudio object

setGeneric("eco.2gstudio", 
           function(eco, type = "separated", ...) {
             
             
             if(type == "separated") {
               dat <- adegenet::df2genind(eco$G)
               dat <- adegenet::genind2df(dat, sep = ":")
               for(i in 1:ncol(dat)) {	
                 dat[, i] = gstudio::locus(dat[, i], type = "separated")
               }
             } else {
               dat<-eco$G
               for(i in 1:ncol(dat)) {
                 dat[, i] = gstudio::locus(dat[, i], type = type)
               }
             }
             
             dat
           })
