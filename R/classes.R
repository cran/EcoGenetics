# Leandro Roser leandroroser@ege.fcen.uba.ar
# June 17, 2015 


# eco.correlog-class

setClass("eco.correlog", 
         
         representation(OUT = "list",
                        IN = "list",
                        BREAKS = "numeric", 
                        CARDINAL = "numeric", 
                        NAMES = "character",
                        METHOD = "character", 
                        DISTMETHOD = "character",
                        TEST = "character",
                        NSIM = "numeric",
                        PADJUST = "character")
)

#################################################

# eco.variogram class

setClass("eco.variogram", contains = "data.frame")

#################################################

# eco.gsa class

setClass("eco.gsa",  
         representation(METHOD = "character",
                        OBS = "numeric",
                        EXP = "numeric",
                        PVAL = "numeric",
                        ADJUST = "character",
                        ALTER = "character",
                        NSIM ="numeric",
                        MULTI = "list")
         
)

#################################################

# eco.lsa class

setClass("eco.lsa",  
         representation(OUT = "list",
                        METHOD = "character",
                        TEST = "character",
                        NSIM ="numeric",
                        PADJ = "character",
                        COND = "logical",
                        XY = "data.frame")
)


#################################################

# eco.weight class

setClass("eco.weight",  
         
         representation(W = "matrix",
                        XY = "data.frame",
                        METHOD = "character",
                        PAR = "character",
                        PAR.VAL = "numeric",
                        ROW.SD = "logical",
                        SELF ="logical",
                        NONZERO = "numeric",
                        NONZEROIND = "numeric",
                        AVG = "numeric")
)

#################################################

# eco.lagweight class
 
setClass("eco.lagweight",  
         
         representation(W = "list",
                        XY = "data.frame",
                        PAR = "character",
                        PAR.VAL = "numeric",
                        ROW.SD = "logical",
                        SELF ="logical",
                        CUMMUL = "logical",
                        MEAN = "numeric",
                        LOGMEAN = "numeric",
                        CARDINAL = "numeric",
                        BREAKS = "numeric",
                        METHOD = "character")
)


# eco.mlm-class

setClass("eco.mlm",
         
         representation( MLM = "list",
                         SUMMARY.MLM = "list", 
                         ANOVA.MLM = "list",
                         PREDICTED = "data.frame",
                         RESIDUALS = "data.frame",
                         df1 = "dataframeORmatrix",
                         df2 = "dataframeORmatrix")
)

# eco.mctree-class

setClass("eco.mctree",
         
         representation( TREES = "list",
                         CLASSPREDICT = "list", 
                         FREQUENCIES = "list",
                         PREDICTED = "data.frame",
                         RESIDUALS = "data.frame",
                         df1 = "dataframeORmatrix",
                         df2 = "dataframeORmatrix") 
)


# eco.detrend class

setClass("eco.detrend",  
         
         representation(POLY.DEG = "numeric",
                        RES = "data.frame",
                        XY = "data.frame",
                        MODEL = "list",
                        ANALYSIS ="eco.mlm")
                       
)


###################################

setClass("int.multiplot")
setClass("eco.IBD", representation(SP = "list"), contains = "eco.correlog")


