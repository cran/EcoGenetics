# Leandro Roser leandroroser@ege.fcen.uba.ar
# June 17, 2015 


# Fitting Multiple Linear Regression models by stepwise AIC selection and
# Multiple Classification and Regression Trees via party

setGeneric("eco.lmtree", 
           function(df1, df2, 
                    analysis = c("mlm", "mctree"), mod.class = "+", 
                    fact = NULL, ...) 	{
             
             
             analysis <- match.arg(analysis)
             df1 <- data.frame(df1)
             df2 <- data.frame(df2)
             
             if((mod.class == "+") || (mod.class == "*")) {
               indep <- paste(colnames(df2), collapse = mod.class, sep = "")
             } else {
               indep <- mod.class
             }
             
             data <- data.frame(df1, df2)
             
             if(analysis == "mlm") {
               
               mlm.mod <- new("eco.mlm")
               mod <- results <- anovas <- list()
               
               capture.output(for(i in 1:ncol(df1)) {
                 dep <- paste(colnames(df1)[i], "~", sep = "")
                 smod <- paste(dep, indep, sep = "")
                 smod <- as.formula(smod)
                 mod[[i]] <- lm(smod, data = data, ...)
                 mod[[i]] <- step(mod[[i]], scope = list(colnames(df1)[i] ~ 1,
                                                         upper = mod[[i]]))
               })
               
               for(i in 1:ncol(df1)) {
                 results[[i]] <- summary(mod[[i]])
                 anovas[[i]] <- anova(mod[[i]])
               }
               predicted <- as.data.frame(sapply(mod, predict))
               residuals <- as.data.frame(sapply(mod, residuals))
               colnames(predicted) <- colnames(df1)
               colnames(residuals) <- colnames(df1)
               
               
               
               names(mod) <- colnames(df1)
               names(results) <- colnames(df1)
               names(anovas) <- colnames(df1)
               
               mlm.mod@MLM <- mod
               mlm.mod@SUMMARY.MLM <- results
               mlm.mod@ANOVA.MLM <- anovas
               mlm.mod@PREDICTED <- predicted
               mlm.mod@RESIDUALS <- residuals
               mlm.mod@df1 <- df1
               mlm.mod@df2 <- df2
               
               mlm.mod
               return(mlm.mod)
             }
             
             else if(analysis == "mctree") {
               
               tre.new <- new("eco.mctree")
               
               tre <- list()
               prediction <- list()
               freq <- list()
               
               for(i in 1:ncol(df1))
               {
                 dep <- paste(colnames(df1)[i], "~", sep = "")
                 smod <- paste(dep, indep, sep = "")
                 smod <- as.formula(smod)
                 tre[[i]] <- party::ctree(smod, data = data, ...)
                 prediction[[i]] <- party::where(tre[[i]])
                 if(!is.null(fact)) {
                   freq[[i]] <- table(fact, prediction[[i]])
                 }
               }
               names(tre) <- colnames(df1)
               names(prediction) <- colnames(df1)
               if(!is.null(fact)) {
                 names(freq) <- colnames(df1)
               } 
               if(is.null(fact)){
                 freq[[1]] <- "factor not provided for computing frequencies"
               }
               
               predicted <- as.data.frame(sapply(tre, predict))
               residuals <- df1 - predicted
               
               tre.new@TREES <- tre
               tre.new@CLASSPREDICT <- prediction
               tre.new@PREDICTED <- predicted
               tre.new@RESIDUALS <- residuals
               tre.new@FREQUENCIES <- freq
               return(tre.new)
             }
           })
