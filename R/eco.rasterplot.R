
#' rasterplot graphs
#' @param x data matrix (raster)
#' @param filter optional data matrix used as filter
#' @param condition condition used to filter data
#' @param limits values limits used for computing the data gradient for the plot
#' @param title plot title
#' @param z.name name for the legend
#' @param vertical should be partitioned the populations on the x axis? Default
#' TRUE. 
#' @description This function generates a multivariate plot for 
#' a data matrix (raster), with an option for filtering the data
#' and to graph using groups. The resterplot graph is a flexible tools
#' for multiple data sources (environmental, genetic, phenotypic, etc.).
#' 
#' 
#' @examples
#' data(eco.test)
#' 
#' # using the ecogen object "eco" to perform a multiple-lsa
#' con <- eco.weight(eco[["XY"]], method = "knearest", k = 4, row.sd = TRUE)
#' test.lsa <- eco.lsa(eco[["P"]], con = con, method = "I", nsim = 99, multi = "matrix")
#' 
#' # the default method for this object, is a resterplot
#' plot(test.lsa)
#' 
#' # adding a factor
#' test.lsa <- eco.lsa(eco[["P"]], con = con, method = "I",
#' nsim = 99, multi = "matrix", pop = eco[["S"]][,1])
#' plot(test.lsa)
#' 
#' # The generic rasterplot method requires a data matrix, and, as option, a condition 
#' # and a filter matrix. The condition is an expression, containing the word "filter" and 
#' # logical elements, e.g., "filter < 50", "filter <50 || filter > 2", etc. ). 
#' # Filter is used as a logical matrix (TRUE-FALSE, in relation to the passed condition),
#' # for filtering the data. If a condition is passed but not a filter matrix, the condition
#' # is applied over the data matrix, also using the word "filter". 
#' # Internally, the multi.lsa plot uses three fundamental elements. 
#' # - a data matrix: in the example, ecoslot.OBS(test.lsa)
#' #  a filter matrix: in the example, ecoslt.PVAL(test.lsa); i.e., the data matrix will be filtered
#' # by P-value using the third element, an expresion.
#' # - an expression: in the example: "filter < 0.05"
#'  
#'  # Combining the three elements, the multivariate plot can be manually constructed:
#'  my.plot <- eco.rasterplot(x= ecoslot.OBS(test.lsa), filter = ecoslot.PVAL(test.lsa), condition = "filter < 0.05")
#'  
#'  
#'  # add population
#'  my.plot <- eco.rasterplot(x= ecoslot.OBS(test.lsa), filter = ecoslot.PVAL(test.lsa), 
#'  condition = "filter < 0.05", grp = ecoslot.POP(test.lsa))
#'  
#'  
#'  # extra manipulation with ggplot2 syntax (ggplot2 commands allowed by rasterplot)
#'  ## rotate plot
#'  my.plot + coord_flip()
#'  
#'  ## change design
#'  my.plot + theme_grey()
#'  
#'  
#'  # using the data as filter
#'  eco.rasterplot(x= ecoslot.OBS(test.lsa), filter = ecoslot.OBS(test.lsa), 
#'  condition = "filter > 0 & filter < 3")
#'  
#'  
#'  # example of bad syntax (incorrect use of && over matrices)
#'  eco.rasterplot(x= ecoslot.OBS(test.lsa), filter = ecoslot.OBS(test.lsa), 
#'  condition = "filter > 0 && filter < 3")
#'  
#' @export



setGeneric("eco.rasterplot",  
           
           function(x, 
                    filter = NULL,
                    condition = NULL,
                    grp = NULL,
                    limits = NULL,
                    title = NULL,
                    z.name = NULL,
                    vertical  = TRUE,
                    ...) {
             
             x <- aue.image2df(x)
             
             # group configuration
             if(!is.null(grp)) {
               grp <- aue.image2df(grp)
               x <- cbind(x, grp[, 3])
               colnames(x)[4] <- "Group"
               if(vertical) {
                 x <- x[order(x[, 4]), ]
               } else {
                 x <- x[order(x[, 4], decreasing = TRUE), ]
               }
             }
    
             colnames(x)[1:2] <- c("Sample", "Variable")
             
             # y limits
             minplot <- min(x[, 2])
             maxplot <- max(x[, 2])
             
             if(!is.null(condition)) {
               
               # filter present /absent. If absent, the filter is the data matrix
               if(!is.null(filter)){
                 filter <- as.vector(filter)
               } else {
                 filter <- x
               }
               
               # conditional expression checkpoint
               ## obtain expression-matrix or character.
               cond <- deparse(substitute(condition))
               cond <- gsub("\\\"", "", cond)
               cond <- gsub(" ", "", cond)
               
               # evaluate expression
               cond <- eval(parse(text = cond))
               
               # check consistency
               if(length(cond) != length(filter)) {
                 stop("bad condition syntax")
               }
               
               # all FALSE -> stop
               if(all(cond == FALSE)) {
                 stop("the condition is not satisfied by any filter value")
               }
               
               x <- x[cond, ]
             }
             
             # na remotion
             na.data <- which(is.na(x[, 3]))
             if(any(na.data)) {
               x <- x[-na.data, , drop = FALSE]
             }
               
             # limits configuration
             if(is.null(limits)) {
               limits = c(min(x[, 3]), max(x[, 3]))
             }
             
             # OUTPUT CREATION
             
             if(is.null(title)) {
               title <- ""
             } 
             
             if(is.null(z.name)) {
               z.name <- "   z"
             } 
             
             out <- ggplot2::ggplot(x, ggplot2::aes(Sample, Variable, fill =  z)) + 
               ggplot2::geom_raster() +
               ggplot2::labs(title = title)+
               scale_fill_gradient2(name= z.name,space = "Lab",na.value = "white",
                                    high= scales::muted("red"),
                                    low = scales::muted("blue"), limits = limits) +
               ggplot2::theme_bw()+
               ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
                              axis.title = ggplot2::element_text(size = 14, face = "bold"), 
                              legend.position = "right") + 
               ggplot2::scale_y_discrete(expand = c(0.1, 0), limits = c(minplot:maxplot), 
                                         breaks = scales::pretty_breaks())
             
             
             if(!is.null(grp)) {
               if(vertical) {
               out <- out + ggplot2::facet_grid(.~ Group , scales = "free") 
               } else {
               out <- out + ggplot2::facet_grid(Group ~., ,scales = "free") 
               }
             }
        
             out

             
           })

#-------------------------------------------------------------------#
#' rasterplot graph for eco.lsa results
#' 
#' @param x eco.multilsa object returned by \code{\link{eco.lsa}} or 
#' @param significant plot only significant results?  Default TRUE
#' @param rescaled plot the rescaled observed values ([-1,1] range)?
#' @param alpha threshold P value for results with permutation tests. default = 0.05.
#' @param limits values limits used for computing the data gradient for the plot
#' @param title plot title
#' @param z.name name for the legend
#' @description Plot method for local spatial analysis
#' @rdname eco.multilsa-method
#' @aliases plot,eco.multilsa-method
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @seealso  \code{\link{eco.lsa}}
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' 
#' @exportMethod plot


setMethod("eco.rasterplot", 
          
          c("eco.multilsa", 
            "missing",
            "missing"),
          
          function(x,
                   grp =  NULL,
                   limits = NULL,
                   title = NULL,
                   z.name = NULL,
                   vertical = TRUE,
                   significant = TRUE,
                   rescaled = FALSE,
                   alpha = 0.05) {
            
            
            if(rescaled) {
              values <- x@OBS.RES
            } else {
            values <- x@OBS
            }
            
            filter <- x@POP
            
            # significant configuration
            if(significant) {
              if(x@TEST == "permutation") {
                values[x@PVAL > alpha] <- NA
              } else {
                values[values > x@LWR & values < x@UPPR] <- NA
              }
              
              if(all(is.na(values))) {
                msg <- paste("No significant results. Use significant = FALSE to view all the results")
                return(message(msg))
              } 
            } 
            
            # legend configuration
            sel <- match(x@METHOD,  c("G", "G*", "I", "C"))
            title <- c("Getis Ord's G", "Getis Ord's G*", 
                        "local Moran's I", "local Geary's C")
            title <- title[sel]
            
            message(paste("plot options: significant =", significant))
            message(paste("plot options: rescaled =", rescaled))
            
            callGeneric(x = values, filter = NULL, condition =  NULL,
                        grp = x@POP, limits = limits,
                        title = title,  z.name = paste("  ", x@METHOD),
                        vertical = vertical)
            
          })
