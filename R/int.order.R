#' Ordering the rows of the data frames contained in an ecogen object
#' 
#' @details This program generates an ecogen object with the rows of all
#' the data frames ordered in reference to the row names of the XY data frame.
#' This is useful when the data frames are loaded into the ecogen object,
#' but were not ordered previously. Also, this tool can be useful 
#' for reorder rows when is needed. 
#' First, the reference data frame in the slot XY will be in the desired row order. 
#' This program then aligns all the data frames by coincidence of row names
#' with those in the slot XY.  
#' 
#' @param eco Object of class "ecogen".
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' eco1 <- eco
#' eco1[["P"]] <- eco[["P"]][sample(1:225), ]  #object with shuffled rows
#' eco1[["E"]] <- eco[["E"]][sample(1:225), ]
#' ordered <- int.order(eco1)
#' head(ordered[["P"]]); head(eco[["P"]])
#' 
#' }
#' 
#' @author Leandro Roser \email{learoser@@gmail.com}
#' 
#' @keywords internal

setGeneric("int.order", 
           function(eco) {
             
             refnames <- eco@ATTR$names
             ord <- 0
             #function to order rows
             .order <- function(.refnames, df_to_order) {
               
               this_order <- match(.refnames, rownames(df_to_order), nomatch = 0)
               if(length(this_order) != length(.refnames)) {
                 stop("invalid object: some row names not matching")
               }
               
               if(!all(this_order == seq_len(length(.refnames))) && !all(this_order == 0)) {
                 df_to_order <- df_to_order[this_order, , drop = FALSE]
                 ord <<- ord + 1
               }
               df_to_order
             }
             
             if(length(refnames) != 0) {
               
               eco@XY <- .order(refnames, eco@XY)
               eco@P <- .order(refnames, eco@P)
               eco@G <- .order(refnames, eco@G)
               eco@A <- .order(refnames, eco@A)
               eco@E <- .order(refnames, eco@E)
               eco@S <- .order(refnames, eco@S)
               eco@C <- .order(refnames, eco@C)
               
               #ordering missing cells
               dum <- as.matrix(eco@A - eco@A)
               dum[eco@INT@missing.cells] <- NA
               dum <- .order(refnames, dum)
               eco@INT@missing.cells <- which(is.na(dum))
               
               if(ord != 0) {
                 message("Note: ordered rows")
               }
               
             }
             
             eco
             
           })