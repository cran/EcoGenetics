#' Combining ecogen objects by row
#' 
#' @param ... "ecogen" objects to combine. 
#' @param check_colnames Check for duplicated column names? Default TRUE.
#' @param check_rownames Check for duplicated row names? Default TRUE.
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' 
#' # split the object "eco" into a list of ecogen objects by population
#' x <- eco.split(eco, "pop", asList = TRUE)
#' 
#' # re-bind the objects
#' eco.r <- eco.rbind(eco)
#' 
#' # create a new objects with the first and second population
#' eco.r <- eco.rbind(x[[1]], x[[3]])
#' 
#' # duplicated row names are not allowed by eco.rbind with default options
#' eco2 <- eco
#' eco.rbind(eco, eco2)
#' 
#' eco.rbind(eco, eco2,check_rownames = FALSE)
#' }
#' 
#' @author Leandro Roser \email{learoser@@gmail.com}
#' 
#' @export

setGeneric("eco.rbind", 
           function(..., check_colnames = TRUE, check_rownames = TRUE)  {
             
             
             #-------ECOGEN OBJECTS OBTENTION---------------------------#
             # unlist dots
             u <- unlist(list(...))
             
             # ecogen objects
             u.ecogen <- u[sapply(u, is.ecogen)]
             # all ecogen objects
             u.ecogen <- c(u.ecogen)
             
             # checkpoint -> all objects passed are of class ecogen
             u.no_ecogen <- sapply(u.ecogen, function(x)!is.ecogen(x))
             if(any(u.no_ecogen)) {
               stop("non ecogen arguments passed to eco.rbind found")
             }
             
             #------PROCESSING THE DATA----------------------------------#
             
             # create a list of ecogen objects as lists, removing the slots A and OUT
             ecolist <- lapply(u.ecogen, function(x) as.int.list(x)[-c(4,8)])
             
             
             # checkpoint -> check column names
             if(check_colnames) {
               for(i in 1:6) {
                 #column names
                 cnames <- lapply(ecolist, function(x) colnames(x[[i]]))
                 cnames <- do.call(rbind, lapply(cnames, toupper))
                 # check null column names
                 cnames.null <- sapply(cnames, is.null)
                 if(any(cnames.null)) {
                   stop("null column names found")
                 }
                 # check different column names- no case sensitive
                 cnames <- apply(cnames, 2, unique)
                 # if non unique names, R generates a list
                 if(is.list(cnames)) {
                   stop("non unique column names found")
                 }
               }
             }
             # checkpoint -> check ploidy and ncod in the data
             areG <- lapply(ecolist, function(x) dim(x[[3]]))
             areG <- sapply(areG, function(x) x[[1]] * x[[2]])
             
             # at least two non empty data frames
             if(sum(areG != 0) > 1) {
               
               cuales <- which(areG != 0)
               
               #checkpoint -> ploidy
               checkG.ploidy <- lapply(ecolist[cuales], function(x) x[[7]]@ploidy)
               if(length(unique(unlist(checkG.ploidy))) != 1) {
                 stop("different ploidy levels found")
               }
               #checkpoint -> type
               checkG.type <- lapply(ecolist[cuales], function(x) x[[7]]@type)
               if(length(unique(unlist(checkG.type))) != 1) {
                 stop("different marker(s) type found")
               }
               #checkpoint -> ncod
               checkG.ncod <- lapply(ecolist[cuales], function(x) x[[7]]@ncod)
               if(length(unique(unlist(checkG.ploidy))) != 1) {
                 stop("different number of digits coding alleles found")
               }
             }
             
             #-------------------------------------------------------------------#
             # bind rows: function that bind rows of ecogen objects as list.
             # Duplicated row names present stops the fuction
             bind.rows <- function(ecolist, i) {

               # checkpoint -> duplicated row names not allowed
               if(check_rownames) {
                 names.x <- do.call(c, lapply(ecolist, function(y) rownames(y[[i]])))
                 if(any(duplicated(names.x))) {
                   stop(aue.formatLine("Duplicated row names found -
                        use check_rownames = FALSE to use data.frame methods for duplicated row names"))
                   # old method:
                   #return(data.frame()) 
                 }
               }
               
               do.call(rbind, lapply(ecolist, function(y) y[[i]]))
               
             }
             #---------------------END BIND ROWS---------------------------------#
             
             # generating the output as list by a call to <bind.rows>
             eco.out <- list()
             for(i in 1:6) {
               eco.out[[i]] <- bind.rows(ecolist, i)
             }
             
             # generating the output as ecogen. Maybe will be faster fill each
             # slot by hand. 
             z <- ecogen(XY = eco.out[[1]], P = eco.out[[2]],
                         G = eco.out[[3]], E = eco.out[[4]], 
                         S = eco.out[[5]], C = eco.out[[6]])
             
             # check validity
             validObject(z)
             
             z
           })
