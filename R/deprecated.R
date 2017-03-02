#' Functions deprecated in EcoGenetics version 1.2.0-2
#'@param ... parameters
#' @export

eco.order <- function(...) {
  stop("eco.order is deprecated for Ecogenetics >= 1.2.1. 
        Ordering is now an option in the 'ecogen' function")
}

aue.filter <- function(...) {
  stop("aue.filter is deprecated for Ecogenetics >= 1.2.0-2. Use eco.slide")
}

aue.idig <- function(...) {
  stop("aue.idig is deprecated for Ecogenetics >= 1.2.0-2. Use eco.format")
}

aue.char2num <- function(...) {
  stop("aue.char2num is deprecated for Ecogenetics >= 1.2.0-2. Use eco.format")
}

eco.2columns <- function(...) {
  stop("eco.2columns is deprecated for Ecogenetics >= 1.2.0-2. Use eco.convert")
}

eco.append <- function(...) {
  stop("eco.append is deprecated for Ecogenetics >= 1.2.0-2. 
        The function has been improved with accessors. Use the 
        accessor <ecoslot.OUT>.
       See help(\"EcoGenetics accessors\"), Details and Examples")
}

eco.2geneland  <- function(...){
  stop("eco.2geneland is deprecated for Ecogenetics >= 1.2.0-4. Use ecogen2geneland")
}

eco.2genepop  <- function(...){
  stop("eco.2geneland is deprecated for Ecogenetics >= 1.2.0-4. Use ecogen2genepop")
}

eco.2genind <- function(...){
  stop("eco.2geneland is deprecated for Ecogenetics >= 1.2.0-4. Use ecogen2genind")
}

eco.2gstudio <- function(...){
  stop("eco.2geneland is deprecated for Ecogenetics >= 1.2.0-4. Use ecogen2gstudio")
}

eco.2hierfstat <- function(...){
  stop("eco.2geneland is deprecated for Ecogenetics >= 1.2.0-4. Use ecogen2hierfstat")
}

eco.2spagedi <- function(...){
  stop("eco.2geneland is deprecated for Ecogenetics >= 1.2.0-4. Use ecogen2spagedi")
}

eco.genepop2df <- function(...){
  stop("eco.2geneland is deprecated for Ecogenetics >= 1.2.0-4. Use genepop2ecogen")
}






