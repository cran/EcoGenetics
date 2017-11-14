
#' Sliding window for matrix data
#' 
#' @description This program applies a function defined by the user, 
#' using a moving window (circle area or square) and assigning
#' the value to the focal pixel.
#' @param mat Input raster matrix.
#' @param r Half a side for square, radius for circle, diagonal length for rhombus.
#' @param slide Number of elements between two focal pixels, for column 
#' and row dimension
#' @param fun Function to apply in each focal pixel.
#' @param window Window type. Default "square".
#' @param within Should the function be computed in focal pixels of the borders, only if 
#' the area is within the matrix? Default TRUE.
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' ras <- matrix(eco[["P"]][, 1], 15, 15)
#' image(ras)
#' ras.square <- eco.slide.matrix(ras, 1, 1, mean, "square")
#' image(ras.square)
#' 
#' # or allowing more control over the function:
#' ras.square <- eco.slide.matrix(ras, r = 3, slide = 1, function(x) mean(x, na.rm = TRUE), "square")
#' image(ras.square)
#'
#'# sliding a circle:
#' ras.circle <- eco.slide.matrix(ras, r = 3, slide = 1, mean, "circle", within = FALSE)
#' image(ras.circle) 
#'
#' }
#' 
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @export

eco.slide.matrix <- function(mat, r, slide, fun, 
                            window = c("square", "circle"),
                            within = TRUE) {
  x <- match.arg(window)
  
  # function selection
  if(x == "square") {
    fun.local <- function(obj, rad, x.ptr, y.ptr) aue.square(obj, rad, x.ptr, y.ptr)
  } else if(x == "circle") {
    fun.local <- function(obj, rad, x.ptr, y.ptr) aue.point(obj, rad, x.ptr, y.ptr) 
  } 
  #------------------------------------------------#
  
  # create a sequence of row / column indices
  iseq <- seq(from = slide, to = ncol(mat), by = slide)
  # create a sequence of column indices
  jseq <- seq(from = slide, to = nrow(mat), by = slide)
  
  # remove min or max value of jseq if are lower/higher than matrix dimension
  if(within) {
    if(min(iseq) - r < 1) {
      iseq <- iseq[-1]
    }
    if(max(iseq)+ r > ncol(mat)) {
      iseq <- iseq[-length(iseq)]
    }
    if(min(jseq) - r < 1) {
      jseq <- jseq[-1]
    }
    if(max(jseq) + r > ncol(mat)) {
      jseq <- jseq[-length(jseq)]
    }
  }
  
  # pre-allocate memory
  out <- matrix(0, length(iseq), length(jseq))
  colnames(out) <- iseq
  rownames(out) <- jseq
  # run loop
  i.temp <- 1
  for(i in  iseq) {
    j.temp <- 1
    for(j in jseq) {
      area <- which(fun.local(mat, r, i, j) != 0)
      out[i.temp, j.temp] <- fun(mat[area])
      j.temp <- j.temp + 1
    }
    i.temp <- i.temp + 1
  }
  t(out)
}
