% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eco.slide.matrix.R
\name{eco.slide.matrix}
\alias{eco.slide.matrix}
\title{Sliding window for matrix data}
\usage{
eco.slide.matrix(
  mat,
  r,
  slide,
  fun,
  window = c("square", "circle"),
  within = TRUE
)
}
\arguments{
\item{mat}{Input raster matrix.}

\item{r}{Half a side for square, radius for circle, diagonal length for rhombus.}

\item{slide}{Number of elements between two focal pixels, for column 
and row dimension}

\item{fun}{Function to apply in each focal pixel.}

\item{window}{Window type. Default "square".}

\item{within}{Should the function be computed in focal pixels of the borders, only if 
the area is within the matrix? Default TRUE.}
}
\description{
This program applies a function defined by the user, 
using a moving window (circle area or square) and assigning
the value to the focal pixel.
}
\examples{

\dontrun{

data(eco.test)
ras <- matrix(eco[["P"]][, 1], 15, 15)
image(ras)
ras.square <- eco.slide.matrix(ras, 1, 1, mean, "square")
image(ras.square)

# or allowing more control over the function:
ras.square <- eco.slide.matrix(ras, r = 3, slide = 1, function(x) mean(x, na.rm = TRUE), "square")
image(ras.square)

# sliding a circle:
ras.circle <- eco.slide.matrix(ras, r = 3, slide = 1, mean, "circle", within = FALSE)
image(ras.circle) 

}

}
\author{
Leandro Roser \email{learoser@gmail.com}
}
