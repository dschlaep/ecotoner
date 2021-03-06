% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilitities_grids.R
\name{rasterToPoints}
\alias{rasterToPoints}
\title{rasterToPoints}
\usage{
rasterToPoints(x, fun = NULL, spatial = FALSE, not_convertible = list(NA),
  ...)
}
\arguments{
\item{x}{A Raster* object }

\item{fun}{Function to select a subset of raster values}

\item{spatial}{Logical. If \code{TRUE}, the function returns a SpatialPointsDataFrame object }

\item{not_convertible}{A list of single elements, e.g., NA, 10, and of vectors of length two, e.g., c(0, 10).}

\item{...}{Additional arguments. Currently only \code{progress} to specify a progress bar. "text", "window", or "" (the default, no progress bar)}
}
\description{
Identical to \code{\link{raster::rasterToPoints}} if \code{not_convertible} is \code{NA}.
}
\details{
rasterToPoints does not convert cells with values that are contained as single values in not_convertible or that are within and including the range of min-max of the length two vectors of not_convertible.
}
\examples{
if (requireNamespace("raster")) {
  r <- raster::raster(nrow = 5, ncol = 5)
  r[] <- 1:(raster::ncell(r))
  r[1:5] <- NA
  
  identical(rasterToPoints(r, not_convertible = NA), raster::rasterToPoints(r))
  rasterToPoints(r, not_convertible = NULL)
  rasterToPoints(r, not_convertible = 25)
  rasterToPoints(r, not_convertible = list(NA, 25, c(7.1, 20.5)))
}
}
\seealso{
\code{\link{raster::rasterToPoints}} which this function simplifies
}

