% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spectral_indices.r
\name{mod_SR}
\alias{mod_SR}
\title{Modified Simple Ratio}
\usage{
mod_SR(df, outp_fname = NULL)
}
\arguments{
\item{df}{A data frame where columns represent measurements in a single wavelength, and columns are named following Quantalab's conventions}

\item{outp_fname}{In case the input is raster data, this is the optional output filename to write the result to}
}
\value{
A vector with the value of the index
}
\description{
Calculate Modified Simple Ratio vegetation index
}
\references{
Chen 1996
}
