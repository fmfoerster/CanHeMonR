% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PWN_functions.R
\name{samp_df}
\alias{samp_df}
\title{Avoid Data Frames With Too Many Rows}
\usage{
samp_df(df, maxsamp)
}
\arguments{
\item{df}{a data frame}

\item{maxsamp}{the maximum}
}
\value{
a dataframe with no more than maxsamp rows
}
\description{
Ensure that a dataframe has <= maxsamp rows.
If it doesn't, sample a subset (n=maxsamp) rows from a data frame
}

