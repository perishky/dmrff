% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.r
\name{dmrff.stats}
\alias{dmrff.stats}
\title{dmrff.stats}
\usage{
dmrff.stats(regions, estimate, se, methylation, chr, pos, verbose = T)
}
\arguments{
\item{regions}{Data frame of genomic regions providing
chromosome (chr), start and end coordinates.}

\item{estimate}{Vector of EWAS effect estimates
(corresponds to rows of \code{methylation}).}

\item{se}{Vector of standard errors of the coefficients.}

\item{methylation}{Methylation matrix (rows=features, columns=samples).}

\item{chr}{Feature chromosome (corresponds to rows of \code{methylation}).}

\item{pos}{Feature chromosome position.}

\item{verbose}{If \code{TRUE} (default), then output status messages.}
}
\value{
A data frame listing summary statistics for all genomic regions.
}
\description{
Calculate statistics for a set of genomic regions.
}
\details{
Warning! Ensure that the order of the CpG sites corresponding to the the rows of `methylation`
match the order of the CpG sites corresponding to the other variables,
e.g. `estimate` and `chr`.
}
