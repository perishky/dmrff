% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cohort.r
\name{dmrff.cohort}
\alias{dmrff.cohort}
\title{dmrff.cohort}
\usage{
dmrff.cohort(object, maxgap = 500, p.cutoff = 0.05, verbose = T)
}
\arguments{
\item{object}{Object generated by \code{\link{dmrff.pre}} for the dataset.}

\item{maxgap}{Maximum distance between consecutive features (Default: 500bp).}

\item{p.cutoff}{Unadjusted p-value cutoff for membership in a candidate DMR
(Default: 0.05).}

\item{verbose}{If \code{TRUE} (default), then output status messages.}
}
\value{
A data frame listing all candidate regions and their summary statistics.
}
\description{
Identify differentially methylated regions 
within an individual dataset with a `pre` object.
}
\examples{
pre <- dmrff.pre(est, se, p, meth, ...)
dmrs <- dmrff.cohort(pre)
dmrs[which(dmrs$p.adjust < 0.05 & dmrs$n > 1), ]

}
