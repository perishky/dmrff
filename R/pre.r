#' dmrff.pre
#'
#' Construct an object for including this dataset in a DMR meta-analysis.
#'
#' @param estimate Vector of association effect estimates
#' (corresponds to rows of \code{methylation}).
#' @param se Vector of standard errors of the effect estimates.
#' @param methylation Methylation matrix (rows=features, columns=samples).
#' @param chr Feature chromosome (corresponds to rows of \code{methylation}).
#' @param pos Feature chromosome position.
#' @param maxsize Maximum number of CpG sites in a DMR.
#' @param verbose If \code{TRUE} (default), then output status messages.
#' @return A list object containing all information needed for inclusion in a meta-analysis.
#'
#' @export
dmrff.pre <- function(estimate, se, methylation, chr, pos, maxsize=20, verbose=T) {
    idx <- order(chr, pos)
    methylation <- impute.matrix(methylation)
    m <- methylation - rowMeans(methylation)
    ss <- sqrt(rowSums(m^2))
    rho <- do.call(cbind, mclapply(1:maxsize, function(size) {
        sapply(1:length(idx), function(i) {
            if (i + size <= length(idx)) {
                numer <- sum(m[idx[i],] * m[idx[i+size],])
                denom <- ss[idx[i]] * ss[idx[i+size]]
                numer/denom
            }
            else NA
        })
    })) ## 2 minutes, 3.5Mb x maxsize
    sites <- paste(chr, pos, sep=":")
    list(sites=sites,
         chr=chr,
         pos=pos,
         estimate=estimate,
         se=se,
         rho=rho)
}
   
