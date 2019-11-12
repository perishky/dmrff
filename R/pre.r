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
    stopifnot(is.vector(estimate))
    stopifnot(is.vector(se))
    stopifnot(is.matrix(methylation))
    stopifnot(is.vector(chr))
    stopifnot(is.vector(pos))
    stopifnot(length(estimate) == length(se))
    stopifnot(length(estimate) == nrow(methylation))
    stopifnot(length(estimate) == length(chr))
    stopifnot(length(estimate) == length(pos))
    
    # sort input by chromosomal position
    idx <- order(chr,pos)
    sorted <- identical(idx, 1:length(idx))
    if (!sorted) {
        estimate <- estimate[idx]
        se <- se[idx]
        chr <- chr[idx]
        pos <- pos[idx]
        methylation <- methylation[idx,,drop=F]
    }

    ## calculate rho (CpG correlation matrix)
    m <- methylation - rowMeans(methylation, na.rm=T)
    mm <- rowSums(m^2, na.rm=T)
    ss <- sqrt(mm)
    rho <- do.call(cbind, parallel::mclapply(1:maxsize, function(size) {
        sapply(1:length(chr), function(i) {
            if (i + size <= length(idx)) {
                numer <- sum(m[i,] * m[i+size,])
                denom <- ss[i] * ss[i+size]
                numer/denom
            }
            else NA
        })
    }))
    n <- rowSums(!is.na(m))
    sd <- sqrt(mm/(n-1))
    
    sites <- paste(chr, pos, sep=":")
    list(sites=sites,
         chr=chr,
         pos=pos,
         estimate=estimate,
         se=se,
         rho=rho,
         sd=sd)
}
   
