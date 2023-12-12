#' dmrff.stats
#'
#' Calculate statistics for a set of genomic regions.
#'
#' Warning! Ensure that the order of the CpG sites corresponding to the the rows of `methylation`
#' match the order of the CpG sites corresponding to the other variables,
#' e.g. `estimate` and `chr`.
#' 
#' @param regions Data frame of genomic regions providing
#' chromosome (chr), start and end coordinates.
#' @param estimate Vector of EWAS effect estimates
#' (corresponds to rows of \code{methylation}).
#' @param se Vector of standard errors of the coefficients.
#' @param methylation Methylation matrix (rows=features, columns=samples).
#' @param chr Feature chromosome (corresponds to rows of \code{methylation}).
#' @param pos Feature chromosome position.
#' @param verbose If \code{TRUE} (default), then output status messages.
#' @return A data frame listing summary statistics for all genomic regions.
#'
#' @export
dmrff.stats <- function(regions, estimate, se, methylation, chr, pos, verbose=T) {
    stopifnot(is.data.frame(regions) && all(c("chr","start","end") %in% colnames(regions)))
    stopifnot(is.vector(estimate))
    stopifnot(is.vector(se))
    stopifnot(is.matrix(methylation))
    stopifnot(is.vector(chr))
    stopifnot(is.vector(pos))
    stopifnot(length(estimate) == length(se))
    stopifnot(length(estimate) == nrow(methylation))
    stopifnot(length(estimate) == length(chr))
    stopifnot(length(estimate) == length(pos))

    ## identify CpG sites in each region
    members <- region.members(regions, data.frame(chr=chr, pos=pos))

    ## scale summary stats as if methylation was standarized
    methylation.sd <- row.sds(methylation,na.rm=T)
    estimate <- estimate/methylation.sd
    se <- se/methylation.sd

    ## calculate region statistics
    stats <- do.call(rbind, parallel::mclapply(members, function(idx) {
        idx <- na.omit(idx)
        if (length(idx) == 0) return(c(B=NA,S=NA))
        ivwfe.stats(estimate[idx], se[idx], methylation[idx,,drop=F])
    }))
    
    regions$estimate <- stats[,"B"]
    regions$se <- stats[,"S"]
    regions$z <- stats[,"B"]/stats[,"S"]
    regions$p.value <- 2*pnorm(-abs(regions$z), lower.tail=T)
    regions$n <- sapply(members, length)
    
    regions
}
