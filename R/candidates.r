#' dmrff.candidates
#'
#' Identify candidate differentially methylated regions from EWAS summary statistics.
#'
#' @param estimate Vector of EWAS effect estimates
#' (corresponds to rows of \code{methylation}).
#' @param p.value Vector of p-values.
#' @param chr Feature chromosome (corresponds to rows of \code{methylation}).
#' @param pos Feature chromosome position.
#' @param p.cutoff Unadjusted p-value cutoff for membership in a candidate DMR
#' (Default: 0.05).
#' @param maxgap Maximum distance between consecutive features (Default: 500bp).
#' @param verbose If \code{TRUE} (default), then output status messages.
#' @return A data frame listing all candidate regions.
#'
#' @export
dmrff.candidates <- function(estimate, p.value, chr, pos,
                             maxgap=500, p.cutoff=0.05, verbose=T) {
    stopifnot(is.vector(estimate))
    stopifnot(is.vector(p.value))
    stopifnot(is.vector(chr))
    stopifnot(is.vector(pos))
    stopifnot(length(estimate) == length(p.value))
    stopifnot(length(estimate) == length(chr))
    stopifnot(length(estimate) == length(pos))
    
    if (is.na(p.cutoff) || !is.numeric(p.cutoff) || p.cutoff > 1 || p.cutoff < 0)
        stop("'p.cutoff' is a p-value cutoff so must be between 0 and 1")

    # sort by chromosomal position
    idx <- order(chr,pos)
    sorted <- identical(idx, 1:length(idx))
    if (!sorted) {
        estimate <- estimate[idx]
        p.value <- p.value[idx]
        chr <- chr[idx]
        pos <- pos[idx]
    }
    
    sig.idx <- which(p.value < p.cutoff)
    if (length(sig.idx) == 0) {
        if (verbose)
            msg("No candidate regions found")
        return(NULL)
    }
    
    candidates <- bh.regionFinder(x = sign(estimate[sig.idx]),
                                  chr = chr[sig.idx],
                                  pos = pos[sig.idx],
                                  maxGap=maxgap,
                                  cutoff = c(-1,1),
                                  verbose = FALSE)
    
    if (verbose)
        msg("Found ", nrow(candidates), " candidate regions.")
    
    candidates$start.idx <- sig.idx[candidates$indexStart]
    candidates$end.idx <- sig.idx[candidates$indexEnd]
    candidates$candidate <- 1:nrow(candidates)

    if (sorted) 
        candidates[,c("chr","start","end","candidate","start.idx","end.idx")]
    else
        candidates[,c("chr","start","end","candidate")]
}
