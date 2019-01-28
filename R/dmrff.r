#' dmrff
#'
#' Identifying differentially methylated regions efficiently with power and control.
#'
#' @param estimate Vector of association effect estimates
#' (corresponds to rows of \code{methylation}).
#' @param se Vector of standard errors of the effect estimates.
#' @param p.value Vector of p-values.
#' @param methylation Methylation matrix (rows=features, columns=samples).
#' @param chr Feature chromosome (corresponds to rows of \code{methylation}).
#' @param pos Feature chromosome position.
#' @param p.cutoff Unadjusted p-value cutoff for membership in a candidate DMR
#' (Default: 0.05).
#' @param maxgap Maximum distance between consecutive features (Default: 500bp).
#' @param verbose If \code{TRUE} (default), then output status messages.
#' @return A data frame listing all candidate bumps and their summary statistics.
#' 
#' @examples
#'
#' dmrs <- dmrff(estimate, ## effect estimate for each CpG site
#'               se,       ## standard error of the estimate for each CpG site
#'               p.value,  ## p-value
#'               methylation, ## methylation matrix
#'               chr,      ## chromosome of each CpG site
#'               pos)      ## position of each CpG site
#'                      
#' dmrs[which(dmrs$p.adjust < 0.05),
#'      c("chr","start","end","n","B","S","z","p.value","p.adjust")]
#'
#' @export
dmrff <- function(estimate, se, p.value, methylation, chr, pos,
                  maxgap=500, p.cutoff=0.05, verbose=T) {
    stopifnot(is.vector(estimate))
    stopifnot(is.vector(se))
    stopifnot(is.vector(p.value))
    stopifnot(is.matrix(methylation))
    stopifnot(length(estimate) == length(se))
    stopifnot(length(estimate) == nrow(methylation))
    stopifnot(length(estimate) == length(p.value))
    stopifnot(length(estimate) == length(chr))
    stopifnot(length(estimate) == length(pos))

    
    candidates <- dmrff.candidates(estimate=estimate,
                                   p.value=p.value,
                                   chr=chr, 
                                   pos=pos,
                                   maxgap=maxgap,
                                   p.cutoff=p.cutoff,
                                   verbose=verbose)
 
    stats <- shrink.candidates(candidates$start.idx, candidates$end.idx,
                               function(start.idx,end.idx) {
                                   idx <- start.idx:end.idx
                                   ivwfe.getz(estimate[idx], se[idx], methylation[idx,,drop=F])
                               })

    full <- do.call(rbind, mclapply(1:nrow(stats), function(i) {
        idx <- stats$start.idx[i]:stats$end.idx[i]
        ivwfe.stats(estimate[idx], se[idx], methylation[idx,,drop=F])
    }))

    stats$B <- full[,"B"]
    stats$S <- full[,"S"]
    
    collate.stats(stats, chr, pos)
}


collate.stats <- function(stats, chr, pos) {   
    stats <- with(stats, data.frame(chr=chr[start.idx],
                                    start=pos[start.idx],
                                    end=pos[end.idx],
                                    n=end.idx-start.idx+1,
                                    start.idx=start.idx,
                                    end.idx=end.idx,
                                    start.orig=start.orig,
                                    end.orig=end.orig,
                                    z.orig=z.orig,
                                    p.orig=2*pnorm(-abs(z.orig), lower.tail=T),
                                    B=B,
                                    S=S,
                                    z=z,
                                    p.value=2*pnorm(-abs(z), lower.tail=T)))
    number.tests <- length(chr) + calculate.number.shrink.tests(stats)
    stats$p.adjust <- pmin(1, stats$p.value * number.tests)
    stats
}
