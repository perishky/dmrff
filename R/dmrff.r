#' dmrff
#'
#' Identifying differentially methylated regions efficiently with power and control.
#'
#' Warning! Ensure that the order of the CpG sites corresponding to the the rows of `methylation`
#' match the order of the CpG sites corresponding to the other variables,
#' e.g. `estimate` and `chr`.
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
#' @return A data frame listing all candidate regions and their summary statistics.
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
#' dmrs[which(dmrs$p.adjust < 0.05 & dmrs$n > 1),]
#'
#' @export
dmrff <- function(estimate, se, p.value, methylation, chr, pos,
                  maxgap=500, p.cutoff=0.05, minmem=T, verbose=T) {
    stopifnot(is.vector(estimate))
    stopifnot(is.vector(se))
    stopifnot(is.vector(p.value))
    stopifnot(is.matrix(methylation))
    stopifnot(is.vector(chr))
    stopifnot(is.vector(pos))
    stopifnot(length(estimate) == length(se))
    stopifnot(length(estimate) == nrow(methylation))
    stopifnot(length(estimate) == length(p.value))
    stopifnot(length(estimate) == length(chr))
    stopifnot(length(estimate) == length(pos))

    # order input by chromosomal position
    ord.idx <- order(chr,pos)
    estimate <- estimate[ord.idx]
    se <- se[ord.idx]
    p.value <- p.value[ord.idx]
    chr <- chr[ord.idx]
    pos <- pos[ord.idx]
    sorted <- identical(ord.idx,1:length(ord.idx))
    
    # identify candidate regions
    candidates <- dmrff.candidates(
        estimate=estimate,
        p.value=p.value,
        chr=chr, 
        pos=pos,
        maxgap=maxgap,
        p.cutoff=p.cutoff,
        verbose=verbose)

    if (is.null(candidates)) {
        return (NULL)
    }

    ## reduce dataset to just cover the candidates
    if (minmem) {
        red.end.idx <- cumsum(candidates$end.idx-candidates$start.idx+1)
        red.start.idx <- c(1,head(red.end.idx,-1)+1)
        red.idx <- unlist(lapply(1:nrow(candidates), function(i) {
            with(candidates, start.idx[i]:end.idx[i])
        }))
        estimate <- estimate[red.idx]
        se <- se[red.idx]
        p.value <- p.value[red.idx]
        methylation <- methylation[ord.idx[red.idx],,drop=F]
        gc()
    } else {
        red.idx <- 1:nrow(methylation)
        red.start.idx <- candidates$start.idx
        red.end.idx <- candidates$end.idx
        methylation <- methylation[ord.idx,,drop=F]
    }

    ## scale summary stats as if methylation was standarized
    methylation.sd <- row.sds(methylation,na.rm=T)
    estimate <- estimate/methylation.sd
    se <- se/methylation.sd
    
    # identify sub-regions that maximize statistical significance
    stats <- shrink.candidates(
        red.start.idx, red.end.idx,
        function(start.idx,end.idx) {
            idx <- start.idx:end.idx
            ivwfe.getz(estimate[idx], se[idx], methylation[idx,,drop=F])
        })
    
    # calculate B and S statistics for each region (recall z=B/S)
    full <- do.call(rbind, parallel::mclapply(1:nrow(stats), function(i) {
        idx <- stats$start.idx[i]:stats$end.idx[i]
        ivwfe.stats(estimate[idx], se[idx], methylation[idx,,drop=F])
    }))

    stats$estimate <- stats$B <- full[,"B"]
    stats$se <- stats$S <- full[,"S"]
    stats$start.idx <- red.idx[stats$start.idx]
    stats$end.idx <- red.idx[stats$end.idx]
    stats$start.orig <- red.idx[stats$start.orig]
    stats$end.orig <- red.idx[stats$end.orig]
    
    collate.stats(stats, chr, pos, simple=!sorted)
}

collate.stats <- function(stats, chr, pos, simple=F) {   
    stats <- with(stats, {
        data.frame(
            chr=chr[start.idx],
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
            estimate=estimate,
            se=se,
            z=z,
            p.value=2*pnorm(-abs(z), lower.tail=T))
    })
    number.tests <- length(chr) + calculate.number.shrink.tests(stats)
    stats$number.tests <- number.tests
    stats$p.adjust <- pmin(1, stats$p.value * number.tests)
    if (simple)
        stats$start.idx <- stats$end.idx <- stats$start.orig <- stats$end.orig <- stats$z.orig <- stats$p.orig <- NULL
    stats
}

