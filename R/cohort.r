#' dmrff.cohort
#'
#' Identify differentially methylated regions 
#' within an individual dataset with a `pre` object.
#'
#' @param object Object generated by \code{\link{dmrff.pre}} for the dataset.
#' @param p.cutoff Unadjusted p-value cutoff for membership in a candidate DMR
#' (Default: 0.05).
#' @param maxgap Maximum distance between consecutive features (Default: 500bp).
#' @param verbose If \code{TRUE} (default), then output status messages.
#' @return A data frame listing all candidate regions and their summary statistics.
#' 
#' @examples
#' pre <- dmrff.pre(est, se, p, meth, ...)
#' dmrs <- dmrff.cohort(pre)
#' dmrs[which(dmrs$p.adjust < 0.05 & dmrs$n > 1), ]
#' 
#' @export
dmrff.cohort <- function(object, maxgap=500, p.cutoff=0.05, verbose=T) {
    if (!is.list(object)
        && all(c("estimate","se","chr","pos","sites","rho") %in% names(object)))
        stop("'object' was not created by dmrff.pre()")
    
    idx <- order(object$chr, object$pos)
    sorted <- identical(idx, 1:length(idx))
    if (!sorted) { ## support for a previous version of dmrff that did not sort
        object$sites <- object$sites[idx]
        object$chr <- object$chr[idx]
        object$pos <- object$pos[idx]
        object$estimate <- as.numeric(object$estimate[idx])
        object$se <- as.numeric(object$se[idx])
    }
    object$p.value <- 2*pnorm(abs(object$estimate/object$se), lower.tail=F)
        
    candidates <- dmrff.candidates(estimate=object$estimate,
                                   p.value=object$p.value,
                                   chr=object$chr, 
                                   pos=object$pos,
                                   maxgap=maxgap,
                                   p.cutoff=p.cutoff,
                                   verbose=verbose)
 
    stats <- shrink.candidates(candidates$start.idx, candidates$end.idx,
                               function(start.idx,end.idx) {
                                   idx <- start.idx:end.idx
                                   if (length(idx) > ncol(object$rho)) return(0) 
                                   ivwfe.getz(object$estimate[idx], object$se[idx],
                                              rho=extract.rho(object$rho[idx,,drop=F]))
                               })

    full <- do.call(rbind, mclapply(1:nrow(stats), function(i) {
        idx <- stats$start.idx[i]:stats$end.idx[i]
        if (length(idx) > ncol(object$rho)) return(c(B=0,S=1)) 
        ivwfe.stats(object$estimate[idx], object$se[idx],
                    rho=extract.rho(object$rho[idx,,drop=F]))
    }))

    stats$estimate <- stats$B <- full[,"B"]
    stats$se <- stats$S <- full[,"S"]
    
    collate.stats(stats, object$chr, object$pos, simple=!sorted)
}

extract.rho <- function(pre) {
    n <- nrow(pre)
    rho <- diag(x=1.05, n)
    if (n > 1) 
        for (i in 1:(n-1)) {
            vals <- pre[i,1:(n-i)]
            rho[i,(i+1):n] <- vals
            rho[(i+1):n,i] <- vals
        }
    rho
}

