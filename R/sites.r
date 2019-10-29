#' dmrff.sites
#'
#' Link genomic locations to a set of genomic regions.
#'
#' @param regions Data frame of genomic regions providing
#' chromosome (chr), start and end coordinates.
#' @param chr A vector providing the chromosome of each location.
#' @param pos A vector providing the chromosomal position of each location.
#' @return A data frame identifying the region containing each genomic location.
#'
#' @examples
#'
#' dmrs <- dmrff(estimate, se, p.value, methylation, chr, pos)
#' dmrs <- dmrs[which(dmrs$p.adjust < 0.05 & dmrs$n > 1),]
#' dmr.sites <- dmrff.sites(dmrs, chr, pos, estimate, se, p.value)
#' dmr.sites$dmr.z <- dmrs$z[dmr.sites$region]
#' 
#' @export
dmrff.sites <- function(regions, chr, pos) {
    stopifnot(is.data.frame(regions) && all(c("chr","start","end") %in% colnames(regions)))
    stopifnot(is.vector(chr))
    stopifnot(is.vector(pos))
    stopifnot(length(chr) == length(pos))

    sites <- data.frame(site=1:length(chr),
                        chr=chr,
                        pos=pos)
                       
    members <- region.members(regions, sites)

    idx <- unlist(members)
    if (length(idx) > 0)
        data.frame(region=rep(1:nrow(regions), sapply(members, length)),
                   sites[idx,],
                   stringsAsFactors=F)
    else
        NULL
}
