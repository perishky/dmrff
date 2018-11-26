#' dmrff.annotate
#'
#' Annotate a set of regions with feature annotations.
#'
#' For example, the regions could be differentially methylated regions
#' and the features could be CpG sites.
#'
#' @param regions Data frame listing the regions to annotate.
#' Must have columns "chr", "start" and "end" to specify genomic regions.
#' @param annotation Data frame listing the features that will be used to
#' annotate the genomic regions.
#' Must have columns "chr" and "pos" to identify genomic locations.
#' @return A data frame list all features contained in the regions along
#' with the regions they belong to.
#'
#' @examples
#'
#' methylation <- ... ## methylation matrix
#' ewas <- data.frame(chr=..., ## chromosome of each CpG site
#'                    pos=..., ## chromosomal position of each CpG site
#'                    symbol=..., ## symbol of gene linked to each CpG site
#'                    ..., ## other CpG site annotations
#'                    estimate=..., ## EWAS effect estimate
#'                    se=...,       ## EWAS standard error of the estimate 
#'                    p.value=...)  ## EWAS p-value
#'                 
#' dmrs <- with(ewas, dmrff(estimate, se, p.value, methylation, chr, pos))
#' dmrs <- dmrs[which(dmrs$p.adjust < 0.05),]
#' dmrs <- dmrff.annotate(dmrs, annotation)
#' 
#' @export
dmrff.annotate <- function(regions, annotation) {
    stopifnot(is.data.frame(regions))
    stopifnot(all(c("chr","start","end") %in% colnames(regions)))
    stopifnot(is.data.frame(annotation))
    stopifnot(all(c("chr","pos") %in% colnames(annotation)))

    regions <- regions[order(regions$chr, regions$start, regions$end),]
    annotation <- annotation[order(annotation$chr, annotation$pos),]
    
    sites.idx <- region.members(regions, annotation)
    regions.idx <- rep(1:length(sites.idx), sapply(sites.idx, length))
    cbind(regions[regions.idx,],
          annotation[unlist(sites.idx), setdiff(colnames(annotation), "chr")])
}
