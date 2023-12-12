#' Plot a region
#'
#' dmrff.plot
#'
#' Calculate statistics for a set of genomic regions.
#'
#' @param region.chr/start/end Genomic region to plot.
#' @param estimate Vector of EWAS effect estimates
#' (corresponds to rows of \code{methylation}).
#' @param se Vector of standard errors of the coefficients.
#' @param methylation Methylation matrix (rows=features, columns=samples).
#' @param chr Feature chromosome (corresponds to rows of \code{methylation}).
#' @param pos Feature chromosome position.
#' @param expand Proportion to expand to include sites beyond the beginning and end of the region (default: 1, i.e. plot 3x the length of the region).
#' @param ci Show confidence intervals rather than standard errors (default: TRUE).
#' @param verbose If \code{TRUE} (default), then output status messages.
#' @return \code{\link{ggplot}} object showing the plot. 
#'
#' @export
dmrff.plot <- function(region.chr, region.start, region.end, estimate, se, chr, pos, expand=1, ci=T, verbose=T) {
    require(ggplot2)
    stopifnot(length(region.chr)==1)
    stopifnot(length(region.start)==1)
    stopifnot(length(region.end)==1)
    stopifnot(is.numeric(region.start))
    stopifnot(is.numeric(region.end))
    stopifnot(region.start < region.end)
    stopifnot(is.vector(estimate))
    stopifnot(is.vector(se))
    stopifnot(is.vector(chr))
    stopifnot(is.vector(pos))
    stopifnot(length(estimate)==length(se))
    stopifnot(length(estimate)==length(chr))
    stopifnot(length(estimate)==length(pos))
    stopifnot(region.chr %in% chr)
    region.len <- region.end-region.start
    plot.start <- max(0,region.start - region.len*expand)
    plot.end <- region.end + region.len*expand
    idx <- which(chr==region.chr & pos >= plot.start & pos <= plot.end)
    if (length(idx) == 0)
        stop("The region does not contain any CpG sites")
    dat <- data.frame(estimate,se,pos)[idx,]
    dat$ymax <- dat$estimate + dat$se*ifelse(ci,1.96,1)
    dat$ymin <- dat$estimate - dat$se*ifelse(ci,1.96,1)
    dat <- dat[order(dat$pos),]
    gap <- max(1, min(tail(dat$pos,-1)-head(dat$pos,-1), na.rm=T)/2)
    (ggplot(dat, aes(x=pos, y=estimate, ymin=ymin, ymax=ymax)) +
     geom_rect(
         aes(xmin=region.start-gap, xmax=region.end+gap, ymin=-Inf, ymax=Inf),
         fill="#BBBBBB", alpha=0.1) +
     geom_hline(yintercept=0, linetype="dashed", color="black") +
     geom_pointrange() +
     labs(
         y=paste0("estimate (", ifelse(ci, "95% CI", "SE"), ")"),
         x=paste0("chromosome ", region.chr, " position")))    
}



#' Remove points from a scatter plot where density is really high
#' @param x x-coordinates vector
#' @param y y-coordinates vector
#' @param resolution number of partitions for the x and y-dimensions.
#' @param max.per.cell maximum number of points per x-y partition.
#' @return index into the points that omits points from x-y partitions
#' so that each has at most \code{max.per.cell} points.
scatter.thinning <- function(x,y,resolution=100,max.per.cell=100) {
    x.cell <- floor((resolution-1)*(x - min(x,na.rm=T))/diff(range(x,na.rm=T))) + 1
    y.cell <- floor((resolution-1)*(y - min(y,na.rm=T))/diff(range(y,na.rm=T))) + 1
    z.cell <- x.cell * resolution + y.cell
    frequency.table <- table(z.cell)
    frequency <- rep(0,max(z.cell, na.rm=T))
    frequency[as.integer(names(frequency.table))] <- frequency.table
    f.cell <- frequency[z.cell]
    
    big.cells <- length(which(frequency > max.per.cell))
    sort(c(which(f.cell <= max.per.cell),
           sample(which(f.cell > max.per.cell),
                  size=big.cells * max.per.cell, replace=F)),
         decreasing=F)
}

#' Manhattan plot
#'
#' @param object Output from the \code{\link{dmrff}()}, \code{\link{dmrff.meta}()}
#' or \code{\link{dmrff.cohort}} functions.
#' @param title Title for the plot (Default: "Manhattan plot").
#' @return \code{\link{ggplot}} object showing the Manhattan plot. 
#' @export
dmrff.manhattan.plot <- function(object, chr, pos, title="Manhattan plot") {
    require(ggplot2)
    stopifnot(class(object) == "data.frame")
    stopifnot(all(c("chr","start","end","p.value","n") %in% colnames(object)))
    chr <- as.character(chr)
    stopifnot(all(object[["chr"]] %in% chr))
    stopifnot(length(chr) == length(pos))
    object$loc <- paste(as.character(object$chr), object$start)
    loc <- paste(chr,pos)
    stopifnot(all(object$loc %in% loc))
    
    chromosomes <- unique(chr)
    chromosomes <- chromosomes[order(as.integer(sub("chr","",tolower(chromosomes))))]

    stats <- data.frame(chr=factor(as.character(chr), levels=chromosomes),
                        pos=pos,
                        chr.colour=0)
    stats$loc <- loc
    stats$chr.colour[stats$chromosome %in% chromosomes[seq(1,length(chromosomes),2)]] <- 1
    object$p.value[which(object$p.value < .Machine$double.xmin)] <- .Machine$double.xmin
    object$logp <- -log(object$p.value,10) * sign(object$estimate)
    object$pos <- object$start
    stats$stat <- 0
    stats$stat[match(object$loc, stats$loc)] <- object$logp
        
    stats <- stats[order(stats$stat, decreasing=T),]
    
    chr.lengths <- sapply(chromosomes, function(chr) max(stats$pos[which(stats$chr == chr)]))
    chr.lengths <- as.numeric(chr.lengths)
    chr.starts <- c(1,cumsum(chr.lengths)+1)
    names(chr.starts) <- c(chromosomes, "NA")
    stats$global <- stats$pos + chr.starts[stats$chr] - 1

    stats <- stats[which(abs(stats$stat) > -log(0.05,10)),]
    
    selection.idx <- scatter.thinning(stats$global, stats$stat,
                                      resolution=100, max.per.cell=100)
    stats <- stats[selection.idx,]
    
    ret <- (ggplot(stats, aes(x=pos, y=stat)) +
            geom_point(aes(colour=chr.colour)) +
            facet_grid(. ~ chr, space="free_x", scales="free_x") +
            theme(strip.text.x = element_text(angle = 90)) +
            guides(colour=FALSE) +
            labs(x="Position",
                 y=bquote(-log[10]("p-value") * sign(beta))) +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
            ggtitle(title))

    if (any(object$p.adjust < 1)) {
        sig.threshold <- 0.05/max(object$p.adjust/object$p.value, na.rm=T)
        ret <- (ret +
                geom_hline(yintercept=log(sig.threshold,10), colour="red") +
                geom_hline(yintercept=-log(sig.threshold,10), colour="red"))
    }

    ret
}


