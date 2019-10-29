## Code in this file from:
## 
## Jaffe AE, Murakami P, Lee H, Leek JT, Fallin DM, Feinberg AP, Irizarry
## RA (2012). Bump hunting to identify differentially methylated regions
## in epigenetic epidemiology studies. International journal of
## epidemiology, 41(1), 200?209. doi: 10.1093/ije/dyr238.
##
## https://github.com/rafalab/bumphunter

greaterOrEqual <- function (x, y) {
    precision <- sqrt(.Machine$double.eps)
    (x >= y) | (abs(x - y) <= precision)
}

## getSegments
bh.getSegments <- function(x, f = NULL, cutoff=quantile(abs(x), 0.99), assumeSorted = FALSE, verbose=FALSE){
    if(is.null(f))
        f <- rep(1L, length(x))
    stopifnot(length(x) == length(f))
    stopifnot(length(cutoff) <= 2)
    if(is.character(f))
        f <- as.factor(f)
    if(is.numeric(f))
        f <- as.integer(f)
    stopifnot(is.factor(f) || is.integer(f))
    if(length(cutoff) == 1)
        cutoff <- c(-cutoff, cutoff)
    cutoff <- sort(cutoff)

    reordered <- FALSE
    if(!assumeSorted && is.unsorted(f)) {
	od <- order(f)
	x <- x[od]
	f <- f[od]
	reordered <- TRUE
    }
        
    if(verbose) message("getSegments: segmenting")
    Indexes <- split(seq(along=x), f)
    direction <- as.integer(greaterOrEqual(x, cutoff[2]))
    direction[x <= cutoff[1]] <- -1L

    ## We now need to convert f into cid
    if(verbose) message("getSegments: splitting")
    segments <- cumsum(c(1, diff(direction) != 0)) +
        cumsum(c(1, diff(f) != 0))
    names(segments) <- NULL

    res <- list(upIndex = split(which(direction>0), segments[direction>0]),
                dnIndex = split(which(direction<0), segments[direction<0]),
                zeroIndex = split(which(direction==0), segments[direction==0]))
    names(res[[1]]) <- NULL
    names(res[[2]]) <- NULL
    names(res[[3]]) <- NULL

    if(reordered) {
        res <- lapply(res, function(sp) lapply(sp, function(xx) od[xx]))
    }
    res
}

## clusterMaker
bh.clusterMaker <- function(chr, pos, assumeSorted = FALSE, maxGap=300){
    nonaIndex <- which(!is.na(chr) & !is.na(pos))
    Indexes <- split(nonaIndex, chr[nonaIndex])
    clusterIDs <- rep(NA, length(chr))
    LAST <- 0
    for(i in seq(along = Indexes)){
        Index <- Indexes[[i]]
        x <- pos[Index]
        if(!assumeSorted){
            Index <- Index[order(x)]
            x <- pos[Index]
        }
        y <- as.numeric(diff(x) > maxGap)
        z <- cumsum(c(1, y))
        clusterIDs[Index] <- z + LAST
        LAST <- max(z) + LAST
    }
    clusterIDs
}

bh.regionFinder <- function(x, chr, pos, cluster=NULL, y=x, summary=mean,
                         ind=seq(along=x),order=TRUE, oneTable=TRUE,
                         maxGap=300, cutoff=quantile(abs(x), 0.99),
                         assumeSorted = FALSE, verbose = TRUE){
    if(any(is.na(x[ind]))){
        warning("NAs found and removed. ind changed.")
        ind <- intersect(which(!is.na(x)),ind)
    } 
    if(is.null(cluster))
        cluster <- bh.clusterMaker(chr, pos, maxGap=maxGap, assumeSorted = assumeSorted)
    Indexes <- bh.getSegments(x = x[ind], f = cluster[ind], cutoff = cutoff,
                           assumeSorted = assumeSorted, verbose = verbose)
    clusterN <- table(cluster)[as.character(cluster)]
    
    res <- vector("list",2)
    for(i in 1:2){
      res[[i]]<-
        data.frame(chr=sapply(Indexes[[i]],function(Index) chr[ind[Index[1]]]),
                   start=sapply(Indexes[[i]],function(Index) min(pos[ind[Index]])),
                   end=sapply(Indexes[[i]], function(Index) max(pos[ind[Index]])),
                   value=sapply(Indexes[[i]],function(Index)summary(y[ind[Index]])),
                   area=sapply(Indexes[[i]],function(Index)abs(sum(y[ind[Index]]))),
                   cluster=sapply(Indexes[[i]],function(Index)cluster[ind[Index]][1]),
                   indexStart=sapply(Indexes[[i]], function(Index) min(ind[Index])),
                   indexEnd = sapply(Indexes[[i]], function(Index) max(ind[Index])),
                   stringsAsFactors=FALSE)

      res[[i]]$L <- res[[i]]$indexEnd - res[[i]]$indexStart+1
      res[[i]]$clusterL <- sapply(Indexes[[i]], function(Index) clusterN[ind[Index]][1])
    }
    names(res) <- c("up","dn")
    if(order & !oneTable){
        if(nrow(res$up)>0) res$up <- res$up[order(-res$up$area),]
        if(nrow(res$dn)>0) res$dn <- res$dn[order(-res$dn$area),]
    }
    if(oneTable){
        res <- rbind(res$up,res$dn)
        if(order & nrow(res)>0) res <- res[order(-res$area),]
    }
    return(res)
}

