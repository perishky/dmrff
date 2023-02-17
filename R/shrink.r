shrink.candidates <- function(start.idx, end.idx, FUN, ...) {
    stopifnot(length(start.idx) == length(end.idx))
    stopifnot(all(start.idx <= end.idx))
    do.call(rbind, parallel::mclapply(1:length(start.idx), function(i) {
        shrink.candidate(start.idx[i], end.idx[i], FUN, ...)
    }))    
}

shrink.candidate <- function(start.idx, end.idx, FUN, ...) {
    n <- end.idx-start.idx + 1
    z <- matrix(0, ncol=n, nrow=n)
    for (i in 1:n) 
        for (j in i:n) 
            z[i,j] <- FUN(start.idx + i - 1, start.idx + j - 1, ...)

    shrink.fun <- function(z, start, end) {
        idx <- start:end
        max.idx <- which(abs(z[idx,idx,drop=F]) >= max(abs(z[idx,idx]), na.rm=T), arr.ind=T)
        if (length(max.idx) == 0) {
            return(NULL)
        }
        best.start <- idx[max.idx[1,1]]
        best.end <- idx[max.idx[1,2]]
        ret <- cbind(start.idx=best.start, end.idx=best.end)
        if (best.start > start)
            ret <- rbind(ret, shrink.fun(z, start, best.start-1))
        if (best.end < end)
            ret <- rbind(ret, shrink.fun(z, best.end+1, end))
        ret
    }

    ret <- shrink.fun(z, 1, n)
    tryCatch({
        data.frame(start.idx=start.idx + ret[,1] - 1,
                end.idx=start.idx + ret[,2] - 1,
                z=z[ret],
                start.orig=start.idx,
                end.orig=end.idx,
                z.orig=z[1,n])
    }, error=function(e) {
        save(z,n,start.idx,end.idx,ret,file="shrink-error-20230217.rda")
        stop(e)
    })
}

calculate.number.shrink.tests <- function(stats) {
    originals <- unique(stats[,c("start.orig","end.orig")])    
    originals$tests <- pmax(choose(originals$end.orig-originals$start.orig+1, 2),1)
    sum(originals$tests)
}

