impute.matrix <- function(x, FUN=function(x) rowMedians(x, na.rm=T)) {
    ##idx <- which(is.na(x), arr.ind=T)
    ## the line above causes a 'long vectors not supported yet' error
    ## when x contains more elements the .Machine$integer.max
    ## the following 5 lines allows x to be much larger without generating this error
    idx <- lapply(1:ncol(x), function(i) {
        idx <- which(is.na(x[,i]))
        cbind(row=idx, col=rep(i,length(idx)))
    })
    idx <- do.call(rbind, idx)
    
    if (length(idx) > 0) {
        v <- FUN(x)
        v[which(is.na(v))] <- FUN(matrix(v, nrow=1))
        x[idx] <- v[idx[,"row"]]
    }
    x
}

