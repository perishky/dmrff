impute.matrix <- function(x, FUN=function(x) rowMedians(x, na.rm=T)) {
    idx <- which(is.na(x), arr.ind=T)
    if (length(idx) > 0) {
        v <- FUN(x)
        v[which(is.na(v))] <- FUN(matrix(v, nrow=1))
        x[idx] <- v[idx[,"row"]]
    }
    x
}
