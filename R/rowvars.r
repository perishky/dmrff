row.vars <- function(x, na.rm=F) {
    x <- x - rowMeans(x, na.rm=na.rm)
    if (na.rm)
        n <- rowSums(!is.na(x))
    else
        n <- ncol(x)
    rowSums(x^2,na.rm=na.rm)/(n-1)
    }
row.sds <- function(...) {
    sqrt(row.vars(...))
}
