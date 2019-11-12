## estimate effect sizes from association tests of features
## se standard errors of the coefficients
## mat input data matrix (features x samples)
ivwfe.stats <- function(estimate, se, mat=NULL, rho=NULL) {
    ## http://onlinelibrary.wiley.com/doi/10.1002/sim.6835/full
    ## calculate rho
    if (is.null(rho)) {
        stopifnot(!is.null(mat))
        rho <- ivwfe.rho(mat)        
    }    
    ## remove missing values
    na <- is.na(estimate) | is.na(se)
    if (sum(na) > 0) {
        if (sum(na) == length(estimate))
            return(c(B=NA, S=NA))

        estimate <- estimate[!na]
        se <- se[!na]
        rho <- rho[!na,!na,drop=F]
    }
    ## calculate statistics
    ivwfe.stats0 <- function(estimate, se, rho) {
        omega <- (se%*%t(se))*rho
        omega.inv <- solve(omega)
        one <- matrix(1,nrow=nrow(omega.inv), ncol=1)
        S2 <- 1/(t(one) %*% omega.inv %*% one)
        c(B=S2 * (t(one) %*% omega.inv) %*% estimate,
          S=sqrt(S2))
    }
    return(ivwfe.stats0(estimate, se, rho))
}

ivwfe.getz <- function(estimate, se, mat=NULL, rho=NULL) {
    stats <- ivwfe.stats(estimate, se, mat, rho)
    as.vector(stats["B"]/stats["S"])
}



ivwfe.rho <- function(mat) {
    mat <- t(mat)
    ## correlation matrix + nudge (to ensure invertible)
    rho <- cor(mat, use="p") + diag(x=0.05,ncol(mat),ncol(mat))
}

ivwfe.ma <- function(estimates, se) {
    weights <- 1/se^2
    se <- sqrt(1/rowSums(weights, na.rm=T))
    estimates <- rowSums(estimates * weights, na.rm=T)/rowSums(weights, na.rm=T)
    z <- estimates/se
    p <- 2*pnorm(-abs(z), lower.tail=T)
    data.frame(estimate=estimates,
               se=se,
               z=z,
               p.value=p)
}

