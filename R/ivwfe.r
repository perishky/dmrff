## coef effect sizes from association tests of features
## se standard errors of the coefficients
## mat input data matrix (features x samples)
ivwfe.stats <- function(coef, se, mat=NULL, rho=NULL) {
    ## From James Staley:
    ##     Here is the link to Stephen Burgess’ paper we discussed on
    ## Thursday last week:
    ## http://onlinelibrary.wiley.com/doi/10.1002/sim.6835/full    
    ##     I have taken a look at the maths and Stouffer’s test statistic
    ## weighted by 1/SE is the same as IVW FE meta-analysis beta/SE. So, the
    ## generalised linear model approach in Stephen’s paper which is based on
    ## IVW FE meta-analysis would probably be very similar to your Stouffer’s
    ## method corrected for correlated Z-statistics.    
    ##     The test statistic would be: T = B/S ~ N(0,1)    
    ## where B = (1^TΩ^-11)^-11^TΩ^-1β and S = sqrt((1^TΩ^-11)^-1) where β are
    ## the effect estimates, Ω is the variance-covariance matrix of the CpGs
    ## and 1 is a vector of 1’s the same length as the number of CpGs.

    if (is.null(rho)) {
        stopifnot(!is.null(mat))
        rho <- ivwfe.rho(mat)        
    }
    ## The second diagonal matrix is the 'nudge' matrix to ensure matrix inversion

    ## remove missing values
    na <- is.na(coef) | is.na(se)
    if (sum(na) > 0) {
        if (sum(na) == length(coef))
            return(c(B=NA, S=NA))

        coef <- coef[!na]
        se <- se[!na]
        rho <- rho[!na,!na,drop=F]
    }
    
    omega <- (se%*%t(se))*rho
    omega.inv <- solve(omega)
    one <- matrix(1,nrow=nrow(omega.inv), ncol=1)
    S2 <- 1/(t(one) %*% omega.inv %*% one)
    c(B=S2 * (t(one) %*% omega.inv) %*% coef,
      S=sqrt(S2))    
}

ivwfe.getz <- function(coef, se, mat=NULL, rho=NULL) {
    stats <- ivwfe.stats(coef, se, mat, rho)
    as.vector(stats["B"]/stats["S"])
}



ivwfe.rho <- function(mat) {
    mat <- t(mat)
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

