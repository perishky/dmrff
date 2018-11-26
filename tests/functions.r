## create variable y, such that cor(x,y) == r
rcor <- function(x,r) {
    e <- rnorm(length(x), mean=0, sd=sqrt(1-r^2))
    mx <- mean(x)
    sx <- sd(x)
    s <- (x-mx)/sx
    (r*s + e)*sx + mx
}

## create a variable that corresponds to a true dmr in a dataset
## 1. select a random cluster of CpG sites of minimum size
##    (defined by arguments cluster.sites, maxgap and cluster.number)
## 2. select CpG sites for DMR in the cluster (defined by dmr.sites,cluster.position)
## 3. sum standardized methylation for CpG sites in the DMR
## 4. generate variable with given correlation with the sum (defined by r).
generate.true.dmr.var <- function(mat, chr, pos, cluster.sites=40, dmr.sites=20, cluster.position=0.5, r=0.5, maxgap=500, cluster.number=NULL) {
    stopifnot(dmr.sites <= cluster.sites)
    stopifnot(cluster.position >= 0 && cluster.position <= 1)
    stopifnot(r >= 0 && r <= 1)
    ## identify clusters
    clusters <- bh.clusterMaker(chr, pos, maxGap=maxgap)
    cluster.size <- table(clusters)
    cluster.size <- cluster.size[which(cluster.size >= cluster.sites)]

    ## select a cluster
    if (is.null(cluster.number))
        cluster <- sample(names(cluster.size), 1)
    else {
        if (cluster.number > length(cluster.size))
            cluster.number <- 1
        cluster <- names(cluster.size)[cluster.number]
    }
        
    ## indices of cpg sites in the cluster
    cluster.idx <- which(clusters==cluster)

    ## indices of cpg sites in the dmr within in the cluster
    dmr.center <- max(1, floor(length(cluster.idx)*cluster.position))
    dmr.start <- max(1, dmr.center - floor(dmr.sites/2))
    dmr.end <- min(length(cluster.idx), dmr.start + dmr.sites -1)
    dmr.idx <- cluster.idx[dmr.start:dmr.end]

    ## sum CpG sites after scaling methylation
    dmr.mat <- t(scale(t(mat[dmr.idx,,drop=F])))
    dmr.avg <- colSums(dmr.mat)

    ## create a random variable correlated with the sum
    var <- rcor(dmr.avg,r)
    
    list(var=var, dmr.idx=dmr.idx, cluster.idx=cluster.idx)
}

## generate a random variable correlated with x
## the correlation is determined by 'distance'
## with correlation decreasing exponentially with increasing distance
generate.spatial <- function(x, distance) {
    mean.cor <- function(distance) 1/exp(distance/200)
    r <- rnorm(1, mean=mean.cor(distance), sd=0.15)
    if (r < -1) r <- -1
    if (r > 1) r <- 1
    rcor(x, r)
}

## generate DNA methylation for dataset of 'n' individuals
## with CpG sites at given positions ('pos') on a chromosome
## CpG sites are correlated with neighboring CpG sites according
## to their distance from the neighbors (see generate.spatial() above)
generate.methylation <- function(n, pos) {
    cpg.mean <- rep(NA, length(pos))
    cpg.mean[1] <- runif(1)
    cpg.change <- rnorm(length(pos), sd=0.1)
    for (i in 2:length(cpg.mean)) {
        new <- cpg.mean[i-1] + cpg.change[i]
        if (new < 0 | new > 1) cpg.change[i] <- -cpg.change[i]
        cpg.mean[i] <- cpg.mean[i-1] + cpg.change[i]
    }
    
    cpg.sd <- rnorm(length(pos), mean=0.2, sd=0.025)
    
    methylation <- matrix(NA, ncol=n, nrow=length(pos))
    methylation[1,] <- rnorm(ncol(methylation))
    distance <- tail(pos,-1) - head(pos,-1)
    for (i in 2:nrow(methylation))
        methylation[i,] <- generate.spatial(methylation[i-1,], distance[i-1])

    methylation <- t(scale(t(methylation)))
    methylation <- methylation*cpg.sd + cpg.mean
    idx <- which(methylation < 0 | methylation > 1, arr.ind=T)[,1]
    for (i in unique(idx)) {
        med.cpg <- median(methylation[i,])
        min.cpg <- min(methylation[i,])
        max.cpg <- max(methylation[i,])
        factor.min <- med.cpg/(med.cpg - min.cpg)
        factor.max <- (1-med.cpg)/(max.cpg - med.cpg)
        factor <- 1
        if (min.cpg < 0) factor <- min(factor.min, factor)
        if (max.cpg > 1) factor <- min(factor.max, factor)        
        methylation[i,] <- (methylation[i,] - med.cpg)*factor + med.cpg
    }
    methylation
}

## generate random positions for 'n' CpG sites across a chromosome
## on average consecutive CpG sites are 200bp apart
generate.manifest <- function(n=1000) {
    manifest <- data.frame(chr=1, pos=sample(1:(200*n), n))
    manifest <- manifest[order(manifest$chr, manifest$pos),]
    manifest
}

## generate a dataset with methylation data and variables for 'n' samples
## with CpG site positions given by 'manifest'
generate.dataset <- function(n,manifest) {
    ## variable of interest and covariates
    data <- data.frame(variable=rnorm(n), ## variable of interest
                       covariate=rnorm(n),             
                       categorical=factor(sample(0:3,n,replace=T)))

    methylation <- generate.methylation(n, manifest$pos)

    list(data=data, methylation=methylation)
}

## test associations of 'data[,varname]' with each row of 'methylation'
## adjusting for other covariates in 'data'.
ewas <- function(methylation, data, varname) {
    ret <- do.call(rbind, lapply(1:nrow(methylation), function(i) {
        fit <- lm(methylation[i,] ~ ., data=data)
        coef(summary(fit))[varname,]
    }))
    ret <- as.data.frame(ret)
    colnames(ret) <- c("estimate", "se", "t", "p.value")
    ret$p.adjust <- p.adjust(ret$p.value, "bonferroni")
    ret
}



