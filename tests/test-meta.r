library(dmrff)

options(mc.cores=4)

source("functions.r")

## generate datasets
set.seed(20180220)
n.datasets <- 10
n.sites <- 1000
n.samples <- 100
manifest <- generate.manifest(n.sites)
datasets <- lapply(1:n.datasets, function(i) {
    dataset <- generate.dataset(n.samples, manifest)                     
    var <- generate.true.dmr.var(dataset$methylation, manifest$chr, manifest$pos,
                                 cluster.sites=20, dmr.sites=10, cluster.position=0.5,
                                 r=0.5, maxgap=500, cluster.number=1)
    dataset$data$variable <- var$var
    dataset
})

## run EWAS and create meta-analysis objects for each
for (i in 1:length(datasets)) {
    datasets[[i]]$ewas <- with(datasets[[i]], {
        ewas(methylation, data, "variable")
    })
    datasets[[i]]$pre <- with(datasets[[i]], {
        dmrff.pre(estimate=ewas$estimate, se=ewas$se,
                  methylation=methylation, chr=manifest$chr, pos=manifest$pos)
    })
}

## make sure that 'sd' calculates correct
stopifnot(all(abs(datasets[[1]]$pre$sd - apply(datasets[[1]]$methylation, 1, sd)) < 2e-16))

## identify DMRs (there should be at least one)
ret <- dmrff.meta(lapply(datasets, function(dataset) dataset$pre),
                  maxgap=500, p.cutoff=0.05, verbose=T)

ret$dmrs[which(ret$dmrs$p.adjust < 0.05 & ret$dmrs$n > 1), ]
##   chr start  end n  B  S   estimate          se       z      p.value
## 1   1  6452 7133 8 NA NA 0.05147351 0.003374736 15.2526 1.582028e-52
##       p.adjust
## 1 1.882614e-49




ret$ewas$p.adjust <- p.adjust(ret$ewas$p.value, "bonferroni")
ret$ewas[with(ret$ewas, min(which(p.adjust < 0.05)):max(which(p.adjust < 0.05))),]
##      estimate          se         z      p.value chr  pos     p.adjust
## 29 0.03495611 0.005587286  6.256366 3.940522e-10   1 6280 3.940522e-07
## 30 0.05537181 0.005421784 10.212840 1.737056e-24   1 6452 1.737056e-21
## 31 0.05636577 0.005410459 10.417928 2.053781e-25   1 6485 2.053781e-22
## 32 0.05699328 0.005383369 10.586917 3.426934e-26   1 6627 3.426934e-23
## 33 0.06032290 0.005354420 11.266001 1.931436e-29   1 6630 1.931436e-26
## 34 0.06226195 0.005303403 11.740001 7.948510e-32   1 6683 7.948510e-29
## 35 0.04192692 0.005571865  7.524755 5.281939e-14   1 6991 5.281939e-11
## 36 0.04481686 0.005533165  8.099679 5.510442e-16   1 7044 5.510442e-13
## 37 0.04274037 0.005563133  7.682788 1.556625e-14   1 7133 1.556625e-11
## 38 0.02776729 0.005602469  4.956259 7.186349e-07   1 7602 7.186349e-04


##########################
## verify that meta-analysis gives the same results if
## inputs are not sorted by chromosomal position
udatasets <- lapply(datasets, function(dataset) {
    idx <- sample(1:nrow(manifest), nrow(manifest), replace=F)
    dataset <- list(data=dataset$data,
                    methylation=dataset$methylation[idx,],
                    manifest=manifest[idx,])
    dataset$ewas <- with(dataset, ewas(methylation, data, "variable"))
    dataset$pre <- with(dataset, dmrff.pre(estimate=ewas$estimate, se=ewas$se,
                                           methylation=methylation,
                                           chr=manifest$chr,
                                           pos=manifest$pos))
    dataset
})
                           
## identify DMRs 
uret <- dmrff.meta(lapply(udatasets, function(dataset) dataset$pre),
                  maxgap=500, p.cutoff=0.05, verbose=T)

uret$dmrs[which(uret$dmrs$p.adjust < 0.05 & uret$dmrs$n > 1), ]
##   chr start  end n  B  S   estimate          se       z      p.value
## 1   1  6452 7133 8 NA NA 0.05147351 0.003374736 15.2526 1.582028e-52
##       p.adjust
## 1 1.882614e-49


stopifnot(all(sapply(colnames(ret$dmrs), function(col)
                     identical(ret$dmrs[[col]], uret$dmrs[[col]]))))


###########################
## check meta-analysis
## by calculating everything 'by hand'

library(parallel)
## meta-analyse EWAS statistics
library(metafor)
ewas.ma.stats <- do.call(rbind, mclapply(1:n.sites, function(i) {
    tryCatch({
        fit <- rma.uni(yi=sapply(datasets, function(dataset) with(dataset$pre, estimate[i]/sd[i])),
                       sei=sapply(datasets, function(dataset) with(dataset$pre, se[i]/sd[i])),
                       method="FE")
        c(estimate=unname(fit$b["intrcpt",1]),
          se=fit$se,
          z=fit$zval,
          p.value=fit$pval)
    }, error=function(e) {
        c(estimate=0,se=1,z=0,p.value=1)
    })
}))
ewas.ma.stats <- as.data.frame(ewas.ma.stats)

## identify candidate regions from meta-analyzed statistics
candidates <- dmrff.candidates(ewas.ma.stats$estimate, ewas.ma.stats$p.value,
                               manifest$chr, manifest$pos)
## generate all candidate sub-regions
candidates <- do.call(rbind, lapply(1:nrow(candidates), function(i) {
    idx <- candidates$start.idx[i]:candidates$end.idx[i]
    candidates <- cbind(start.idx=idx, end.idx=idx)
    if (length(idx) > 1)
        candidates <- rbind(candidates, t(combn(idx, 2)))
    candidates
}))
candidates <- data.frame(chr=manifest$chr[candidates[,"start.idx"]],
                         start=manifest$pos[candidates[,"start.idx"]],
                         end=manifest$pos[candidates[,"end.idx"]])

## calculate DMR stats for each candidate region in each dataset
dmr.stats <- lapply(datasets, function(dataset)
                    dmrff.stats(candidates, dataset$pre$estimate, dataset$pre$se,
                                dataset$methylation, manifest$chr, manifest$pos))

## meta-analyse candidate region statistics
dmr.ma.stats <- do.call(rbind, lapply(1:nrow(candidates), function(i) {
    fit <- rma.uni(yi=sapply(dmr.stats, function(stats)stats$estimate[i]),
                   sei=sapply(dmr.stats, function(stats) stats$se[i]),
                   method="FE")
    c(estimate=unname(fit$b["intrcpt",1]),
      se=fit$se,
      z=fit$zval,
      p.value=fit$pval)
}))
dmr.ma.stats <- cbind(candidates, dmr.ma.stats)
dmr.ma.stats$p.adjust <- p.adjust(dmr.ma.stats$p.value, "bonferroni")

## cover all candidate regions with the 'most significant' sub-regions
dmr.ma.stats <- dmr.ma.stats[order(dmr.ma.stats$p.value),]
rownames(dmr.ma.stats) <- 1:nrow(dmr.ma.stats)
dmr.ma.stats$best <- T
for (i in 1:nrow(dmr.ma.stats)) {
    if (dmr.ma.stats$best[i]) {
        idx <- with(dmr.ma.stats,
                    which(chr == chr[i]
                          & (start[i] <= start & start <= end[i]
                             | start[i] <= start & end <= end[i]
                             | start <= start[i] & start[i] <= end
                             | start <= end[i] & end[i] <= end)))
        idx <- idx[which(idx > i)]
        if (length(idx) > 0)
            dmr.ma.stats$best[idx] <- F
    }
}

## the same dmrs should be identified by both methods
cols <- c("chr","start","end")
x <- dmr.ma.stats[which(dmr.ma.stats$best & dmr.ma.stats$p.adjust < 0.05),cols]
y <- ret$dmrs[which(ret$dmrs$p.adjust < 0.05),cols]
stopifnot(all(sort(apply(x, 1, paste))
              == sort(apply(y, 1, paste))))

## z-scores for the same regions should be identical
dmr.ma.stats <- dmr.ma.stats[match(with(ret$dmrs, paste(chr, start, end)),
                                   with(dmr.ma.stats, paste(chr, start, end))),]
stopifnot(diff(range(dmr.ma.stats$z - ret$dmrs$z)) < 1e-14)


###################################################
## test dmrff.cohort

dmrs1 <- with(datasets[[1]], dmrff(estimate=ewas$estimate,
                                   se=ewas$se,
                                   p.value=ewas$p.value,
                                   methylation=methylation,
                                   chr=pre$chr,
                                   pos=pre$pos))

dmrs2 <- dmrff.cohort(datasets[[1]]$pre)

## the exact set of candidates may not be exactly the same
## because p-values were calculated from estimate/se in 'pre'
dmrs1$id <- with(dmrs1, paste(chr, start, end))
dmrs2$id <- with(dmrs2, paste(chr, start, end))
common.candidates <- intersect(dmrs1$id, dmrs2$id)
length(common.candidates) ## 47
nrow(dmrs1) ## 47
nrow(dmrs2) ## 48
dmrs1 <- dmrs1[match(common.candidates, dmrs1$id),]
dmrs2 <- dmrs2[match(common.candidates, dmrs2$id),]
stopifnot(cor(dmrs1$estimate, dmrs2$estimate) >= 0.99)
stopifnot(cor(dmrs1$se, dmrs2$se) >= 0.99)
stopifnot(all((dmrs1$p.adjust < 0.05) == (dmrs2$p.adjust < 0.05)))
