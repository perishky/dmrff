## for (file in list.files("../R", ".r$", full.names=T))
##     source(file)

library(dmrff)
library(parallel)

options(mc.cores=4)

source("functions.r")
source("bumphunter.r") ## just for the data simulation

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
                           
## identify DMRs (there should be one)
ret <- drmff.meta(lapply(datasets, function(dataset) dataset$pre),
                  maxgap=500, p.cutoff=0.05, verbose=T)
ret$dmrff[which(ret$dmrff$p.adjust < 0.05), ]
## > ret$dmrff[which(ret$dmrff$p.adjust < 0.05), ]
##   chr start  end n start.idx end.idx start.orig end.orig   z.orig       p.orig
## 1   1  6485 6627 2        31      32         28       38 4.652635 3.277197e-06
## 2   1  6452 6452 1        30      30         28       38 4.652635 3.277197e-06
## 3   1  6280 6280 1        29      29         28       38 4.652635 3.277197e-06
## 5   1  6683 7044 3        34      36         28       38 4.652635 3.277197e-06
## 6   1  6630 6630 1        33      33         28       38 4.652635 3.277197e-06
## 7   1  7133 7602 2        37      38         28       38 4.652635 3.277197e-06
##           z      p.value     p.adjust
## 1 11.159733 6.418373e-29 7.432475e-26
## 2  9.292509 1.506934e-20 1.745029e-17
## 3  6.173627 6.674092e-10 7.728598e-07
## 5 10.323678 5.507193e-25 6.377329e-22
## 6 10.295473 7.385555e-25 8.552473e-22
## 7  4.862225 1.160736e-06 1.344133e-03

ret$ewas$p.adjust <- p.adjust(ret$ewas$p.value, "bonferroni")
ret$ewas[with(ret$ewas, min(which(p.adjust < 0.05)):max(which(p.adjust < 0.05))),]
##       estimate           se         z      p.value chr  pos     p.adjust
## 29 0.004265129 0.0006742131  6.326085 2.514591e-10   1 6280 2.514591e-07
## 30 0.006215383 0.0006527400  9.521988 1.698969e-21   1 6452 1.698969e-18
## 31 0.005596188 0.0005731201  9.764424 1.600180e-22   1 6485 1.600180e-19
## 32 0.006125614 0.0006093820 10.052175 8.986206e-24   1 6627 8.986206e-21
## 33 0.005464965 0.0005180199 10.549721 5.094832e-26   1 6630 5.094832e-23
## 34 0.003130973 0.0003476549  9.005979 2.137506e-19   1 6683 2.137506e-16
## 35 0.003438597 0.0005115889  6.721408 1.799771e-11   1 6991 1.799771e-08
## 36 0.003072825 0.0004192734  7.328930 2.319973e-13   1 7044 2.319973e-10

########################### check meta-analysis
########################### by calculating everything 'by hand'

## meta-analyse EWAS statistics
library(metafor)
ewas.ma.stats <- do.call(rbind, mclapply(1:n.sites, function(i) {
    fit <- rma.uni(yi=sapply(datasets, function(dataset) dataset$pre$estimate[i]),
                   sei=sapply(datasets, function(dataset) dataset$pre$se[i]),
                   method="FE")
    c(estimate=unname(fit$b["intrcpt",1]),
      se=fit$se,
      z=fit$zval,
      p.value=fit$pval)
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
x <- dmr.ma.stats[which(dmr.ma.stats$best & dmr.ma.stats$p.adjust < 0.05),
             c("chr","start","end","z")]
y <- ret$dmrff[which(ret$dmrff$p.adjust < 0.05),
               c("chr","start","end","z")]
stopifnot(all(sort(apply(x, 1, paste)) == sort(apply(y, 1, paste))))

## z-scores for the same regions should be identical
dmr.ma.stats <- dmr.ma.stats[match(with(ret$dmrff, paste(chr, start, end)),
                                   with(dmr.ma.stats, paste(chr, start, end))),]
stopifnot(diff(range(dmr.ma.stats$z - ret$dmrff$z)) < 1e-14)

