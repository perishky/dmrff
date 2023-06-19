library(dmrff)

options(mc.cores=4)

source("functions.r")

## construct a random dataset
set.seed(20180220)
n.sites <- 1000
n.samples <- 100
manifest <- generate.manifest(n.sites)
dataset <- generate.dataset(n.samples, manifest)

## show methylation correlation structure
r <- sapply(2:nrow(dataset$methylation), function(i)
            cor(dataset$methylation[i,], dataset$methylation[i-1,]))
distance <- tail(manifest$pos,-1) - head(manifest$pos,-1)
distance <- 1/exp(distance/200)
coef(summary(lm(r~distance)))["distance",]
##     Estimate    Std. Error       t value      Pr(>|t|) 
## 9.580629e-01  1.950377e-02  4.912193e+01 1.783062e-268 

## run EWAS
ewas.stats <- ewas(dataset$methylation, dataset$data, "variable")

## identify DMRs (there should be none)
ret <- dmrff(ewas.stats$estimate, ewas.stats$se, ewas.stats$p.value,
             dataset$methylation,
             manifest$chr, manifest$pos)
stopifnot(sum(ret$p.adjust < 0.05) == 0)

## generate a variable that has a dmr in the data
var <- generate.true.dmr.var(dataset$methylation, manifest$chr, manifest$pos,
                             cluster.sites=20, dmr.sites=10, cluster.position=0.5,
                             r=0.5, maxgap=500)
dataset$data$variable <- var$var

## run EWAS
ewas.stats <- ewas(dataset$methylation, dataset$data, "variable")

## look at the stats for the DMR CpG sites
ewas.stats[var$dmr.idx,]
##        estimate           se        t      p.value   p.adjust
## 594 0.004540030 0.0023754879 1.911199 5.902572e-02 1.00000000
## 595 0.010686986 0.0027294084 3.915495 1.708092e-04 0.17080922
## 596 0.009171777 0.0023424307 3.915495 1.708092e-04 0.17080922
## 597 0.007929622 0.0017906810 4.428271 2.562602e-05 0.02562602
## 598 0.003468908 0.0010697681 3.242673 1.639467e-03 1.00000000
## 599 0.001634432 0.0008076755 2.023624 4.584787e-02 1.00000000
## 600 0.003870655 0.0012709909 3.045384 3.015081e-03 1.00000000
## 601 0.006311263 0.0016498583 3.825337 2.349228e-04 0.23492276
## 602 0.008159530 0.0024571669 3.320707 1.279400e-03 1.00000000
## 603 0.008063059 0.0024281156 3.320707 1.279400e-03 1.00000000

## identify DMRs (there should be one)
ret <- dmrff(ewas.stats$estimate, ewas.stats$se, ewas.stats$p.value,
             dataset$methylation,
             manifest$chr, manifest$pos)
stopifnot(sum(ret$p.adjust < 0.05) == 1)

ret[which(ret$p.adjust < 0.05),]
##   chr  start    end n start.idx end.idx start.orig end.orig   z.orig
## 1   1 122502 122796 7       595     601        595      604 5.025078
##         p.orig        B          S estimate         se        z      p.value
## 1 5.032271e-07 0.054384 0.01079676 0.054384 0.01079676 5.037065 4.727241e-07
##      p.adjust
## 1 0.000510542

## test unordered input
idx <- sample(1:nrow(ewas.stats), size=nrow(ewas.stats), replace=F)
ret <- dmrff(ewas.stats$estimate[idx],
             ewas.stats$se[idx],
             ewas.stats$p.value[idx],
             dataset$methylation[idx,,drop=F],
             manifest$chr[idx],
             manifest$pos[idx])

stopifnot(sum(ret$p.adjust < 0.05) == 1)

ret[which(ret$p.adjust < 0.05),]
##   chr  start    end n        B          S estimate         se        z
## 1   1 122502 122796 7 0.054384 0.01079676 0.054384 0.01079676 5.037065
##        p.value    p.adjust
## 1 4.727241e-07 0.000510542

