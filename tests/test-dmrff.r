library(dmrff)
library(parallel)

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
 ## 9.471601e-01  2.174437e-02  4.355887e+01 5.733639e-233 

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
##          estimate           se          t     p.value p.adjust
## 316 -0.0005470468 0.0013763417 -0.3974644 0.691926174        1
## 317  0.0040301904 0.0027890924  1.4449828 0.151788140        1
## 318  0.0085432014 0.0032910487  2.5958903 0.010947997        1
## 319  0.0091552368 0.0031276493  2.9271942 0.004289582        1
## 320  0.0074271524 0.0022735968  3.2666972 0.001519592        1
## 321  0.0074153650 0.0023855728  3.1084212 0.002488674        1
## 322  0.0084965433 0.0026289360  3.2319324 0.001695868        1
## 323  0.0032488870 0.0014932333  2.1757397 0.032082156        1
## 324  0.0017993648 0.0008530498  2.1093316 0.037573678        1
## 325  0.0065404556 0.0019777913  3.3069492 0.001336953        1

## identify DMRs (there should be one)
ret <- dmrff(ewas.stats$estimate, ewas.stats$se, ewas.stats$p.value,
             dataset$methylation,
             manifest$chr, manifest$pos)
stopifnot(sum(ret$p.adjust < 0.05) == 1)

ret[which(ret$p.adjust < 0.05),]
##   chr start   end n start.idx end.idx start.orig end.orig   z.orig     p.orig
## 1   1 63431 64019 5       318     322        318      325 2.322495 0.02020632
##          z      p.value    p.adjust
## 1 4.568012 4.923719e-06 0.005307769

