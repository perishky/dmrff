library(dmrff)

options(mc.cores=4)

source("functions.r")

library(peakRAM)

stats <- read.csv(text="n.sites,n.samples,dat,peak,coverage
10000,100,NA,NA,NA
10000,200,NA,NA,NA
10000,400,NA,NA,NA
20000,800,NA,NA,NA
40000,1600,NA,NA,NA
40000,3200,NA,NA,NA")

for (i in 1:nrow(stats)) {
    cat(date(), stats$n.sites[i], stats$n.samples[i], "\n")
    ## construct a random dataset
    set.seed(20180220)
    n.sites <- stats$n.sites[i]
    n.samples <- stats$n.samples[i]
    manifest <- generate.manifest(n.sites)
    dataset <- generate.dataset(n.samples, manifest)
    ## generate a variable that has a dmr in the data
    var <- generate.true.dmr.var(
        dataset$methylation, manifest$chr, manifest$pos,
        cluster.sites=20, dmr.sites=10, cluster.position=0.5,
        r=0.5, maxgap=500)
    dataset$data$variable <- var$var
    ## run EWAS
    ewas.stats <- ewas(dataset$methylation, dataset$data, "variable")
    gc(reset=T)
    mem <- peakRAM({
        ## identify DMRs (there should be one)
        ret <- dmrff(
            ewas.stats$estimate, ewas.stats$se, ewas.stats$p.value,
            dataset$methylation,
            manifest$chr, manifest$pos, minmem=T)
    })
    stats$coverage[i] <- sum(ret$end.idx-ret$start.idx+1)
    stats$peak[i] <- mem$Peak_RAM_Used_MiB
    stats$dat[i] <- format(object.size(dataset), units="MB")
    cat(stats$peak[i], stats$dat[i], stats$coverage[i], "\n")
    rm(dataset)
    gc(reset=T)
}

stats
##   n.sites n.samples      dat  peak coverage
## 1   10000       100   7.8 Mb   3.5      547
## 2   10000       200  15.4 Mb   4.4      519
## 3   10000       400  30.7 Mb   7.9      545
## 4   20000       800 122.4 Mb  27.5     1066
## 5   40000      1600 488.9 Mb 102.9     2069
## 6   40000      3200 977.2 Mb 219.2     2225

## if minmem=F, then peak is about 11x larger!
##   n.sites n.samples      dat   peak coverage
## 1   10000       100   7.8 Mb   26.6      547
## 2   10000       200  15.4 Mb   48.6      519
## 3   10000       400  30.7 Mb   94.5      545
## 4   20000       800 122.4 Mb  248.2     1066
## 5   40000      1600 488.9 Mb 1223.6     2069
## 6   40000      3200 977.2 Mb 2444.6     2225


stats$n.size <- stats$n.sites * stats$n.samples

stats$dat <- as.numeric(sub(" Mb", "", stats$dat))

fit <- lm(peak ~ dat, data=stats)
coef(fit)
##(Intercept)         dat 
##  0.3682071   0.2211342 
cor(stats$peak, predict(fit, new=stats))
## [1] 0.9994415

fit <- lm(peak ~ n.size, data=stats)
coef(fit)
##  (Intercept)       n.size 
## 4.121876e-01 1.688032e-06 
cor(stats$peak, predict(fit, new=stats))
## [1] 0.9994488
predict(fit, new=data.frame(n.size=20000*850000))
## 28697 Mb used by dmrff

fit <- lm(dat ~ n.size, data=stats)
coef(fit)
##  (Intercept)       n.size 
## 2.008734e-01 7.633464e-06 
cor(stats$dat, predict(fit, new=stats))
predict(fit, new=data.frame(n.size=20000*850000))
## 129 769 Mb used by the dataset

## Conclusions:
## This trick of immediately reducing the methylation matrix
## to just sites in candidate regions reduces memory requirements
## to the memory required to load the methylation matrix.

        
