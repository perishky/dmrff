---
title: "Identifying differentially methylated regions using dmrff"
output:
  pdf_document: default
  word_document:
    highlight: tango
  html_document:
    toc: true
---

# Identifying differentially methylated regions using dmrff

## Download and prepare an example DNA methylation dataset

We'll use a small cord blood DNA methylation dataset
available on GEO: [GSE79056](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79056).


```r
library(geograbi) ## https://github.com/yousefi138
```

```
## Loading required package: XML
```

```
## Loading required package: data.table
```

```r
samples <- geograbi.get.samples("GSE79056")
```

```
## Tue Dec 12 19:02:34 2023 reading
```

```r
vars <- geograbi.extract.characteristics(samples)
methylation <- geograbi.get.data("GSE79056") ## 86Mb
```

```
## Tue Dec 12 19:02:35 2023 reading
```

We'll be investigating associations with gestational age.

```r
colnames(vars)[which(colnames(vars) == "ga weeks")] <- "ga"
vars$ga <- as.numeric(vars$ga)
```

## Test DNA methylation associations 

Prepare surrogate variables to handle unknown confounding.

```r
library(sva, quietly=T)
```

```
## This is mgcv 1.8-41. For overview type 'help("mgcv-package")'.
```

```r
## construct EWAS model
mod <- model.matrix(~ gender + ga, vars)
## construct null model
mod0 <- mod[,1]
## to save time, SVA will be applied to a random selection of 5000 CpG sites
set.seed(20191029)
random.idx <- sample(1:nrow(methylation), 5000)
methylation.sva <- methylation[random.idx,]
## missing methylation values are replaced with mean values
methylation.mean <- rowMeans(methylation.sva, na.rm=T)
idx <- which(is.na(methylation.sva), arr.ind=T)
if (nrow(idx) > 0)
    methylation.sva[idx] <- methylation.mean[idx[,"row"]]
sva.fit <- sva(methylation.sva, mod=mod, mod0=mod0)
```

```
## Number of significant surrogate variables is:  3 
## Iteration (out of 5 ):1  2  3  4  5
```

Add surrogate variables to the model.

```r
design <- cbind(mod, sva.fit$sv)
```

Test the associations using `limma`. 

```r
library(limma)
fit <- lmFit(methylation, design)
fit <- eBayes(fit)
```

Save the summary statistics for gestational age.

```r
stats <- data.frame(estimate=fit$coefficients[,"ga"],
                    se=sqrt(fit$s2.post) * fit$stdev.unscaled[,"ga"],
                    p.value=fit$p.value[,"ga"])
```

## Add genomic coordinates for CpG sites (Illumina Beadchips)

`dmrff` and all other methods for identifying differentially methylated regions
require information about the genomic locations of all CpG sites. 

If your data was generated using one of the Illumina BeadChips, then it is possible
to annotate your EWAS summary statistics using the appropriate
[Bioconductor](https://www.bioconductor.org/) annotation package.

For the Illumina 450K microaray, we use the
`IlluminaHumanMethylation450kanno.ilmn12.hg19` package.


```r
if (!require(IlluminaHumanMethylation450kanno.ilmn12.hg19, quietly=T)) {
    install.packages("BiocManager")
    BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
}
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following object is masked from 'package:limma':
## 
##     plotMA
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unsplit, which.max, which.min
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following objects are masked from 'package:data.table':
## 
##     first, second
```

```
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
```

```
## 
## Attaching package: 'IRanges'
```

```
## The following object is masked from 'package:nlme':
## 
##     collapse
```

```
## The following object is masked from 'package:data.table':
## 
##     shift
```

```
## 
## Attaching package: 'matrixStats'
```

```
## The following objects are masked from 'package:genefilter':
## 
##     rowSds, rowVars
```

```
## 
## Attaching package: 'MatrixGenerics'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##     colWeightedMeans, colWeightedMedians, colWeightedSds,
##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##     rowWeightedSds, rowWeightedVars
```

```
## The following objects are masked from 'package:genefilter':
## 
##     rowSds, rowVars
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## 
## Attaching package: 'Biobase'
```

```
## The following object is masked from 'package:MatrixGenerics':
## 
##     rowMedians
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     anyMissing, rowMedians
```

```
## 
## Attaching package: 'Biostrings'
```

```
## The following object is masked from 'package:base':
## 
##     strsplit
```

```
## locfit 1.5-9.6 	 2022-07-11
```

```
## Setting options('download.file.method.GEOquery'='auto')
```

```
## Setting options('GEOquery.inmemory.gpl'=FALSE)
```

```r
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
```

If we were using the more recent Illumina MethylationEPIC BeadChip,
then we would load the `IlluminaHumanMethylationEPICanno.ilm10b2.hg19` package instead.

We construct an annotation data frame as follows:

```r
data(list="IlluminaHumanMethylation450kanno.ilmn12.hg19")
data(Locations)
data(Other)
annotation <- cbind(as.data.frame(Locations), as.data.frame(Other))
```

We then add the annotation to the `stats` object
(**Warning**: we assume that the row names of the `stats` object
are the Illumina CpG site identifiers, e.g. cg05775921).


```r
common <- intersect(rownames(methylation), rownames(annotation))
annotation <- annotation[match(common, rownames(annotation)),]
stats <- stats[match(common, rownames(stats)),]
methylation <- methylation[match(common, rownames(methylation)),]

stats <- cbind(stats, annotation)
```

## Apply `dmrff` to summary statistics

If `dmrff` is not already installed, we install it and then load it:


```r
if (!require(dmrff, quietly=T)) {
    library(devtools)
    install_github("perishky/dmrff")
}
library(dmrff)
```

Below we assume that you have already performed an epigenome-wide association analysis
and have loaded your summary statistics in R as a data frame `stats`
which has the following columns:
- `estimate` (regression coefficient),
- `se` (standard error of the coefficient),
- `p.value`,
- `chr` (chromosome of the CpG site),
- `pos` (position of the CpG site on the chromosome).

We also assume that you have loaded your DNA methylation dataset in R as matrix `methylation`
for which rows correspond to CpG sites and columns to samples.
The DNA methylation levels are necessary for `dmrff` to calculate
and adjust for dependencies between CpG sites.

`dmrff` can then be applied as follows.

```r
dmrs <- dmrff(estimate=stats$estimate,
              se=stats$se,
              p.value=stats$p.value,
              methylation=methylation,
              chr=stats$chr,
              pos=stats$pos,
              maxgap=500,
              verbose=T)
```

```
## [dmrff.candidates] Tue Dec 12 19:03:40 2023 Found  74588  candidate regions.
```

The algorithm will then identify differentially methylated regions by:

* Identifying regions composed of CpG sites with EWAS p-values < 0.05 and consist direction of effect (at most 500bp between CpG sites -- see `maxgap` parameter).
* Meta-analysing EWAS statistics within each region as well as sub-regions to identify the sub-region with lowest meta-analysed p-value.
* Adjusting meta-analysed p-values for multiple tests (all EWAS tests + sub-region meta-analysis tests).
* Finally, returning the DMR results.

Our output `dmrs` is a data frame listing all genomic regions tested
along with summary statistics for each.

We just keep regions with > 1 CpG site and Bonferroni adjusted p < 0.05.

```r
dmrs <- dmrs[which(dmrs$p.adjust < 0.05 & dmrs$n > 1),]
```
There are 1368 such DMRs.

Here are the first 10 regions:


```r
kable(dmrs[1:10,])
```



|   |chr   |    start|      end|  n|          B|         S|   estimate|        se|          z| p.value|  p.adjust|
|:--|:-----|--------:|--------:|--:|----------:|---------:|----------:|---------:|----------:|-------:|---------:|
|1  |chr6  | 32119616| 32121420| 39|  0.0594703| 0.0036198|  0.0594703| 0.0036198|  16.429200|       0| 0.0000000|
|3  |chr6  | 33245474| 33245779| 17| -0.1546681| 0.0117290| -0.1546681| 0.0117290| -13.186791|       0| 0.0000000|
|5  |chr6  | 33245893| 33246094|  7| -0.1404528| 0.0162813| -0.1404528| 0.0162813|  -8.626631|       0| 0.0000000|
|10 |chr6  | 31867847| 31869088| 32| -0.0591950| 0.0041669| -0.0591950| 0.0041669| -14.205932|       0| 0.0000000|
|12 |chr6  | 32036179| 32038958| 53|  0.0392814| 0.0050978|  0.0392814| 0.0050978|   7.705565|       0| 0.0000000|
|14 |chr11 | 14993378| 14993642|  4|  0.1616448| 0.0199695|  0.1616448| 0.0199695|   8.094569|       0| 0.0000000|
|15 |chr11 | 14994230| 14994606|  5|  0.1211192| 0.0196781|  0.1211192| 0.0196781|   6.155039|       0| 0.0004818|
|23 |chr6  | 30458158| 30459295| 11| -0.1580182| 0.0116069| -0.1580182| 0.0116069| -13.614168|       0| 0.0000000|
|24 |chr6  | 30459540| 30460600|  9| -0.0923506| 0.0112549| -0.0923506| 0.0112549|  -8.205369|       0| 0.0000000|
|27 |chr17 | 41278141| 41278425| 10| -0.1191146| 0.0189142| -0.1191146| 0.0189142|  -6.297623|       0| 0.0001940|

We can also list the information for each CpG site in these regions.


```r
sites <- dmrff.sites(dmrs, stats$chr, stats$pos)
sites <- cbind(sites, stats[sites$site, c("UCSC_RefGene_Name", "estimate", "p.value")])
sites <- cbind(sites, dmr=dmrs[sites$region,c("start","end","z","p.adjust")])
```

Here are the sites in the largest DMR:


```r
max.region <- which.max(dmrs$n)
kable(sites[which(sites$region == max.region),],row.names=F)
```



| region|   site|chr  |      pos|UCSC_RefGene_Name |   estimate|   p.value| dmr.start|  dmr.end|    dmr.z| dmr.p.adjust|
|------:|------:|:----|--------:|:-----------------|----------:|---------:|---------:|--------:|--------:|------------:|
|      5| 403116|chr6 | 32036179|TNXB              |  0.0010315| 0.0109635|  32036179| 32038958| 7.705565|            0|
|      5| 418198|chr6 | 32036258|TNXB              |  0.0000935| 0.8952558|  32036179| 32038958| 7.705565|            0|
|      5| 409882|chr6 | 32036278|TNXB              |  0.0005737| 0.4059512|  32036179| 32038958| 7.705565|            0|
|      5|  70316|chr6 | 32036356|TNXB              |  0.0076859| 0.0305838|  32036179| 32038958| 7.705565|            0|
|      5|  44051|chr6 | 32036404|TNXB              |  0.0007714| 0.2085628|  32036179| 32038958| 7.705565|            0|
|      5| 465690|chr6 | 32036449|TNXB              |  0.0009136| 0.0543140|  32036179| 32038958| 7.705565|            0|
|      5| 356956|chr6 | 32036530|TNXB              |  0.0008525| 0.1268730|  32036179| 32038958| 7.705565|            0|
|      5| 406883|chr6 | 32036532|TNXB              |  0.0012313| 0.0300202|  32036179| 32038958| 7.705565|            0|
|      5|  90566|chr6 | 32036611|TNXB              |  0.0006175| 0.2440241|  32036179| 32038958| 7.705565|            0|
|      5|  14425|chr6 | 32036677|TNXB              |  0.0024381| 0.2242297|  32036179| 32038958| 7.705565|            0|
|      5| 338599|chr6 | 32036709|TNXB              |  0.0000897| 0.8567340|  32036179| 32038958| 7.705565|            0|
|      5| 351435|chr6 | 32036724|TNXB              |  0.0000125| 0.9839411|  32036179| 32038958| 7.705565|            0|
|      5| 301883|chr6 | 32036743|TNXB              |  0.0006430| 0.2947373|  32036179| 32038958| 7.705565|            0|
|      5| 231708|chr6 | 32036822|TNXB              |  0.0037987| 0.0000617|  32036179| 32038958| 7.705565|            0|
|      5| 462157|chr6 | 32036848|TNXB              |  0.0008652| 0.2164398|  32036179| 32038958| 7.705565|            0|
|      5| 390942|chr6 | 32036897|TNXB              |  0.0010422| 0.1459314|  32036179| 32038958| 7.705565|            0|
|      5| 475183|chr6 | 32036911|TNXB              | -0.0001132| 0.8998336|  32036179| 32038958| 7.705565|            0|
|      5| 330755|chr6 | 32036954|TNXB              |  0.0028732| 0.0002115|  32036179| 32038958| 7.705565|            0|
|      5| 240567|chr6 | 32037018|TNXB              |  0.0024086| 0.0001065|  32036179| 32038958| 7.705565|            0|
|      5| 291118|chr6 | 32037097|TNXB              | -0.0003098| 0.9055158|  32036179| 32038958| 7.705565|            0|
|      5| 352756|chr6 | 32037148|TNXB              |  0.0009316| 0.1774843|  32036179| 32038958| 7.705565|            0|
|      5| 103583|chr6 | 32037290|TNXB              |  0.0020771| 0.0009626|  32036179| 32038958| 7.705565|            0|
|      5| 138987|chr6 | 32037339|TNXB              |  0.0012636| 0.0084907|  32036179| 32038958| 7.705565|            0|
|      5|  38523|chr6 | 32037342|TNXB              |  0.0017828| 0.0004507|  32036179| 32038958| 7.705565|            0|
|      5| 417928|chr6 | 32037421|TNXB              |  0.0006383| 0.2333117|  32036179| 32038958| 7.705565|            0|
|      5| 160017|chr6 | 32037426|TNXB              |  0.0009770| 0.1394028|  32036179| 32038958| 7.705565|            0|
|      5| 179334|chr6 | 32037453|TNXB              |  0.0011397| 0.0251615|  32036179| 32038958| 7.705565|            0|
|      5| 379817|chr6 | 32037474|TNXB              |  0.0004407| 0.4323696|  32036179| 32038958| 7.705565|            0|
|      5|  52129|chr6 | 32037587|TNXB              |  0.0009903| 0.1194331|  32036179| 32038958| 7.705565|            0|
|      5|  60243|chr6 | 32037602|TNXB              |  0.0004996| 0.4616771|  32036179| 32038958| 7.705565|            0|
|      5| 272610|chr6 | 32037632|TNXB              |  0.0002757| 0.6319254|  32036179| 32038958| 7.705565|            0|
|      5| 386420|chr6 | 32037653|TNXB              |  0.0010359| 0.0942346|  32036179| 32038958| 7.705565|            0|
|      5| 300711|chr6 | 32037735|TNXB              |  0.0015226| 0.0277306|  32036179| 32038958| 7.705565|            0|
|      5| 428514|chr6 | 32037751|TNXB              |  0.0003291| 0.5704414|  32036179| 32038958| 7.705565|            0|
|      5| 110400|chr6 | 32037800|TNXB              | -0.0004129| 0.7202410|  32036179| 32038958| 7.705565|            0|
|      5| 449994|chr6 | 32037847|TNXB              |  0.0028153| 0.0000569|  32036179| 32038958| 7.705565|            0|
|      5|  76776|chr6 | 32037913|TNXB              |  0.0006452| 0.1327238|  32036179| 32038958| 7.705565|            0|
|      5|  87741|chr6 | 32037916|TNXB              |  0.0014233| 0.0020986|  32036179| 32038958| 7.705565|            0|
|      5| 444167|chr6 | 32037978|TNXB              |  0.0012293| 0.0209276|  32036179| 32038958| 7.705565|            0|
|      5|   5081|chr6 | 32038011|TNXB              | -0.0010160| 0.2819558|  32036179| 32038958| 7.705565|            0|
|      5| 448813|chr6 | 32038027|TNXB              |  0.0017132| 0.0062513|  32036179| 32038958| 7.705565|            0|
|      5| 298918|chr6 | 32038045|TNXB              |  0.0009735| 0.3071656|  32036179| 32038958| 7.705565|            0|
|      5| 416617|chr6 | 32038097|TNXB              |  0.0008797| 0.1731802|  32036179| 32038958| 7.705565|            0|
|      5| 434939|chr6 | 32038155|TNXB              |  0.0018101| 0.0032325|  32036179| 32038958| 7.705565|            0|
|      5| 470733|chr6 | 32038177|TNXB              |  0.0023991| 0.0048218|  32036179| 32038958| 7.705565|            0|
|      5| 337001|chr6 | 32038185|TNXB              |  0.0025098| 0.0009209|  32036179| 32038958| 7.705565|            0|
|      5| 389831|chr6 | 32038230|TNXB              |  0.0016509| 0.0134805|  32036179| 32038958| 7.705565|            0|
|      5| 276348|chr6 | 32038462|TNXB              |  0.0010248| 0.0718441|  32036179| 32038958| 7.705565|            0|
|      5| 374532|chr6 | 32038541|TNXB              |  0.0027082| 0.0002163|  32036179| 32038958| 7.705565|            0|
|      5| 272153|chr6 | 32038747|TNXB              |  0.0028779| 0.0028012|  32036179| 32038958| 7.705565|            0|
|      5| 411763|chr6 | 32038881|TNXB              |  0.0020330| 0.0027048|  32036179| 32038958| 7.705565|            0|
|      5| 254265|chr6 | 32038887|TNXB              |  0.0006585| 0.3466269|  32036179| 32038958| 7.705565|            0|
|      5|  13100|chr6 | 32038958|TNXB              |  0.0026444| 0.0000209|  32036179| 32038958| 7.705565|            0|

Here we generate a simple plot of the largest DMR:


```r
dmrff.plot(
    dmrs$chr[max.region], dmrs$start[max.region], dmrs$end[max.region],
    stats$estimate, stats$se, stats$chr, stats$pos, ci=T, expand=1)
```

```
## Loading required package: ggplot2
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)
