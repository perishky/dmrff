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

```r
samples <- geograbi.get.samples("GSE79056")
```

```r
vars <- geograbi.extract.characteristics(samples)
methylation <- geograbi.get.data("GSE79056") ## 86Mb
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
## [dmrff.candidates] Tue Oct 29 17:01:56 2019 Found  41834  candidate regions.
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
There are 440 such DMRs.

Here are the first 10 regions:


```r
kable(dmrs[1:10,])
```



|    |chr   |    start|      end|  n|          B|         S|   estimate|        se|         z| p.value|  p.adjust|
|:---|:-----|--------:|--------:|--:|----------:|---------:|----------:|---------:|---------:|-------:|---------:|
|1   |chr6  | 32120584| 32120878|  6|  0.0133915| 0.0014515|  0.0133915| 0.0014515|  9.225734|       0| 0.0000000|
|2   |chr6  | 32119944| 32120203|  3|  0.0047073| 0.0006831|  0.0047073| 0.0006831|  6.891180|       0| 0.0000021|
|4   |chr6  | 32121143| 32121156|  2|  0.0089995| 0.0012446|  0.0089995| 0.0012446|  7.230775|       0| 0.0000002|
|21  |chr6  | 32013699| 32017224| 47| -0.0004568| 0.0000769| -0.0004568| 0.0000769| -5.943456|       0| 0.0010452|
|23  |chr6  | 33245488| 33245537|  3| -0.0065765| 0.0007459| -0.0065765| 0.0007459| -8.816650|       0| 0.0000000|
|28  |chr6  | 33245893| 33246094|  5| -0.0038212| 0.0005975| -0.0038212| 0.0005975| -6.395753|       0| 0.0000598|
|46  |chr11 | 14994230| 14994561|  4|  0.0052970| 0.0008988|  0.0052970| 0.0008988|  5.893703|       0| 0.0014143|
|76  |chr6  | 32117292| 32117377|  2|  0.0030610| 0.0005236|  0.0030610| 0.0005236|  5.846178|       0| 0.0018838|
|98  |chr6  | 30653191| 30653242|  2| -0.0032880| 0.0004738| -0.0032880| 0.0004738| -6.940174|       0| 0.0000015|
|103 |chr6  | 30653407| 30653732|  4| -0.0033657| 0.0006041| -0.0033657| 0.0006041| -5.571241|       0| 0.0094726|

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



| region|   site|chr  |      pos|UCSC_RefGene_Name |   estimate|   p.value| dmr.start|  dmr.end|     dmr.z| dmr.p.adjust|
|------:|------:|:----|--------:|:-----------------|----------:|---------:|---------:|--------:|---------:|------------:|
|      4| 206488|chr6 | 32013699|TNXB;TNXB;TNXB    |  0.0033248| 0.0013753|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 259872|chr6 | 32013710|TNXB;TNXB;TNXB    |  0.0039750| 0.0001142|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 253108|chr6 | 32013974|TNXB;TNXB         |  0.0007077| 0.1580097|  32013699| 32017224| -5.943456|    0.0010452|
|      4|  76731|chr6 | 32014003|TNXB;TNXB         |  0.0014911| 0.0152387|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 132640|chr6 | 32014148|TNXB;TNXB         |  0.0016923| 0.0819380|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 300117|chr6 | 32014189|TNXB;TNXB         |  0.0020784| 0.0073160|  32013699| 32017224| -5.943456|    0.0010452|
|      4|  83714|chr6 | 32014300|TNXB;TNXB         |  0.0016903| 0.0294514|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 246136|chr6 | 32014476|TNXB;TNXB         |  0.0013224| 0.5403594|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 271379|chr6 | 32014484|TNXB;TNXB         |  0.0007101| 0.5998312|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 145379|chr6 | 32014510|TNXB;TNXB         |  0.0016457| 0.2082371|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 174594|chr6 | 32014605|TNXB;TNXB         | -0.0001531| 0.8101197|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 180049|chr6 | 32014663|TNXB;TNXB         |  0.0029876| 0.0015277|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 172312|chr6 | 32014674|TNXB;TNXB         |  0.0009082| 0.3807045|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 144588|chr6 | 32014893|TNXB;TNXB         |  0.0021443| 0.0006677|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 254656|chr6 | 32014926|TNXB;TNXB         |  0.0007943| 0.0784213|  32013699| 32017224| -5.943456|    0.0010452|
|      4|  96182|chr6 | 32015083|TNXB;TNXB         |  0.0012410| 0.1573622|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 251549|chr6 | 32015175|TNXB;TNXB         |  0.0004095| 0.4306034|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 228450|chr6 | 32015215|TNXB;TNXB         |  0.0004607| 0.3710192|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 214234|chr6 | 32015297|TNXB;TNXB         |  0.0030980| 0.0001455|  32013699| 32017224| -5.943456|    0.0010452|
|      4|  53331|chr6 | 32015618|TNXB              |  0.0001046| 0.8517118|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 190142|chr6 | 32015688|TNXB              |  0.0010612| 0.0066559|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 141947|chr6 | 32015713|TNXB              |  0.0014876| 0.0432615|  32013699| 32017224| -5.943456|    0.0010452|
|      4|  98162|chr6 | 32015737|TNXB              |  0.0013062| 0.0218474|  32013699| 32017224| -5.943456|    0.0010452|
|      4|  46174|chr6 | 32015773|TNXB              |  0.0021020| 0.0396521|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 253689|chr6 | 32016070|TNXB              |  0.0004271| 0.4273317|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 289210|chr6 | 32016100|TNXB              |  0.0001844| 0.5002354|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 288177|chr6 | 32016115|TNXB              | -0.0002369| 0.4819618|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 312500|chr6 | 32016172|TNXB              |  0.0004049| 0.5155801|  32013699| 32017224| -5.943456|    0.0010452|
|      4|  49395|chr6 | 32016214|TNXB              |  0.0016431| 0.0583541|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 215071|chr6 | 32016239|TNXB              |  0.0025661| 0.1157791|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 166791|chr6 | 32016247|TNXB              |  0.0019620| 0.0495010|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 108179|chr6 | 32016257|TNXB              |  0.0024417| 0.0048269|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 291014|chr6 | 32016288|TNXB              |  0.0033794| 0.0000426|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 243825|chr6 | 32016290|TNXB              |  0.0019575| 0.0083271|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 104944|chr6 | 32016360|TNXB              |  0.0000404| 0.9559929|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 180966|chr6 | 32016368|TNXB              |  0.0012041| 0.2480942|  32013699| 32017224| -5.943456|    0.0010452|
|      4|  85924|chr6 | 32016426|TNXB              |  0.0021802| 0.0132794|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 198624|chr6 | 32016473|TNXB              |  0.0009756| 0.0614283|  32013699| 32017224| -5.943456|    0.0010452|
|      4|  15233|chr6 | 32016520|TNXB              |  0.0002999| 0.5368667|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 283755|chr6 | 32016535|TNXB              | -0.0005487| 0.4551498|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 155224|chr6 | 32016678|TNXB              |  0.0003520| 0.5391198|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 216005|chr6 | 32016690|TNXB              | -0.0002819| 0.6021935|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 237173|chr6 | 32016803|TNXB              |  0.0015205| 0.0448445|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 296767|chr6 | 32017024|TNXB              |  0.0015014| 0.1116998|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 204967|chr6 | 32017079|TNXB              |  0.0010082| 0.0341457|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 290092|chr6 | 32017089|TNXB              |  0.0009824| 0.0197532|  32013699| 32017224| -5.943456|    0.0010452|
|      4| 233696|chr6 | 32017224|TNXB              |  0.0006594| 0.2466544|  32013699| 32017224| -5.943456|    0.0010452|

