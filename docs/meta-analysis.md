---
title: "Meta-analysis using dmrff" 
output:
  pdf_document: default
  word_document:
    highlight: tango
  html_document:
    toc: true
---

# Meta-analysis using dmrff

## Download and prepare DNA methylation datasets

We'll use a couple of small cord blood DNA methylation datasets
available on GEO: GSE79056, GSE62924.

*Loading each dataset takes about 1 minute.*

```r
accessions <- c("GSE79056", "GSE62924")
library(geograbi) ## https://github.com/yousefi138
samples <- sapply(accessions, geograbi.get.samples, simplify=F)
vars <- sapply(samples, geograbi.extract.characteristics, simplify=F)
methylation <- sapply(accessions, geograbi.get.data, simplify=F)
```

Remove non-CpG probes from the datasets.

```r
for (acc in accessions)
    methylation[[acc]] <- methylation[[acc]][grepl("^c", rownames(methylation[[acc]])),]
```

Extract the variable of interest from each dataset.

```r
variable <- sapply(accessions, function(acc) {
    idx <- grep("(ga weeks|gestational_age|gestational age)", colnames(vars[[acc]]))
    stopifnot(length(idx) == 1)
    as.numeric(vars[[acc]][,idx])
}, simplify=F)
```

Remove samples with incomplete data.

```r
for (acc in accessions) {
    idx <- which(!is.na(variable[[acc]]))
    samples[[acc]] <- samples[[acc]][idx,]
    vars[[acc]] <- vars[[acc]][idx,]
    variable[[acc]] <- variable[[acc]][idx]
    methylation[[acc]] <- methylation[[acc]][,idx]
}
```

## Calculate surrogate variables to handle unknown confounding

*SVA takes a few seconds for each dataset.*

```r
library(sva, quietly=T)
```

```r
set.seed(20191029)
covariates <- sapply(accessions, function(acc) {
    mod <- model.matrix(~ var, data.frame(var=variable[[acc]]))
    mod0 <- mod[,1,drop=F]
    random.idx <- sample(1:nrow(methylation[[acc]]), 5000)
    meth.sva <- methylation[[acc]][random.idx,]
    meth.mean <- rowMeans(meth.sva, na.rm=T)
    idx <- which(is.na(meth.sva), arr.ind=T)
    if (nrow(idx) > 0)
        meth.sva[idx] <- meth.mean[idx[,1]]
    sva(meth.sva, mod=mod, mod0=mod0)$sv
}, simplify=F)
```

## Test the associations with gestational age

*Testing associations takes a few seconds for each dataset.*

```r
library(limma)
stats <- sapply(accessions, function(acc) {
    mod <- model.matrix(~ ., data.frame(ga=variable[[acc]],
                                        covariates[[acc]]))
    fit <- lmFit(methylation[[acc]], mod)
    fit <- eBayes(fit)
    data.frame(estimate=fit$coefficients[,"ga"],
               se=sqrt(fit$s2.post) * fit$stdev.unscaled[,"ga"],
               p.value=fit$p.value[,"ga"])
},simplify=F)
```

Now we check the agreement between the datasets.

```r
kable(sapply(stats, function(a) sapply(stats, function(b) {
    sites <- intersect(rownames(a), rownames(b))
    sites <- sites[order(a[sites,"p.value"])[1:100]]
    cor(a[sites,"estimate"], b[sites,"estimate"])
})))
```



|         |  GSE79056|  GSE62924|
|:--------|---------:|---------:|
|GSE79056 | 1.0000000| 0.5645889|
|GSE62924 | 0.8821911| 1.0000000|

A previous study reported thousands of associations with gestational age:
> Bohlin J, et al.
> Prediction of gestational age based on genome-wide
> differentially methylated regions.
> Genome Biol. 2016;17(1):207.

The effect estimates are strongly associated.

```r
bohlin <- read.csv("bohlin.csv",stringsAsFactors=F,row.names=1)
kable(sapply(stats, function(stats) {
    cor(bohlin$estimate, stats[rownames(bohlin),"estimate"], use="p")
}))
```



|         |         x|
|:--------|---------:|
|GSE79056 | 0.7668754|
|GSE62924 | 0.5416417|

## Annotate summary statistics with genomic locations

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
stats <- sapply(stats, function(stats) {
    annotation <- annotation[match(rownames(stats), rownames(annotation)),]
    cbind(stats, annotation)
}, simplify=F)
```

## Construct 'pre' objects for DMR meta-analysis

*Takes about 2 minutes per dataset to calculate CpG site correlations.*

```r
library(dmrff)
pre <- sapply(accessions, function(acc) {
    cat(date(), "pre", acc, "\n")
    with(stats[[acc]], dmrff.pre(estimate, se, methylation[[acc]], chr, pos))
}, simplify=F)
```

## Perform DMR meta-analysis

*There are hundreds of DMRs for this phenotype,
so the final meta-analysis takes about 15 minutes
with a single processor.*

```r
options(mc.cores=20) ## if you have 20 processors available
meta <- dmrff.meta(pre)
```

The output contains a data frame of the meta-analysed regions
and the EWAS meta-analysis on which it was based.


```r
kable(meta$dmrs[1:2,])
```



|chr  |    start|      end|  n|B  |S  |  estimate|        se|         z| p.value| p.adjust|
|:----|--------:|--------:|--:|:--|:--|---------:|---------:|---------:|-------:|--------:|
|chr6 | 32120863| 32120878|  2|NA |NA | 0.0192222| 0.0019000| 10.117026|       0|        0|
|chr6 | 32120625| 32120625|  1|NA |NA | 0.0114568| 0.0014242|  8.044579|       0|        0|


```r
kable(meta$ewas[1:2,])
```



|   estimate|        se|         z|   p.value|chr  |   pos|
|----------:|---------:|---------:|---------:|:----|-----:|
|  0.0005157| 0.0008855|  0.582327| 0.5603465|chr1 | 15865|
| -0.0005600| 0.0010924| -0.512686| 0.6081710|chr1 | 18827|

Here we add the CpG sites and gene names to the EWAS sites:

```r
idx <- match(with(meta$ewas, paste(chr,pos)),
             with(annotation, paste(chr,pos)))
meta$ewas$cpg <- rownames(annotation)[idx]
meta$ewas$gene <- annotation$UCSC_RefGene_Name[idx]
```

## Meta-analyzed DMRs


```r
dmrs <- meta$dmrs[which(meta$dmrs$p.adjust < 0.05 & meta$dmrs$n >= 2), ]
```

We've identified 887 DMRs with Bonferroni adjusted p < 0.05.

For convenience, we'll just look at DMRs on chromosome 6 with at least 5 CpG sites.

```r
dmrs6 <- meta$dmrs[which(meta$dmrs$p.adjust < 0.05
                         & meta$dmrs$n >= 5
                         & meta$dmrs$chr == "chr6"), ]
dmrs6 <- dmrs6[order(dmrs6$start),]
kable(dmrs6[,c("chr","start","end","n","estimate","se","p.value","p.adjust")])
```



|     |chr  |    start|      end|  n|   estimate|        se| p.value|  p.adjust|
|:----|:----|--------:|--------:|--:|----------:|---------:|-------:|---------:|
|3485 |chr6 |  1635611|  1635847|  5|  0.0070901| 0.0006276|       0| 0.0000000|
|3508 |chr6 | 29576422| 29577348|  7|  0.0013228| 0.0001990|       0| 0.0000198|
|1765 |chr6 | 30950666| 30951917|  9|  0.0018605| 0.0003175|       0| 0.0030337|
|3563 |chr6 | 31081905| 31082534|  5|  0.0025786| 0.0004701|       0| 0.0271405|
|458  |chr6 | 31127159| 31127863|  8| -0.0018759| 0.0002991|       0| 0.0002347|
|161  |chr6 | 32042947| 32044104|  7|  0.0032088| 0.0003913|       0| 0.0000000|
|170  |chr6 | 32044254| 32044869|  5|  0.0018827| 0.0002917|       0| 0.0000716|
|432  |chr6 | 32048449| 32049373| 12|  0.0016577| 0.0002922|       0| 0.0092113|
|2408 |chr6 | 32135715| 32136052|  6|  0.0043727| 0.0005362|       0| 0.0000000|

Below we use the `dmrff.sites` function to show the CpG sites in each DMR and their
meta-analysed summary statistics. 

```r
sites <- dmrff.sites(dmrs6, meta$ewas$chr, meta$ewas$pos)
sites <- cbind(sites[,c("region","chr","pos")],
               meta$ewas[sites$site,c("cpg", "estimate","se","p.value")])
kable(sites,row.names=F)
```



| region|chr  |      pos|cpg        |   estimate|        se|   p.value|
|------:|:----|--------:|:----------|----------:|---------:|---------:|
|      1|chr6 |  1635611|cg07714812 |  0.0055709| 0.0007271| 0.0000000|
|      1|chr6 |  1635640|cg19295314 |  0.0055100| 0.0008831| 0.0000000|
|      1|chr6 |  1635808|cg13394216 |  0.0073506| 0.0007248| 0.0000000|
|      1|chr6 |  1635818|cg26908825 |  0.0063406| 0.0008468| 0.0000000|
|      1|chr6 |  1635847|cg10707788 |  0.0079302| 0.0007472| 0.0000000|
|      2|chr6 | 29576422|cg00758854 |  0.0014915| 0.0005716| 0.0090670|
|      2|chr6 | 29576818|cg14111380 |  0.0013964| 0.0004779| 0.0034759|
|      2|chr6 | 29576987|cg15154411 |  0.0021600| 0.0003884| 0.0000000|
|      2|chr6 | 29577006|cg04960880 | -0.0002030| 0.0004939| 0.6810776|
|      2|chr6 | 29577082|cg20642417 |  0.0011591| 0.0004367| 0.0079565|
|      2|chr6 | 29577223|cg03258475 |  0.0015119| 0.0007193| 0.0355717|
|      2|chr6 | 29577348|cg05972518 |  0.0016953| 0.0007676| 0.0272032|
|      3|chr6 | 30950666|cg11062798 |  0.0020931| 0.0003904| 0.0000001|
|      3|chr6 | 30950713|cg15167736 |  0.0009189| 0.0004394| 0.0364842|
|      3|chr6 | 30950779|cg13464738 |  0.0016306| 0.0005485| 0.0029488|
|      3|chr6 | 30951084|cg06183469 |  0.0015434| 0.0005918| 0.0091123|
|      3|chr6 | 30951225|cg20759486 |  0.0004405| 0.0007117| 0.5360141|
|      3|chr6 | 30951377|cg15442792 |  0.0014674| 0.0006269| 0.0192352|
|      3|chr6 | 30951394|cg07538160 |  0.0009127| 0.0008493| 0.2825756|
|      3|chr6 | 30951604|cg24311704 |  0.0009916| 0.0004889| 0.0425183|
|      3|chr6 | 30951917|cg04230397 |  0.0012643| 0.0005311| 0.0172935|
|      4|chr6 | 31081905|cg21351647 |  0.0025893| 0.0005321| 0.0000011|
|      4|chr6 | 31082187|cg24926791 | -0.0069040| 0.0077963| 0.3758575|
|      4|chr6 | 31082196|cg11774057 |  0.0025701| 0.0011739| 0.0285675|
|      4|chr6 | 31082200|cg17849733 | -0.0016192| 0.0022597| 0.4736383|
|      4|chr6 | 31082534|cg01016122 |  0.0038909| 0.0008195| 0.0000021|
|      5|chr6 | 31127159|cg23679615 | -0.0034964| 0.0008678| 0.0000560|
|      5|chr6 | 31127173|cg21863888 | -0.0030393| 0.0008101| 0.0001757|
|      5|chr6 | 31127178|cg09045681 | -0.0030642| 0.0009141| 0.0008022|
|      5|chr6 | 31127271|cg10094358 | -0.0031383| 0.0007385| 0.0000214|
|      5|chr6 | 31127357|cg09329266 | -0.0026323| 0.0007889| 0.0008480|
|      5|chr6 | 31127379|cg21799473 | -0.0030396| 0.0008470| 0.0003327|
|      5|chr6 | 31127527|cg23654219 | -0.0016154| 0.0004102| 0.0000820|
|      5|chr6 | 31127863|cg16095155 | -0.0021249| 0.0004718| 0.0000067|
|      6|chr6 | 32042947|cg20170809 |  0.0028046| 0.0005033| 0.0000000|
|      6|chr6 | 32043029|cg14173662 |  0.0022512| 0.0008114| 0.0055315|
|      6|chr6 | 32043191|cg06878311 |  0.0011375| 0.0005462| 0.0372982|
|      6|chr6 | 32043561|cg22825871 |  0.0026882| 0.0004321| 0.0000000|
|      6|chr6 | 32043676|cg25962829 |  0.0014741| 0.0005817| 0.0112724|
|      6|chr6 | 32043739|cg04492496 |  0.0035224| 0.0004740| 0.0000000|
|      6|chr6 | 32044104|cg15318957 |  0.0007906| 0.0005847| 0.1763435|
|      7|chr6 | 32044254|cg05423418 |  0.0021740| 0.0006335| 0.0005995|
|      7|chr6 | 32044384|cg05569328 |  0.0020299| 0.0006655| 0.0022868|
|      7|chr6 | 32044404|cg22851080 |  0.0017799| 0.0004356| 0.0000439|
|      7|chr6 | 32044496|cg27428104 |  0.0023844| 0.0009032| 0.0082886|
|      7|chr6 | 32044869|cg23164535 |  0.0022082| 0.0007485| 0.0031747|
|      8|chr6 | 32048449|cg18431489 |  0.0019926| 0.0006965| 0.0042220|
|      8|chr6 | 32048463|cg15212833 |  0.0011700| 0.0004530| 0.0098059|
|      8|chr6 | 32048632|cg06696874 |  0.0014025| 0.0005259| 0.0076543|
|      8|chr6 | 32048679|cg02550722 |  0.0022442| 0.0005690| 0.0000800|
|      8|chr6 | 32049053|cg27624229 |  0.0019113| 0.0008324| 0.0216721|
|      8|chr6 | 32049177|cg00661399 |  0.0018632| 0.0010408| 0.0734411|
|      8|chr6 | 32049196|cg13199127 |  0.0022885| 0.0008889| 0.0100375|
|      8|chr6 | 32049235|cg17031787 |  0.0024231| 0.0007054| 0.0005924|
|      8|chr6 | 32049263|cg15712520 |  0.0006509| 0.0006902| 0.3456530|
|      8|chr6 | 32049322|cg09975576 |  0.0022240| 0.0004347| 0.0000003|
|      8|chr6 | 32049354|cg07306331 |  0.0006079| 0.0006300| 0.3345328|
|      8|chr6 | 32049373|cg15048806 |  0.0008268| 0.0005814| 0.1549715|
|      9|chr6 | 32135715|cg12305588 |  0.0046133| 0.0006954| 0.0000000|
|      9|chr6 | 32135718|cg08759957 |  0.0045660| 0.0006664| 0.0000000|
|      9|chr6 | 32135728|cg12827283 |  0.0058240| 0.0010218| 0.0000000|
|      9|chr6 | 32135776|cg11356887 |  0.0041585| 0.0009169| 0.0000058|
|      9|chr6 | 32135803|cg26785234 |  0.0049980| 0.0007865| 0.0000000|
|      9|chr6 | 32136052|cg00529757 |  0.0047038| 0.0008191| 0.0000000|

We can use the same function to show the summary statistics for these sites
from one of the original studies.

```r
acc <- "GSE79056"
sites <- dmrff.sites(dmrs6, stats[[acc]]$chr, stats[[acc]]$pos)
sites <- cbind(sites[,c("region","chr","pos")],
               stats[[acc]][sites$site,c("estimate","se","p.value")])
kable(sites,row.names=F)
```



| region|chr  |      pos|   estimate|        se|   p.value|
|------:|:----|--------:|----------:|---------:|---------:|
|      1|chr6 |  1635611|  0.0051696| 0.0007719| 0.0000001|
|      1|chr6 |  1635640|  0.0058897| 0.0009275| 0.0000003|
|      1|chr6 |  1635808|  0.0074477| 0.0007734| 0.0000000|
|      1|chr6 |  1635818|  0.0065229| 0.0008696| 0.0000000|
|      1|chr6 |  1635847|  0.0084053| 0.0007929| 0.0000000|
|      2|chr6 | 29576422|  0.0015431| 0.0005857| 0.0127382|
|      2|chr6 | 29576818|  0.0012675| 0.0004897| 0.0142397|
|      2|chr6 | 29576987|  0.0021279| 0.0003930| 0.0000054|
|      2|chr6 | 29577006| -0.0000880| 0.0005168| 0.8658845|
|      2|chr6 | 29577082|  0.0009461| 0.0004556| 0.0457304|
|      2|chr6 | 29577223|  0.0015131| 0.0007329| 0.0469049|
|      2|chr6 | 29577348|  0.0017563| 0.0007934| 0.0338776|
|      3|chr6 | 30950666|  0.0020634| 0.0003983| 0.0000118|
|      3|chr6 | 30950713|  0.0008106| 0.0004478| 0.0793752|
|      3|chr6 | 30950779|  0.0016826| 0.0005570| 0.0048440|
|      3|chr6 | 30951084|  0.0011679| 0.0006296| 0.0728457|
|      3|chr6 | 30951225|  0.0003523| 0.0007366| 0.6356345|
|      3|chr6 | 30951377|  0.0016004| 0.0006378| 0.0171952|
|      3|chr6 | 30951394|  0.0011478| 0.0008681| 0.1952200|
|      3|chr6 | 30951604|  0.0010759| 0.0005054| 0.0408556|
|      3|chr6 | 30951917|  0.0012252| 0.0005371| 0.0291383|
|      4|chr6 | 31081905|  0.0026277| 0.0005462| 0.0000345|
|      4|chr6 | 31082187| -0.0088934| 0.0080690| 0.2783791|
|      4|chr6 | 31082196|  0.0026329| 0.0012948| 0.0501345|
|      4|chr6 | 31082200| -0.0016841| 0.0023681| 0.4820069|
|      4|chr6 | 31082534|  0.0037538| 0.0008602| 0.0001188|
|      5|chr6 | 31127159| -0.0033843| 0.0008846| 0.0005510|
|      5|chr6 | 31127173| -0.0029682| 0.0008271| 0.0010652|
|      5|chr6 | 31127178| -0.0029469| 0.0009392| 0.0035780|
|      5|chr6 | 31127271| -0.0029435| 0.0007600| 0.0004827|
|      5|chr6 | 31127357| -0.0024225| 0.0008110| 0.0052850|
|      5|chr6 | 31127379| -0.0029386| 0.0008708| 0.0019063|
|      5|chr6 | 31127527| -0.0015538| 0.0004264| 0.0009136|
|      5|chr6 | 31127863| -0.0022210| 0.0004918| 0.0000764|
|      6|chr6 | 32042947|  0.0027819| 0.0005165| 0.0000059|
|      6|chr6 | 32043029|  0.0022362| 0.0008547| 0.0133117|
|      6|chr6 | 32043191|  0.0010619| 0.0005604| 0.0668866|
|      6|chr6 | 32043561|  0.0024802| 0.0004462| 0.0000036|
|      6|chr6 | 32043676|  0.0016935| 0.0006061| 0.0086099|
|      6|chr6 | 32043739|  0.0035517| 0.0004893| 0.0000000|
|      6|chr6 | 32044104|  0.0009878| 0.0006038| 0.1113933|
|      7|chr6 | 32044254|  0.0020807| 0.0006542| 0.0031964|
|      7|chr6 | 32044384|  0.0021768| 0.0006997| 0.0039069|
|      7|chr6 | 32044404|  0.0019054| 0.0004670| 0.0002685|
|      7|chr6 | 32044496|  0.0024133| 0.0009298| 0.0139925|
|      7|chr6 | 32044869|  0.0020486| 0.0007735| 0.0123092|
|      8|chr6 | 32048449|  0.0022535| 0.0007329| 0.0042130|
|      8|chr6 | 32048463|  0.0012539| 0.0004629| 0.0106282|
|      8|chr6 | 32048632|  0.0013326| 0.0005593| 0.0231078|
|      8|chr6 | 32048679|  0.0024140| 0.0005817| 0.0002200|
|      8|chr6 | 32049053|  0.0019729| 0.0008879| 0.0332687|
|      8|chr6 | 32049177|  0.0015578| 0.0010768| 0.1577390|
|      8|chr6 | 32049196|  0.0021966| 0.0009272| 0.0238531|
|      8|chr6 | 32049235|  0.0022975| 0.0007260| 0.0033314|
|      8|chr6 | 32049263|  0.0003256| 0.0007372| 0.6616342|
|      8|chr6 | 32049322|  0.0021926| 0.0004535| 0.0000301|
|      8|chr6 | 32049354|  0.0008003| 0.0006745| 0.2438739|
|      8|chr6 | 32049373|  0.0004626| 0.0006743| 0.4974989|
|      9|chr6 | 32135715|  0.0045975| 0.0006997| 0.0000002|
|      9|chr6 | 32135718|  0.0045876| 0.0006733| 0.0000001|
|      9|chr6 | 32135728|  0.0056907| 0.0010408| 0.0000047|
|      9|chr6 | 32135776|  0.0040638| 0.0009346| 0.0001244|
|      9|chr6 | 32135803|  0.0050073| 0.0008067| 0.0000005|
|      9|chr6 | 32136052|  0.0046870| 0.0008265| 0.0000026|

Going back to full set of DMRs,
we can ask which regions are novel,
i.e. contain no CpG sites identified by Bohlin et al. 

```r
## consider only autosomal regions
dmrs <- dmrs[which(dmrs$chr %in% paste0("chr", 1:22)),]
## collect CpG sites in the regions
sites <- dmrff.sites(dmrs, meta$ewas$chr, meta$ewas$pos)
sites <- cbind(sites, meta$ewas[sites$site,c("cpg","estimate","se","p.value")])
## identify sites with associations according to Bohlin et al.
sites$bohlin <- sites$cpg %in% rownames(bohlin)
## novel regions do not contain any of the Bohlin sites
dmrs$novel <- !(1:nrow(dmrs) %in% sites$region[sites$bohlin])
```

Of 863 DMRs identified, 
28.62%
are novel.

They appear on the following chromosomes:

```r
freq.novel <- table(dmrs$chr[dmrs$novel])
freq <- table(dmrs$chr)
prop <- 100*freq.novel/freq
prop <- as.data.frame(prop)
colnames(prop) <- c("chr","pct")
prop$novel <- freq.novel
prop$total <- freq
prop <- prop[order(prop$pct),]
kable(prop,row.names=F,digits=2)
```



|chr   |   pct| novel| total|
|:-----|-----:|-----:|-----:|
|chr22 | 11.11|     3|    27|
|chr10 | 16.22|     6|    37|
|chr3  | 17.65|     6|    34|
|chr9  | 18.75|     3|    16|
|chr21 | 22.22|     2|     9|
|chr5  | 22.73|     5|    22|
|chr19 | 22.92|    11|    48|
|chr14 | 25.00|     7|    28|
|chr13 | 26.32|     5|    19|
|chr2  | 26.32|    10|    38|
|chr17 | 28.24|    24|    85|
|chr15 | 29.17|     7|    24|
|chr6  | 30.77|    32|   104|
|chr11 | 31.17|    24|    77|
|chr7  | 31.82|    14|    44|
|chr16 | 32.56|    14|    43|
|chr1  | 32.91|    26|    79|
|chr12 | 33.33|    14|    42|
|chr20 | 33.33|     7|    21|
|chr8  | 34.38|    11|    32|
|chr4  | 46.15|    12|    26|
|chr18 | 50.00|     4|     8|
|chrX  |   NaN|     0|     0|
|chrY  |   NaN|     0|     0|



