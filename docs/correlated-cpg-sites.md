## dmrff and CpG site correlation

Here is common question about `dmrff`:
> My EWAS identifies several regions filled with
> CpG sites all weakly associated with my variable of interest.
> Why does dmrff identify some but not all of them
> as differentially methylated regions?

The short answer is that the regions that fail the dmrff
test will tend to have CpG sites that are more dependent on one
another, particularly in how they relate to the variable interest.
In other words, the association of a CpG site with the variable
will tend to disappear when one more more other CpG sites are included
in the model.
As a result, the region as a whole will not explain much more variation than is
explained by any individual CpG site.
The regions that pass the dmrff test will tend to
be composed of CpG sites that each explain a somewhat different
proportion of the variability of the variable of interest.
As a result, the region as a whole will explain more variation
than any individual CpG site.

The simulations below illustrates this behavior.
 
### Generate a dataset with a differentially methylated region of interdependent CpG sites
We'll be generating some random data,
so set the random seed so the output
is reproducible.

```r
set.seed(20200109)
```

Start with a completely random DNA methylation dataset.

```r
meth <- t(sapply(1:1000,function(i) runif(100)))
rownames(meth) <- paste("cg", 1:nrow(meth), sep="")
colnames(meth) <- paste("s", 1:ncol(meth), sep="")
sites <- data.frame(chr=chr <- rep("chr1", nrow(meth)),
                    pos=seq(1,nrow(meth)*50, 50),
                    stringsAsFactors=F)
rownames(sites) <- rownames(meth)
```

Insert a differentially methylated region
into the dataset composed of highly correlated CpG sites 
that are weakly correlated with the variable of interest.

The region is composed of 5 CpG sites.

```r
region <- paste("cg", 501:505, sep="")
region.chr <- sites[region[1],"chr"]
region.start <- min(sites[region,"pos"])
region.end <- max(sites[region,"pos"])
```

The following function can be used to generate a variable
of interest and methylation levels for a requested number of CpG sites
with a given correlation structure.

```r
library(MASS)
generate.random.data <- function(var.r, sites.r, n.samples, n.sites) {
    sigma <- sites.r + matrix(rnorm((n.sites+1)^2, sd=0.02), ncol=n.sites+1)
    sigma[1,] <- var.r + rnorm(ncol(sigma), sd=0.02)
    sigma[,1] <- var.r + rnorm(nrow(sigma), sd=0.02)
    diag(sigma) <- 1
    random.vars <- mvrnorm(n=n.samples, mu=rep(0,ncol(sigma)), Sigma=sigma, empirical=TRUE)
    list(var=random.vars[,1],
         meth=t(random.vars[,-1]))
}
```

CpG sites will have correlation of about R=0.95 between them
and a correlation of R=0.3 with the variable of interest.

```r
random.data <- generate.random.data(0.3, 0.95, ncol(meth), length(region))
var <- random.data$var
meth[region,] <- random.data$meth
```

The correlations between sites are indeed 0.95.

```r
cor(t(meth[region,]))
```

|      | cg501| cg502| cg503| cg504| cg505|
|:-----|-----:|-----:|-----:|-----:|-----:|
|cg501 | 1.000| 0.942| 0.979| 0.944| 0.920|
|cg502 | 0.942| 1.000| 0.951| 0.955| 0.944|
|cg503 | 0.979| 0.951| 1.000| 0.953| 0.935|
|cg504 | 0.944| 0.955| 0.953| 1.000| 0.934|
|cg505 | 0.920| 0.944| 0.935| 0.934| 1.000|

The correlations between the sites and the variable of interest are 0.3.

```r
t(sapply(region, function(site) cor(meth[site,],var)))
```

| cg501| cg502| cg503| cg504| cg505|
|-----:|-----:|-----:|-----:|-----:|
| 0.307| 0.311| 0.286| 0.307| 0.318|

### Perform an EWAS

```r
stats <- t(sapply(rownames(meth), function(site) coef(summary(lm(meth[site,] ~ var)))["var",]))
stats <- as.data.frame(stats)
colnames(stats) <- c("estimate","se","t","p.value")
```

As expected,
the CpG sites are all similarly associated with the variable of interest.

```r
stats[region,]
```

|      |  estimate|        se|        t|   p.value|
|:-----|---------:|---------:|--------:|---------:|
|cg501 | 0.3066443| 0.0961488| 3.189269| 0.0019157|
|cg502 | 0.3113416| 0.0959946| 3.243324| 0.0016159|
|cg503 | 0.2862668| 0.0967878| 2.957676| 0.0038849|
|cg504 | 0.3070926| 0.0961342| 3.194417| 0.0018850|
|cg505 | 0.3180351| 0.0957704| 3.320806| 0.0012618|

### Apply dmrff

```r
library(dmrff, quietly=T)
dmrs <- dmrff(stats$estimate, stats$se, stats$p.value, meth, sites$chr, sites$pos)
```

```
## [dmrff.candidates] Tue Dec 12 18:42:51 2023 Found  50  candidate regions.
```

The following function checks if two regions overlap.
It will be used to identify any regions tested by dmrff
that overlap with the simulated differentially methylated region.

```r
overlaps  <- function(chr1,s1,e1,chr2,s2,e2) {
    chr1 == chr2 & (s1 >= s2 & s1 <= e2
        | e1 >= s2 & e1 <= e2
        | s2 >= s1 & s2 <= e1
        | e2 >= s1 & e2 <= e1)
}
```

Notice how the p-value for any region considered is quite similar to
EWAS p-values for the individual CpG sites in the region.
This is because the CpG sites all explain a similar portion
of the variability of the variable of interest.

```r
overlaps.region <- sapply(1:nrow(dmrs), function(i) {
    overlaps(dmrs$chr[i], dmrs$start[i], dmrs$end[i],
             region.chr, region.start, region.end)
})
dmrs[overlaps.region,c("chr","start","end","n","estimate","se","p.value","p.adjust")]
```

```
##    chr start   end n  estimate         se     p.value p.adjust
## 2 chr1 25001 25201 5 0.3101704 0.09437005 0.001013522        1
```

### Modify the region to make the CpG sites more independent
Replace the methylation levels of the CpG sites in the region
so that each CpG site is weakly associated with the
variable of interest but mostly independent of one another.

```r
random.data <- generate.random.data(0.3, 0.1, ncol(meth), length(region))
var <- random.data$var
meth[region,] <- random.data$meth
```

As before, the correlations between the sites and the variable of interest are roughly 0.3.

```r
t(sapply(region, function(site) cor(meth[site,],var)))
```

| cg501| cg502| cg503| cg504| cg505|
|-----:|-----:|-----:|-----:|-----:|
| 0.274| 0.333| 0.312| 0.309|  0.26|

However, the correlations between the CpG sites are quite low.

```r
cor(t(meth[region,]))
```

|      | cg501| cg502| cg503| cg504| cg505|
|:-----|-----:|-----:|-----:|-----:|-----:|
|cg501 | 1.000| 0.090| 0.109| 0.102| 0.121|
|cg502 | 0.090| 1.000| 0.094| 0.064| 0.125|
|cg503 | 0.109| 0.094| 1.000| 0.106| 0.125|
|cg504 | 0.102| 0.064| 0.106| 1.000| 0.113|
|cg505 | 0.121| 0.125| 0.125| 0.113| 1.000|

### Repeat the EWAS and dmrff analyses

Repeat the EWAS and DMR analyses with the new dataset.

```r
stats <- t(sapply(rownames(meth), function(site) coef(summary(lm(meth[site,] ~ var)))["var",]))
stats <- as.data.frame(stats)
colnames(stats) <- c("estimate","se","t","p.value")
dmrs <- dmrff(stats$estimate, stats$se, stats$p.value, meth, sites$chr, sites$pos)
```

```
## [dmrff.candidates] Tue Dec 12 18:42:52 2023 Found  42  candidate regions.
```

As before,
the CpG sites are all similarly associated with the variable of interest.

```r
stats[region,]
```

|      |  estimate|        se|        t|   p.value|
|:-----|---------:|---------:|--------:|---------:|
|cg501 | 0.2742205| 0.0971430| 2.822854| 0.0057649|
|cg502 | 0.3333426| 0.0952378| 3.500110| 0.0007015|
|cg503 | 0.3117340| 0.0959816| 3.247851| 0.0015929|
|cg504 | 0.3092785| 0.0960626| 3.219551| 0.0017419|
|cg505 | 0.2595833| 0.0975525| 2.660960| 0.0091055|

The p-value for the region is much lower than previously.
This is because the CpG sites are independently associated with the
variable of interest so, combined, they explain much more variation than
any single CpG site.

```r
overlaps.region <- sapply(1:nrow(dmrs), function(i) {
    overlaps(dmrs$chr[i], dmrs$start[i], dmrs$end[i],
             region.chr, region.start, region.end)
})
dmrs[overlaps.region,c("chr","start","end","n","estimate","se","p.value","p.adjust")]
```

```
##    chr start   end n  estimate         se      p.value     p.adjust
## 1 chr1 25001 25201 5 0.2991924 0.05221058 1.001385e-08 1.212677e-05
```
