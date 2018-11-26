# dmrff

Identifying differentially methylated regions (DMRs) efficiently with power and control

*Background.* An epigenome-wide association study (EWAS) tests associations between epigenetic marks
such as DNA methylation and a given phenotype or exposure.
A popular strategy for increasing power is to test associations of epigenetic measurements taken
across genomic regions rather at individual loci.
This strategy has seen some success because patterns of epigenetic marks at neighboring loci
tend to be under similar regulatory control.
Unfortunately, the most commonly used implementations either fail to control
false positive rates (e.g. comb-p) or suffer from low power (e.g. bumphunter).

*dmrff* is a DMR-finding tool that:
- controls false positive rates
- is more powerful than EWAS
- is fast
- can be applied to the summary statistics of any EWAS
- can be used in the context of meta-analysis

A manuscript describing *dmrff* and evaluating its performance will
appear soon on bioRxiv.
