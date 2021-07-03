# CVETool
CVETool(composite variant expression tool) is an automatic workflow for performing eQTL analysis for composite variant, such as the two haplotype configuration of one with an SNV and the other with an SV.

It can accurately genotype the configuration and assumed an additive linear model to evaluate the association between genotype and quantitative traits.

## Installation

CVETool does not need to be installed. You need to replace the software path in the parameter.txt file with your own software path.

## Requirements

Python 3

Perl with List::Util

R with ggplot2, dplyr and forcats

CNVnator

freebayes

## CVETool usage

```
usage:   CVETool.sh [options]
```
Required arguments
```
   -b FILE  sample ID & bam file (tab separated)
   -e FILE  gene expression values (tab separated)
   -c FILE  covariates values (tab separated)
   -r FILE  reference file
   -v chrN:pos1-pos2   SV region
   -p chrN:pos1        SNP position
```
Optional arguments
```
   -r FILE  sample ID & bam file (tab separated)
   -n       taskname
   -t       parallel tasks
```
