# QTL mapping of prepulse inhibition phenotypes in F1 panel derived from common inbred mouse strains

This repository contains code and data to accompany publication of the
manuscript (in review), "Genome-wide association for prepulse
inhibition in a panel of inbred F1 mice."

###License

Copyright (c) 2015, Laura Sittig, Kyle Engel and Peter Carbonetto

The ppi-f1-gwas project repository is free software: you can redistribute
it and/or modify it under the terms of the
[GNU General Public License](http://www.gnu.org/licenses/gpl.html) as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
**without any warranty**; without even the implied warranty of
**merchantability** or **fitness for a particular purpose**. See
[LICENSE](LICENSE) for more details.

###Overview of data files

Here is a brief summary of the files in the [data](data) directory:

+ [pheno.csv](data/pheno.csv) Physiological and behavioural phenotype
data collected on F1 mice that are crosses of inbred lab strains.

+ [mda.F1.csv.gz](data/mda.F1.csv.gz) Genotype data at 302,625 SNPs
for F1 crosses of 29 inbred lab strains obtained from Mouse Diversity
Genotyping Array web resource at Jackson Labs.

+ [unc.F1.csv.gz](data/unc.F1.csv.gz) Genotype data at 2,543,560 SNPs
for F1 crosses of 29 inbred lab strains. These data were obtained from
the UNC resource ([link](http://csbio.unc.edu/imputation)).

###Overview of R source code files

Here is a brief summary of the key files in the [code](code) directory:

+ [map.qtls.gemma.R](data/map.qtls.gemma.R) R script for mapping QTLs
in wild-type F1 mice using a linear mixed model (GEMMA) that corrects
for possible confounding to due cryptic relatedness.

+ [finemap.qtls.gemma.R](data/finemap.qtls.gemma.R) R script for
mapping QTLs using the higher resolution UNC SNP data.

+ [functions.gemma.R](data/functions.gemma.R) Defines functions used
to map QTLs separately on each chromosome using GEMMA.

###Credits

The R code implementing the analysis procedures was developed by:<br>
Laura Sittig, Kyle Engel and Peter Carbonetto<br>
Department of Human Genetics<br>
University of Chicago<br>
June 2015

