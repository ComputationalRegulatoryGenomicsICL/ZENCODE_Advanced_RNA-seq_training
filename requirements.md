Requirements
====

All that is needed to download and install prior to the workshop!

<br>

## Data, presentations, tutorial


Can be found on our group github page [here](https://github.com/ComputationalRegulatoryGenomicsICL/ZENCODE_Advanced_RNA-seq_training)

## R

Download or please update R to the latest version 3.4.2) from [here](https://www.r-project.org/).

## RStudio


Download and install RStudio from the following link (if you already have RStudio, update to latest version): [RStudio weblink](https://www.rstudio.com/products/RStudio/#Desktop).

## R Packages


If asked to update packages `Update all/some/none [a/s/n]` - update all packages (answer a).

### Bioconductor


-   Update Bioconductor

``` r
source("https://bioconductor.org/biocLite.R")
biocLite()
```

-   The rest of the packages

``` r
biocLite("DESeq2")
biocLite("DEXSeq")
biocLite("maSigPro")
biocLite("clusterProfiler")
biocLite("org.Dr.eg.db")
```

### CRAN


``` r
install.packages("RColorBrewer")
install.packages("gplots")
install.packages("ggplot2")
install.packages("venn")
install.packages("tidyverse")
```
