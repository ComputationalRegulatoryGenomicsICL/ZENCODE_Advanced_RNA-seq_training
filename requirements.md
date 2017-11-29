Data
====

### First part of Friday:

Can be found on [github](https://github.com/ComputationalRegulatoryGenomicsICL/ZENCODE_Advanced_RNA-seq_training). Please click on the `clone/download` button to download the directory or if you're git savy clone it!

### Second part of Friday:

Can be found [here](https://www.dropbox.com/sh/xk5dc6poh8oosoo/AACd-UOXXpY8d7mrXHGXDwtia?dl=0). Please download the whole directory. Scripts and presentations are also on [github](https://github.com/ComputationalRegulatoryGenomicsICL/ZENCODE_Advanced_RNA-seq_training)

R
=

Download or please update R to the latest version 3.4.2) from [here](https://www.r-project.org/).

RStudio
=======

Download and install RStudio from the following link (if you already have RStudio, update to latest version): [RStudio weblink](https://www.rstudio.com/products/RStudio/#Desktop).

R Packages
==========

If asked to update packages `Update all/some/none [a/s/n]` - update all packages (answer a).

Bioconductor
------------

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

CRAN
----

``` r
install.packages("RColorBrewer")
install.packages("gplots")
install.packages("ggplot2")
install.packages("venn")
install.packages("tidyverse")
```
