Advanced RNA-seq training - Timecourse
========================================================
author: Damir Baranasic
date: 01.12.2017
autosize: true

First Slide
========================================================

For more details on authoring R presentations please visit <https://support.rstudio.com/hc/en-us/articles/200486468>.

- Bullet 1
- Bullet 2
- Bullet 3

Installed packages
========================================================

Installed with `install.packages`:
- tidyverse
- ggplot2

Installed with `biocInstaller`:
- maSigPro
- DESeq2

Slide With Code
========================================================


```r
summary(cars)
```

```
     speed           dist       
 Min.   : 4.0   Min.   :  2.00  
 1st Qu.:12.0   1st Qu.: 26.00  
 Median :15.0   Median : 36.00  
 Mean   :15.4   Mean   : 42.98  
 3rd Qu.:19.0   3rd Qu.: 56.00  
 Max.   :25.0   Max.   :120.00  
```

Slide With Plot
========================================================

![plot of chunk unnamed-chunk-2](Advanced_RNA-seq_training_Timecourse-figure/unnamed-chunk-2-1.png)

RNA-seq
========================================================

- allele specific expression
- gene fusion
- lncRNA
- eRNA
- alterantively spliced variants

Time course experiments
========================================================

1. Single time series
  - one condition
  - all time points compared to the first one (control)
2. Multi time series
  - several conditions simultaneously
  - controls are sampled over time with the samples
3. Periodicity and cyclic time series
  - sigle or multiple conditions
  - reoccuring exoression patterns and their difference between conditions
  - complex --> a lot of samples needed and synchronisation
  
Experimental worklflow
========================================================

Similar to static algorithms

// insert image here

Experimental design
========================================================

Critical:
  - number of time points
  - number of replicates
  
cost ~ statistical power

- there are tools to estimate this parameters, but they don't consider multi-factor experiments

Generaly, when in doubt:
  - more replicates better than greater sequencing depth
  
Bad design:
  - - statistical power
  - + number of FP
  
Analysis
========================================================

Static tools:
  - sequencing depth and library size
  - batch effect --> protocol, sequencing platform, technical variability etc.
  
Do not consider Correlation between neighbouring time points!

Questions for analysis
========================================================

Number of replicates?

Experimantal design? (two-way or multi-factor)

Differential expression of RNA isoforms?
