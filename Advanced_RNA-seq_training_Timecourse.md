Advanced RNA-seq training - Timecourse
========================================================
author: Damir Baranasic
date: 01.12.2017
autosize: true

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

Data used in this case study
========================================================

## Christelle Etard, Olivier Armant, Urmas Roostalu, Victor Gourain, Marco Ferg and Uwe Strähle

# Loss of function of myosin chaperones triggers Hsf1-mediated transcriptional response in skeletal muscle cells

*Genome Biology* 2015 **16**:267 <https://doi.org/10.1186/s13059-015-0825-8>

```
RNA-seq_Strahle_Lab_0005AS.<SequencingID>.USERvgourain.R.ReadsPerGene.out.tab
```

|Hpf |      wt     |    unc45b   |
|:--:|:-----------:|:-----------:|
| 24 | DCD001548SQ | DCD001560SQ |
|    | DCD001559SQ | DCD001554SQ |
| 48 | DCD001546SQ | DCD001564SQ |
|    | DCD001558SQ | DCD001555SQ |
| 72 | DCD001547SQ | DCD001565SQ |
|    | DCD001545SQ | DCD001551SQ |

All the libraries were unstranded paired-ended sequenced on Illumina HiSeq 2000 producing 50bp long reads

Read the reads per gene counts
========================================================


```r
seq_ids = c("DCD001548SQ", "DCD001559SQ", "DCD001546SQ", "DCD001558SQ", "DCD001547SQ", "DCD001545SQ",
            "DCD001560SQ", "DCD001554SQ", "DCD001564SQ", "DCD001555SQ", "DCD001565SQ", "DCD001551SQ")
files = sapply(seq_ids, function(x) paste("data/RNA-seq_Strahle_Lab_0005AS.", x, ".USERvgourain.R.ReadsPerGene.out.tab", sep = ""), USE.NAMES = FALSE)
head(files, n = 4)
```

```
[1] "data/RNA-seq_Strahle_Lab_0005AS.DCD001548SQ.USERvgourain.R.ReadsPerGene.out.tab"
[2] "data/RNA-seq_Strahle_Lab_0005AS.DCD001559SQ.USERvgourain.R.ReadsPerGene.out.tab"
[3] "data/RNA-seq_Strahle_Lab_0005AS.DCD001546SQ.USERvgourain.R.ReadsPerGene.out.tab"
[4] "data/RNA-seq_Strahle_Lab_0005AS.DCD001558SQ.USERvgourain.R.ReadsPerGene.out.tab"
```

```r
readTagsPerGene <- function(x){
  out <- read.csv(x, skip = 4, header = FALSE, sep = "\t", row.names = 1,
                  col.names = c("","totalCount", "forwardCount", "reverseCount"))
  out
}
allReadsCounts <- lapply(files, readTagsPerGene)
geneNames <- row.names(allReadsCounts[[1]])
totalReadsCount <- sapply(allReadsCounts, function(x) x$totalCount)
rownames(totalReadsCount) <- geneNames
head(totalReadsCount, n = 2)
```

```
                   [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
ENSDARG00000104632  100  100   67   30   48   27   75   69   44    26
ENSDARG00000100660  153  153  110   81   74   98  107  131   72    42
                   [,11] [,12]
ENSDARG00000104632    44    54
ENSDARG00000100660   118   167
```

Make design matrix
========================================================


```r
type <- rep(c("wt", "mut"), each = 6)
type
```

```
 [1] "wt"  "wt"  "wt"  "wt"  "wt"  "wt"  "mut" "mut" "mut" "mut" "mut"
[12] "mut"
```

```r
repl <- rep(c("1", "2"), times = 6)
repl
```

```
 [1] "1" "2" "1" "2" "1" "2" "1" "2" "1" "2" "1" "2"
```

```r
time = rep(c("24", "48", "72"), each = 2, times = 2)
time
```

```
 [1] "24" "24" "48" "48" "72" "72" "24" "24" "48" "48" "72" "72"
```

```r
sample <- Map(function(x, y, z) paste(x, y, z, sep = "_"), type, time, repl)
head(sample)
```

```
$wt
[1] "wt_24_1"

$wt
[1] "wt_24_2"

$wt
[1] "wt_48_1"

$wt
[1] "wt_48_2"

$wt
[1] "wt_72_1"

$wt
[1] "wt_72_2"
```

```r
#files <- setNames(files, sample)
coldata <- data.frame(condition = type, time = time)
row.names(coldata) <- sample
head(coldata)
```

```
        condition time
wt_24_1        wt   24
wt_24_2        wt   24
wt_48_1        wt   48
wt_48_2        wt   48
wt_72_1        wt   72
wt_72_2        wt   72
```
