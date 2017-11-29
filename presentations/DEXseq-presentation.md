Inferring differential exon usage in RNA-Seq data with the DEXSeq package
========================================================
author: Leonie Roos & Nejc Haberman
date: 1st of December 2017
autosize: true

Alternative splicing
========================================================

- In higher eukaryotes, most genes have several isoforms.
- So far, we counted reads in genes. To study alternative splicing, 
we need read counts per exon across all isoforms. 
-  To assess the significance of differences
to isoform ratios between conditions, the
assignment uncertainty has to be combined with
the noise estimates.
- This opens the possibility to study regulation of
isoform abundance ratios. For example: Is a given exon
spliced out more often under one condition than in
control? 

<div align="left">
<img src="./DEXseq-presentation-figure/alternative_splicing.jpeg" width=600 height=400>
</div>
*Nature Reviews Genetics 3, 285–298 (2002)*

DEXseq Bioconductor package
========================================================
- combination of Python scripts and an R package
- Python script to get counting bins from a GTF file
- Python script to get count table from SAM/BAM files
- R functions to set up model frames and perform
 generalized linear models (GLMs)
  - Provides statistical testing which allows us to see 
whether changes in isoform abundances are just random
variation or may be attributed to changes in tissue
type or experimental condition.
  - Testing on the level of individual exons gives power
to study the mechanisms of alternative isoform regulation.
- R functions to visualize results and compile an
HTML report

<div align="left">
<img src="./DEXseq-presentation-figure/DEXseq-example.png" width=500 height=300>
</div>

<small>*http://bioconductor.org/packages/release/bioc/html/DEXSeq.html*</small>


Data preparation 
========================================================
- Annotation conversion into counting exons from GTF to GFF format.


```r
$ python dexseq_prepare_annotation.py <ANNOTATION.GTF> <ANNOTATION.GFF>
```

- Counting reads within exons.

```r
$ python dexseq_count.py -p yes -s no <ANNOTATION.GFF> <INPUT.BAM> <OUTPUT.tab>

  Options:
    -s, whether the data is from a strand-specific assay 
    -p, whether the data is paired-end
    <ANNOTATION.GFF>, pre-processed annotation file using dexseq_prepare_annotation.py python script
    <INPUT.BAM>, mapped reads
```

*http://bioconductor.org/packages/release/bioc/html/DEXSeq.html*

Alternative tools
========================================================
- JunctionSeq (https://bioconductor.org/packages/release/bioc/html/JunctionSeq.html)
  - it is a Bioconductor package for detection and visualization of differential usage of exons and splice junctions in High-Throughput RNA-Seq datasets. 
  - it is heavily based on DEXSeq.
  - characterizing differential usage of transcript exons and/or splice junctions. 
- MAGIQ (https://majiq.biociphers.org/)
  - quantifies the relative inclusion levels of local splicing variations (LSVs) in a given 
experimental condition (also known as “percent spliced in”, relative abundance PSI), 
and quantifying changes of LSVs inclusion levels between two experimental conditions.


Data used in this case study
========================================================
<small>
# Loss of function of myosin chaperones triggers Hsf1-mediated transcriptional response in skeletal muscle cells
## Christelle Etard, Olivier Armant, Urmas Roostalu, Victor Gourain, Marco Ferg and Uwe Strähle
*Genome Biology* 2015 16:267 <https://doi.org/10.1186/s13059-015-0825-8>

```
RNA-seq_Strahle_Lab_0005AS.<SequencingID>.USERvgourain.R.ReadsPerGene.out.tab
```

|Hpf |      wt     |    unc45b   |
|:--:|:-----------:|:-----------:|
| 24 | DCD001548SQ | DCD001560SQ |
|    | DCD001559SQ | DCD001554SQ |
| <strong>48</strong> | <strong>DCD001546SQ</strong> | <strong>DCD001564SQ</strong>|
|    | <strong>DCD001558SQ</strong>| <strong>DCD001555SQ</strong>| 
| 72 | DCD001547SQ | DCD001565SQ |
|    | DCD001545SQ | DCD001551SQ |

All the libraries were unstranded paired-ended sequenced on Illumina HiSeq 2000 producing 50bp long reads.
</small>


Tutorial
========================================================
type: section
<br>
The tutorial that is linked with this presentation:
<br>

[__Tutorials dir__](https://www.dropbox.com/sh/p4tnruoximdieii/AADAzrUz4FDzYRQd02-K10poa?dl=0)

[DESeq-workshop-Tutorial](https://www.dropbox.com/s/jzst8bbkrmy7sa0/DEXseq_tutorial.html?dl=0)

