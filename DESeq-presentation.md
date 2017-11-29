Differential expression analysis for sequencing count data with the DESeq2 package
========================================================
author: Leonie Roos & Nejc Haberman
date: 1st of December 2017
autosize: true

RNA-seq analysis
========================================================

- Discovery
  - detect novel transcripts, including messenger RNAs, and long noncoding RNAs, along with other untranslated regions
  - find transcript boundaries
  - find splice junctions
  - comparison
- Given samples from different experimental conditions, find effects
of the treatment on
  - gene expression across the transcriptome
  - isoform abundance ratios, splice patterns, transcript boundaries

RNA-seq analysis
========================================================

<div align="left">
<img src="./DESeq-presentation-figure/RNA-seq.jpeg" width=800 height=800>
</div>
<small>*Nat Rev Genet 10(1):57-63 (2009)*</small>

Normalisation for library size
========================================================

- Assumption: most genes are not differentially expressed between samples. 
- DESeq2 calculates for every sample the 'effective library size' by a scale factor.

```
    Original library size * scale factor = effective library size 
```

DESeq2 will multiply original counts by the sample scaling factor.

Replicates
========================================================

Another essential factor in designing an RNA-seq experiment is the number of replicates. 

- Two replicates permit to
  - globally estimate variation
- Sufficiently many replicates permit to
  - estimate variation for each gene
  - randomize out unknown covariates
  - spot outliers
  - improve precision of expression and fold-change
estimates

Increasing the number of replicates minimizes the false positives and usually leads to more robust outcomes, ensuring meaningful biological interpretation of the results 


DESseq2 Bioconductor package overview
========================================================

- Normalisation
- Differential expression analysis using negative binomial distribution
- Exploring results
- Data transformation and visualisation

<div align="left">
<img src="./DESeq-presentation-figure/MA-plot.png" width=300 height=300>
<img src="./DESeq-presentation-figure/Dispersion-plot.png" width=300 height=300>
<img src="./DESeq-presentation-figure/Heatmap.png" width=400 height=300>
<img src="./DESeq-presentation-figure/PCA-plot.png" width=300 height=300>
</div>

<small>*http://bioconductor.org/packages/release/bioc/html/DESeq2.html*</small>


Data preparation
========================================================
Count tables of reads per gene:
- using STAR alignment with the *'--quantMode GeneCounts'* option (https://github.com/alexdobin/STAR) 

example of *ReadsPerGene.out.tab* table:

```r
$ head ReadsPerGene.out.tab
N_unmapped      1100914 1100914 1100914
N_multimapping  11275176        11275176        11275176
N_noFeature     6090090 14560190        14670407
N_ambiguous     478755  112331  109422
ENSDARG00000104632      26      14      13
ENSDARG00000100660      42      19      27
ENSDARG00000098417      10      8       3
ENSDARG00000100422      2406    1201    1219
ENSDARG00000102128      250     133     118
```

select columns (gene IDs, unstranded reads) and ignore the first 4 lines

```r
$ tail +5 ReadsPerGene.out.tab | awk '{print $1 "\t" $2}' > ReadsPerGene.DESeq.tab
```


Data preparation
========================================================
Count tables of reads per gene:
- using htseq tool (http://htseq.readthedocs.io/en/master/count.html)

```r
$ htseq-count [options] <alignment_files> <gff_file>

  Options:
    -s, whether the data is from a strand-specific assay 
    -r, the alignment have to be sorted either by read name or by alignment position
    -f, input format (SAM/BAM)
```


Alternative tools
========================================================
- Edge R (http://bioconductor.org/packages/release/bioc/html/edgeR.html)
  - based on the negative binomial distribution (same as DESeq)
  - trimmed mean normalisation is used for size factors normalisation
  - filters low count genes and removes genes with uninformative log fold changes
  
- Cuffdiff (http://cole-trapnell-lab.github.io/cufflinks/manual/) 
  - Beta negative binomial distribution
  - normalisation for gene length and library size 
  - reads/fragments per kilobase per milion mapped reads (RPKM/FPKM)


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

[DESeq-workshop-Tutorial](https://www.dropbox.com/s/iwb9sera4fft748/DEseq_tutorial.html?dl=0)

