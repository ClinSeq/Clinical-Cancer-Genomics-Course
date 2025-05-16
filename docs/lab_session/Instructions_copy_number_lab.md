---
editor_options:
  chunk_output_type: editor
execute:
  echo: true
  eval: false
title: Copy number analysis
toc-title: Table of contents
---

## Introduction

There are plenty of methods and code packages available for copy number
analysis, and you will likely be able to find one for free that suits
your particular application. Although tailored for different types of
sequence data (or microarrays), general processing steps of all copy
number analysis tools include:

1.  Quantification of sequence read depth (or array intensity)
    throughout the reference genome as a measure of DNA abundance in the
    sample
2.  Removal of systematic noise using features such as GC content and
    mappability
3.  Removal of systematic noise using normal (for example non-cancerous)
    reference samples
4.  Segmentation - partitioning of the reference genome so that each
    segment can be assigned one copy number
5.  Copy number calling - assigning each segment some estimate of the
    number of copies per cell
6.  Combine the measure of DNA abundance with single nucleotide polymorphism         
    (SNP) allele ratios to support copy number estimates

In this exercise you will perform the above steps using R, with a BAM
file of aligned reads and a VCF file of SNPs as input.

## Download data and prepare R

The data used in this exercise are available for download
[here](https://course-cg-5534.s3.amazonaws.com/cnvs/copy_number_files.tar.bz2).
Extract the content on your computer. Before starting, make sure you
have a recent version of R, RStudio, and the following packages
installed:

### From CRAN

-   data.table
-   ggplot2
-   patchwork
-   PSCBS
-   stringr

To install these packages, open RStudio and type in the R console for each one:

::: cell
``` {.r .cell-code}
install.packages("package_name")
```
:::

### From Bioconductor

-   aroma.light 
-   bamsignals
-   biomaRt
-   Biostrings
-   BSgenome
-   BSgenome.Hsapiens.UCSC.hg19
-   DNAcopy
-   GenomicRanges
-   Rsamtools
-   VariantAnnotation

To install these packages, type this once in the R console:

::: cell
``` {.r .cell-code}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```
:::

Then, for each package (note that the installation of some of the packages may
take a few minutes):

::: cell
``` {.r .cell-code}
BiocManager::install("package_name")
```
:::

> Several packages from the
> [Bioconductor](https://www.bioconductor.org/) repository are used in
> this exercise, for example
> [GenomicRanges](http://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesHOWTOs.pdf)
> which facilitates working with anything located on a genome.

## Parse sequence target definition

We use targeted sequencing of a prostate cancer sample in this exercise,
as the resulting data files are of a more convenient size than whole
genome or exome. This hybrid-capture panel covers a few hundred cancer
genes and a few thousands of additional targets intended to improve copy
number and structural variant detection.

In RStudio, navigate to the folder containing the exercise files:

::: cell
``` {.r .cell-code}
setwd('your/path/to/exercise_folder')
```
:::

> Create a new Quarto Document in the exercise folder. Keep your code as
> shown in the template, and write your comments and answers to
> questions between code chunks. Make sure you run the code, investigate
> the results and answer the questions. You can click **Render** to run
> all code and generate a html report. This way, when you are through
> this exercise, your report is done too.

The [BED](https://en.wikipedia.org/wiki/BED_(file_format)) (Browser
Extensible Data) format is a text file format commonly used to store
genomic regions as coordinates and associated annotations. The data are
presented in the form of columns separated by spaces or tabs. This
format was developed during the Human Genome Project and then adopted by
other sequencing projects.

We begin by parsing and investigating the target definition BED file.

::: cell
``` {.r .cell-code}
library(data.table)

# Input file
bed_file <- 'targets.bed'

# Check that it exists
file.exists(bed_file)

# Read it 
targets <- fread(bed_file)

# Check content
targets

# Let the first column be the target number
targets <- cbind(1:nrow(targets),targets)

# Set useful column names
colnames(targets) <- c('target','chromosome','start','end')

# Calculate target length
targets[,length:=end-start]
summary(targets)

# Check target distribution over chromosomes
targets[,table(chromosome)]
```
:::

The R code in this exercise features `data.table` syntax, for example
using `:=` to modify the content of a table. If you are more familiar
with base R or tidyverse syntax, a comparison between them can be found
[here](https://mgimond.github.io/rug_2019_12/Index.html).

> Feel free to complete the exercise using the coding style you prefer.

## Annotate with gene symbols

Before parsing the sequence data, let's assign gene symbols to targets
that overlap a gene. We'll take advantage of databases accessible
through Bioconductor packages. First, we need to create a GenomicRanges
object representing our targets:

::: cell
``` {.r .cell-code}
library(GenomicRanges)

# Make a GRanges object
target_ranges <- makeGRangesFromDataFrame(targets)
target_ranges

# Check the sequence naming style
seqlevelsStyle(target_ranges)
```
:::

As you may recognize, our targets have the NCBI/EMBL sequence naming
style (*1* and not *chr1*, etc). This matches the reference genome used
with our sequence data.

### Download gene coordinates

Let's use the `biomaRt` package to download Ensembl gene coordinates for
reference genome b37/Hg19.

> To learn more about a function, just type ?function_name to open its
> documentation in RStudio's bottom right window.

::: cell
``` {.r .cell-code}
library(biomaRt)

?useEnsembl

# Open a connection
biomart_connection <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)

?getBM

# Get gene coordinates
genes <- getBM(
    attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),
    mart = biomart_connection)

# Convert to data.table
genes <- as.data.table(genes)

genes
```
:::

We can now make a GenomicRanges object for the genes:

::: cell
``` {.r .cell-code}
colnames(genes)[3:5] <- c('chromosome','start','end')

gene_ranges <- makeGRangesFromDataFrame(genes)

seqlevelsStyle(gene_ranges)
```
:::

### Overlap targets and genes

The chromosome naming style of the genes is also NCBI/EMBL. We can now
go ahead and match them up to targets using the findOverlaps function:

::: cell
``` {.r .cell-code}
# Overlap the two GRanges sets
overlap <- findOverlaps(target_ranges,gene_ranges)

overlap
```
:::

Overlaps are presented as pairs of query (first argument) and subject
(second argument) rows for which there is any overlap. We can simply
assign targets their corresponding gene name for each overlap and check
the result:

::: cell
``` {.r .cell-code}
targets$gene <- '' # default if no overlap

# Put a gene symbol on each target
targets[queryHits(overlap),gene:=genes[subjectHits(overlap)]$hgnc_symbol]

# Count targets per gene
as.data.table(table(targets$gene))

as.data.table(sort(table(targets$gene),decreasing = T))[N>25]
```
:::

> There are many ways to achieve the same thing. You should focus mostly
> on *what* we do, and not worry too much about exactly how. Google,
> ChatGTP, etc, are there for you as long as you know what you want to
> achieve.

## Parse sequence read counts

We are now ready to access the BAM file and count the number of reads
(actually read pairs, corresponding to DNA fragments) that map to each
target. For this we use the `bamsignals` package and the `GenomicRanges`
object representing the targets. To count a read pair, we require its
midpoint to fall within a target, and that it is not flagged as a
[duplicate](http://samtools.github.io/hts-specs/SAMv1.pdf).

::: cell
``` {.r .cell-code}
library(bamsignals)
library(Rsamtools)

# Imput bam file
bamfile <- 'sample1.bam'
file.exists(bamfile)
file.exists(paste0(bamfile,'.bai'))

# Count reads in targets
targets$count <- bamCount(bamfile, target_ranges, paired.end="midpoint", filteredFlag=1024)
targets

summary(targets)
```
:::

## Plot sequence read count

Let's investigate the raw sequence read count across targets. You are
probably already familiar with
[ggplot2](https://ggplot2.tidyverse.org/).

::: cell
``` {.r .cell-code}
library(ggplot2)
theme_set(theme_bw())

# Simple ggplot
ggplot(data = targets) + geom_point(mapping = aes(x = target, y = count))

# Some adjustments
ggplot(data = targets) + ylim(c(0,3500)) +
    geom_point(mapping = aes(x = target, y = count),size=.5)
```
:::

Targets are just shown in order. To plot genomic positions without
mixing up chromosomes, we can separate the plot by chromosome:

::: cell
``` {.r .cell-code}
# Plot start pos by chromosome
ggplot(data = targets) + ylim(c(0,3500)) +
    geom_point(mapping = aes(x = start, y = count),size=.5) + 
    facet_wrap(facets = vars(chromosome),ncol = 3)
```
:::

The figure may now not look a bit squished in RStudio's bottom-right
plots tab. Click the **Zoom** button and enlarge the new window.

Let's choose a few genes to highlight in the plot.

::: cell
``` {.r .cell-code}
# Genes of interest
my_genes <- c('AR','ATM','BRCA1','BRCA2','PTEN','TMPRSS2','ERG')
targets[gene %in% my_genes,label:=gene]

ggplot() + ylim(c(0,3500)) +
    geom_point(data = targets, mapping = aes(x = target, y = count, col = label), size=.5)

ggplot(data = targets) + ylim(c(0,3500)) +
    geom_point(mapping = aes(x = start, y = count, col = label), size=.5) +
    facet_wrap(facets = vars(chromosome),ncol = 3)
```
:::

## Model and correct for GC content bias

Sequence GC content is known to affect PCR amplification and sequence
coverage. Let's retrieve the GC content for the targets and investigate
if there is a bias in our data set. `BSgenome` is a Bioconductor package
that provides a framework for working with genomic sequences in R. The
package provides a set of objects, functions, and methods for working
with different reference genomes and their annotations, such as gene
locations, exons, and transcription start sites.

First we check whether a matching reference genome is available:

::: cell
``` {.r .cell-code}
library(BSgenome)

# Check which genomes are available
available.genomes()
```
:::

The best match is `BSgenome.Hsapiens.UCSC.hg19`. But that means our
target `GenomicRanges` object needs to have UCSC-style sequence naming:

::: cell
``` {.r .cell-code}
ucsc_ranges <- target_ranges

# Swap sequence naming style
seqlevelsStyle(ucsc_ranges) <- "UCSC"
ucsc_ranges

seqlevelsStyle(ucsc_ranges)
```
:::

We can now assign GC content to targets:

::: cell
``` {.r .cell-code}
library(BSgenome.Hsapiens.UCSC.hg19)
library(stringr)

# Get the target sequences
dna <- getSeq(Hsapiens, ucsc_ranges)

# Calculate GC content for the target sequences
targets$gc <- str_count(as.character(dna),pattern = '[GC]') / width(dna)

# Plot the GC content and read count
ggplot(data = targets) + ylim(c(0,3500)) +
    geom_point(mapping = aes(x = gc, y = count),alpha=.1)
```
:::

There appears to be a modest effect of target GC content on read count.
Let's build a [loess](https://en.wikipedia.org/wiki/Local_regression)
model of the effect.

::: cell
``` {.r .cell-code}
# Make a loess model
loess_model=loess(count ~ gc, data = targets)
```
:::

To plot the model, we predict the read count from GC:

::: cell
``` {.r .cell-code}
# Add predicted count (model y value) to targets
targets[,gc_predicted_count:=predict(loess_model,gc)]

# Plot with the model prediction
ggplot(data=targets) + ylim(c(0,3500)) +
    geom_point(mapping = aes(x = gc, y = count),alpha=.1,size=1) +
    geom_line(mapping = aes(x = gc,y = gc_predicted_count),col='blue')
```
:::

To remove the GC content bias, we divide the observed count with the
predicted. This new count *ratio* should be a slightly better measure of
DNA abundance in our sample.

::: cell
``` {.r .cell-code}
# Correct for the effect of GC content
targets[,count_ratio:=count/gc_predicted_count]
targets
```
:::

To compare the new metric with the previous, we can plot them
side-by-side. With [patchwork](https://patchwork.data-imaginist.com/), we
can use arithmetic operators to combine and align multiple plots.

::: cell
``` {.r .cell-code}
library(patchwork)

# Plot and align
p1 <- ggplot(data = targets) + 
    geom_point(mapping = aes(x = gc, y = count),alpha=.1)

p2 <- ggplot(data = targets) + 
    geom_point(mapping = aes(x = gc, y = count_ratio),alpha=.1)

p1+p2
```
:::

The GC content correction seems to have had a small but positive effect.
It is common to also use either a matched normal or a pool of normal
reference samples (sequenced similarly) to remove a little bit more
systematic noise. We have omitted that step here.

> **Question 1:** An important reason to perform copy number analysis of
> cancer samples is to find homozygous deletions of tumor suppressor
> genes. Let's assume a patient has a small hemizygous germline deletion
> (loss of one copy) affecting one tumor suppressor gene only. In the
> tumor cell population there is also a somatic deletion of the other,
> homologous, chromosome, including the remaining copy of the gene.
> Would the resulting homozygous deletion still be visible if the read
> count ratio of the tumor sample were to be divided by that of the
> normal sample?

> Write answers to all questions as text in your Quarto document, so
> that they are stored in the report with your data analysis. When done,
> render the report (including all analysis steps you have included) by
> clicking the **Render** button. Submit your rendered report in Canvas
> under the Copy number analysis module.

## Segment the data

Although we can already spot some apparent gains and losses, statistical
tools can help us better define copy number segments, for which we can
then calculate the most likely copy number given the observation.
[Circular binary
segmentation](https://pubmed.ncbi.nlm.nih.gov/15475419/) (**CBS**) is
probably the most commonly used segmentation method. Here we use the
`PSCBS` ("Parent specific" CBS, as it can also use SNP allele data) R
package as a wrapper to perform basic CBS on the GC content-adjusted
count ratio.

CBS requires a DNA abundance measure to be log-transformed, making its
distribution more normal-like. If you have worked with copy number
analysis before, you may recognize the term **log ratio**.

::: cell
``` {.r .cell-code}
library(PSCBS)

# Log transform the read count
targets[,log_ratio:=log2(count_ratio)]
summary(targets)

# Segmentation requires values to be finite
targets[is.infinite(log_ratio),log_ratio:=NA]
summary(targets)

# Segment the log ratio
segments <- segmentByCBS(y=targets$log_ratio)
segments
```
:::

After segmentation we can transform the segment mean values back and
plot the segments with the targets.

::: cell
``` {.r .cell-code}
# Linear space segment means
segments <- as.data.table(segments)[,count_ratio:=2^mean]
segments

ggplot() + ylim(c(0,2.5)) +
    geom_point(data = targets, mapping = aes(x = target, y = count_ratio,fill=chromosome),shape=21) + 
    geom_segment(data=segments,col='green',size=2,
                 mapping = aes(x=start,xend=end,y=count_ratio,yend=count_ratio))
```
:::

Note that the segmentation is based only on the target log ratio and
order, and that segment start and end refer to target number. To avoid
segments spanning across chromosomes, and breakpoints just missing the
chromosome boundary, we can add chromosomes to the segmentation call.
The chromosome vector is required to be numeric.

::: cell
``` {.r .cell-code}
# Segment with chromosomes specified
segments <- segmentByCBS(y=targets$log_ratio,
                         chromosome=as.numeric(str_replace(targets$chromosome,'X','23')))
segments

# Tidy up
segments <- as.data.table(segments)[,count_ratio:=2^mean][!is.na(chromosome),-1]
segments

ggplot() + ylim(c(0,2.5)) +
    geom_point(data = targets, mapping = aes(x = target, y = count_ratio,fill=chromosome),shape=21) + 
    geom_segment(data=segments,col='green',size=2,
                 mapping = aes(x=start,xend=end,y=count_ratio,yend=count_ratio))
```
:::

CBS can also take a vector of chromosomal positions as input, in which
case the resulting segment start and end positions are based on them. We
use the target midpoints:

::: cell
``` {.r .cell-code}
# Segment with position and chromosome positions
segments_pos <- segmentByCBS(y=targets$log_ratio,
                             chromosome=as.numeric(str_replace(targets$chromosome,'X','23')),
                             x=targets$start+60)

segments_pos <- as.data.table(segments_pos)[,count_ratio:=2^mean][!is.na(chromosome),-1]
segments_pos


# convert chromosomes back to NCBI
segments_pos[,chromosome:=str_replace(as.character(chromosome),'23','X')]

ggplot(data = targets) + ylim(c(0,2.5)) +
    geom_point(mapping = aes(x = start, y = count_ratio, fill=label),shape=21) +
    geom_segment(data=segments_pos,col='green',
                 mapping = aes(x=start,xend=end,y=count_ratio,yend=count_ratio)) +
    facet_wrap(facets = vars(chromosome),ncol = 2) +
    theme(panel.spacing = unit(0, "lines"),strip.text.x = element_text(size = 6))
```
:::

Let's now take another look at segmented targets plotted in order, as
they are very unevenly distributed over the genome. We can spot some
apparent copy number losses and gains of genes in our selection. Take a
closer look at BRCA2 and part of PTEN.

::: cell
``` {.r .cell-code}
p1 <- ggplot() + ylim(c(0,2.5)) + 
    geom_point(data = targets, mapping = aes(x = target, y = count_ratio,fill=chromosome),shape=21) + 
    geom_segment(data=segments,col='green',size=2,
                 mapping = aes(x=start,xend=end,y=count_ratio,yend=count_ratio))

p2 <- ggplot() + ylim(c(0,2.5)) + 
    geom_point(data = targets, mapping = aes(x = target, y = count_ratio,fill=label),shape=21) + 
    geom_segment(data=segments,col='green',size=2,
                 mapping = aes(x=start,xend=end,y=count_ratio,yend=count_ratio))

p1/p2
```
:::

> **Question 2:** Target density is only about one per megabase, with
> much higher density in some cancer genes. How many targets do you
> think a deletion would have to cover for us to be able to find it?
> What else might influence sensitivity?

For BRCA2 and (part of) PTEN, count ratio (remember, it is the observed
read count (DNA abundance) relative to the expected, given target GC
content), is near 0.5, or 50%. For a largely diploid genome this is
exactly what to expect if these genomic segments were hemizygously
deleted, i.e., one out of the two homologous copies would have been
lost. This would also fit reasonably well with gain of 1 copy of 8q
(resulting in a count ratio of 1.5) and 3 extra copies of AR (being on
the X chromosome, tumor cells would have quadrupled its single original
AR copy, resulting in the count ratio going from about 0.5 to near 2.).

But this solution would require the tumor cell content ("**purity**") in
the sample to be near 100%, which is rare unless the sample comes from a
cell line. Note also that several other segments have a mean count ratio
that would not fit this solution well. There are actually multiple other
combinations of copy number and purity that would fit the observation
reasonably well.

Our measure of relative DNA abundance - count ratio - is automatically
centered around 1, regardless of the average copy number ("**ploidy**")
of the sequenced genome.

> Let's assume the cancer cells sampled had an average of 4 instead of 2
> copies. This would neither make us sequence to twice the amount, nor
> observe twice the sequence count relative to our GC content model for
> most of the genome.

> The average DNA abundance measured throughout the genome, whether it
> is raw sequence count, count ratio, or a logarithm of count ratio,
> does not increase in samples with a higher average copy number (more
> DNA per cell). Instead the average measured DNA abundance will be
> observed for the average copy number, whatever that is.

We are also likely to have some normal DNA content in the sample. If
both the tumor genome and normal genome are near diploid on average, and
the tumor cell content is near 50%, homozygous deletion (loss of both
homologous copies) would be observed at a coverge ratio near 0.5, as for
every tumor cell that would contribute zero copies to the DNA
extraction, one normal cell would contribute its normal two. As loss of
one copy in the tumor cell fraction would then appear at about 0.75 and
gains would appear near 1.25, 1.5, ..., this solution also appears to
fit our data reasonably well.

> **Question 3:** If a cancer genome has an average of 4 copies, and the
> tumor cell content is near 100%, at what count ratio should we ovserve
> a copy number of 3?

Fortunately there is another tool we can use to inform the copy number
estimates. We can use SNPs from a small variant caller and investigate
their variant allele ratios.

## Parse SNP allele ratio

Let's use the `VariantAnnotation` package to parse a VCF file containing
SNPs.

::: cell
``` {.r .cell-code}
library(VariantAnnotation)

# Parse the vcf file
vcf <- readVcf('sample1.vcf')
vcf

# Genotype column
g <- geno(vcf)
g

# Allelic and total read depth 
as.data.table(g$AD)
as.data.table(g$DP)
```
:::

The AD column of the VCF file contains the number of reads supporting
the reference and alternative allele. The DP column contains the total
read depth. Let's make a table of SNPs containing the alt-allele ratio.

::: cell
``` {.r .cell-code}
# A table of SNPs
snp_table <- data.table(id=names(vcf),
                        AD=sapply(g$AD, "[[", 2),
                        DP=unlist(g$DP[,1]))
snp_table

# Compute allele ratio
snp_table[,allele_ratio:=round(AD/DP,3)]
snp_table
```
:::

A GenomicRanges object defining the SNPs' positions on the reference
genome can be accessed with the `rowRanges` function. We can now overlap
that with the GenomicRanges representing our targets and use the
resulting `Hits` object to assign allele ratio to targets.

::: cell
``` {.r .cell-code}
# Match targets to SNPs
overlaps <- findOverlaps(target_ranges,rowRanges(vcf))
overlaps

# Set SNP allele ratio for targets
targets[queryHits(overlaps),allele_ratio:=snp_table$allele_ratio[subjectHits(overlaps)]]
targets
```
:::

Let's plot the SNP allele ratio with the count ratio and see if the copy
number status becomes more clear. Note that many targets contain no SNP,
and that many SNPs are homozygous. The allele ratio of a homozygous SNP
is either 0 or 1, depending on whether the alleles match the reference
genome or not. Only heterozygous SNP allele ratio is affected by copy
number alteration.

::: cell
``` {.r .cell-code}
p1 <- ggplot() + ylim(c(0,2.5)) +
    geom_point(data = targets, mapping = aes(x = target, y = count_ratio,fill=label),shape=21) + 
    geom_segment(data=segments,col='green',size=2,
                 mapping = aes(x=start,xend=end,y=count_ratio,yend=count_ratio))

p2 <- ggplot() +
    geom_point(data = targets, mapping = aes(x = target, y = allele_ratio,fill=label),shape=21)

p1 / p2 + plot_layout(guides = 'collect')
```
:::

We can see that for most segments, particularly segments with a count
ratio near 1, heterozygous SNPs have allele ratios near 0.5. This
indicates that the average copy number may be 2.

> **Question 4:** What do you think is the tumor cell fraction, and what
> has most likely happened to PTEN and BRCA2?

## Estimating genome-wide copy number

Although the ploidy and purity of this example are relatively clear,
that is not always the case. There are several
[methods](https://academic.oup.com/nar/article/44/16/e131/2460163)
available that fits your data to potential solutions of ploidy and
purity, and computes the copy numbers corresponding to the best fit.
Beware, this is somewhat prone to error, especially when the purity is
below 50%, if there are few copy number alterations, and if there is
some tumor cell heterogeneity. Also, reporting the total copy number of
every segment does not make sense in most settings. You should always
review the result and curate if necessary.

## Final questions

> Some questions require a bit of exploration and perhaps more code and
> visualisation. Motivate your answers to the extent you think is
> reasonable, and feel free to collaborate with the other students.

> **Question 5:** Homozygous deletions are rarely larger than one,
> occationally a few, megabases. What appears to be the sizes of the
> deletions affecting PTEN and BRCA2?

> **Question 6:** Let's assume you suspect that this patient may have a
> TMPRSS2-ERG fusion. Is that supported at all by copy number data?

> **Question 7:** Is there a MYC amplification?

> **Question 8:** How many copies of AR would you say there are in each
> tumor cell?

> **Question 9:** Let's assume there is a genomic segment with somatic
> copy-neutral loss of heterozogosity, i.e. one of the two homologous
> copies was lost and the other duplicated. What count ratio would we
> observe in this sample, and at what allele ratio(s) would the
> (germline-heterozygous) SNPs likely appear?

Once you have completed the exercise, make sure the resulting report
(generated by clicking **Render**) is a self-contained HTML file by
specifying [HTML format
options](https://quarto.org/docs/reference/formats/html.html) at the top
of your Quarto document:

    author: "Your Name"
    format:
      html:
        self-contained: true

You can add more options to make the result look better. Also check that
the report is readable, for example that it does not include full
print-outs of large data structures. When done, submit your rendered report
in Canvas under the Copy Number Analysis module.
