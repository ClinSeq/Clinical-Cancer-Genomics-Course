# Analysis of copy number variation

### Copy number analysis tools
There are plenty of tools available for copy number analysis today, and you will likely be able to find one for free that suits your particular application. Although tailored for different types of sequence data, general processing steps of all copy number analysis tools include:
1. Quantification of sequence read depth throughout the reference genome as a measure of DNA abundance in the sample
1. Removal of sample-specific systematic noise using features such as GC content and mappability
1. Removal of assay-specific systematic noise using normal (non-cancer) reference samples
1. Segmentation - partitioning of the reference genome so that each segment can be assigned one copy number
1. Copy number calling, assigning each segment some estimate of the number of copies per cell
1. Combine normalized coverage with SNP allele frequencies to support copy number estimates

In this exercise we will perform the above steps using R, and a BAM file of aligned reads and a VCF file of SNPs as input.


### Download data and prepare R

The data used in this exercise is available for download here: [copy_number_files.tar.bz2](https://course-cg-5534.s3.amazonaws.com/cnvs/copy_number_files.tar.bz2)

Extract the content into a local folder on your computer.

```bash
tar xjvf copy_number_files.tar.bz2
```

Several packages from the [Bioconductor](https://www.bioconductor.org/) repository are used in this exercise, most importantly **GenomicRanges** which facilitates working with anything located on the reference genome. [Here](http://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesHOWTOs.pdf) is a good overview.


### Parse sequence target definition
We use targeted sequencing of a prostate cancer sample in this exercise, as the resulting data files are of a more convenient size than whole genome or exome. This hybrid-capture panel covers a few hundred cancer genes and a few thousands of additional targets intended to improve copy number and structural variant detection.

In RStudio, navigate to the folder containing the exercise files:
```r
setwd('your/path/to/exercise_folder')

```

> Tip: Open a new Rmarkdown file and save it in the exercise data folder. Keep your code as shown in the template, and write your comments and answers to questions between code chunks. Click **knit** any time to run all code and generate a pdf or html report. This way, when you are through this exercise, your report is done too.

The **BED** (Browser Extensible Data) format is a text file format used to store genomic regions as coordinates and associated annotations. The data are presented in the form of columns separated by spaces or tabs. This format was developed during the Human Genome Project and then adopted by other sequencing projects.

Parse and investigate the target definition BED file.
```r
library(data.table)

bed_file <- 'targets.bed'
file.exists(bed_file)

targets <- fread(bed_file)
targets

# Let the first column be the target number
targets <- cbind(1:nrow(targets),targets)

# Set useful column names
colnames(targets) <- c('target','chromosome','start','end')

# Check target length
targets[,length:=end-start]
summary(targets)

# Check target distribution over chromosomes
targets[,table(chromosome)]

```

> The R code in this exercise features some **data.table** syntax. If you are more familiar with base R or tidyverse syntax, a comparison between them can be found [here](https://mgimond.github.io/rug_2019_12/Index.html).

### Annotate with gene symbol
Before parsing the sequence data, let's assign gene symbols to targets that overlap a gene. We'll take advantage of databases accessible through Bioconductor packages. First, we need to create a GenomicRanges object representing our targets:

```r
library(GenomicRanges)

target_ranges <- makeGRangesFromDataFrame(targets)
target_ranges

seqlevelsStyle(target_ranges)

```

As you may recognize, chromosomes have the NCBI/EMBL sequence naming style (1,...). This matches the reference genome used with our sequence data. Let's create a second GenomicRanges object featuring the UCSC names (chr1,...) for use with some database queries.

```r
ucsc_ranges <- target_ranges
seqlevelsStyle(ucsc_ranges) <- "UCSC"

ucsc_ranges

```

We then use the **detailRanges** function from the **csaw** package to overlap our target ranges with known genes.

```r
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(csaw)

d <- detailRanges(ucsc_ranges, orgdb=org.Hs.eg.db,
    txdb=TxDb.Hsapiens.UCSC.hg19.knownGene)
d

as.data.table(d)

```

> There are many ways to achieve the same thing. You should focus mostly on **what** we do, and not worry too much about exactly how. Google is your friend and if you know what you want do do, you can always find a way.

We only care about some level of overlap at this point, let's just retrieve the gene symbols and add as a column to the table of targets. We can then see how targets allocate to genes.

```r
library(stringr)

targets$gene <- str_remove(string = d$overlap,pattern = ':.*')
targets

table(targets$gene)

as.data.table(table(targets$gene))

as.data.table(table(targets$gene))[N>25]

```


### Parse sequence read counts
We are now ready to access the BAM file and count the number of read pairs (i.e. fragments) that map to each target. For this we use the [bamsignals](https://bioconductor.org/packages/release/bioc/vignettes/bamsignals/inst/doc/bamsignals.html) package and the GenomicRanges object with NCBI-style sequence names. To count a read pair, we require it's midpoint to fall within a target, and that it is not flagged as a [duplicate](http://samtools.github.io/hts-specs/SAMv1.pdf).

```r
library(bamsignals)
library(Rsamtools)

tumor_bam <- 'Sample1.bam'
file.exists(tumor_bam)
file.exists(paste0(tumor_bam,'.bai'))

targets$coverage <- bamCount(tumor_bam, target_ranges, paired.end="midpoint", filteredFlag=1024, verbose=F)
targets

summary(targets)

```

> Strictly speaking, read count is not equal to sequence coverage. But given the fragment length and target size in this example, we can safely call the result sequence coverage.

Zero coverage will be a problem later as it log-transforms to -Inf. We set those to 1:

```r
targets[coverage==0, coverage:=1]

```

### Plot sequence coverage
Let's investigate the raw sequence coverage across targets. You are probably already familiar with [ggplot2](https://ggplot2.tidyverse.org/).
```r
library(ggplot2); theme_set(theme_bw())

# Simple ggplot
ggplot(data = targets) + geom_point(mapping = aes(x = target, y = coverage))

# Some adjustments
ggplot(data = targets) + ylim(c(0,3500)) +
  geom_point(mapping = aes(x = target, y = coverage),size=.2)

# To use start position, we can separate the plot by chromosome.
ggplot(data = targets) + ylim(c(0,3500)) +
  geom_point(mapping = aes(x = start, y = coverage),size=.5) + 
  facet_wrap(facets = vars(chromosome),ncol = 2)

```

The figure may now look squished in RStudio's bottom-right plots pane. Click the **Zoom** button and maximize the new window.

Choose a few genes to highlight in the plot.

```r
myGenes <- c('AR','ATM','BRCA1','BRCA2','PTEN','TMPRSS2','ERG')

targets$label=''
targets[gene %in% myGenes,label:=gene][label=='',label:=NA]


ggplot() + ylim(c(0,3500)) +
  geom_point(data = targets, mapping = aes(x = target, y = coverage, col = label), size=.2)

ggplot(data = targets) + ylim(c(0,3500)) +
  geom_point(mapping = aes(x = start, y = coverage, col = label), size=.5) +
  facet_wrap(facets = vars(chromosome),ncol = 2)

```

### Model and correct for GC content bias
Sequence GC content is known to affect PCR amplification and sequence coverage. Let's retrieve the GC content for the targets and investigate if there is a bias in our data set.

```r
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome)
library(Repitools)

targets$gc <- gcContentCalc(ucsc_ranges , organism=Hsapiens)

ggplot(data = targets) +
  geom_point(mapping = aes(x = gc, y = coverage),alpha=.2)

```

There seems to be a modest effect of target GC content on sequence coverage. Let's build a [loess](https://en.wikipedia.org/wiki/Local_regression) model of the effect. We leave the X chromosome out.

```r
loess_tumor <- loess(coverage ~ gc, data = targets, 
                    subset = chromosome!='X', family="symmetric",
                    control = loess.control(surface = "direct"))

```

To plot the model, we predict the sequence coverage for GC content ranging from 1% to 100%:

```r
tumor_gc_line <- data.table(x=1:100/100,
                        y=predict(loess_tumor,data.table(gc=1:100/100)))
 
ggplot(data=targets) +
  geom_point(mapping = aes(x = gc, y = coverage),alpha=.2) +
  geom_line(data=tumor_gc_line, mapping = aes(x=x,y=y),col='blue')

```

We seem to predict negative coverage at some high GC content. We'll try the same model on log2 of coverage. It will now return the predicted log2 coverage, which we can either keep using, or transform back to coverage. We'll keep using coverage for now.

```r
loess_tumor <- loess(log2(coverage) ~ gc, data = targets, 
                    subset = chromosome!='X', family="symmetric",
                    control = loess.control(surface = "direct"))

# prediction is log2(coverage), 2^prediction equals predicted coverage
tumor_gc_line <- data.table(x=1:100/100,
                        y=2^predict(loess_tumor,data.table(gc=1:100/100)))

ggplot(data=targets) + ylim(c(0,3500)) +
  geom_point(mapping = aes(x = gc, y = coverage),alpha=.2) +
  geom_line(data=tumor_gc_line, mapping = aes(x=x,y=y),col='blue')

```

To adjust for GC content bias, we divide the observed coverage with the predicted (using targets' GC content). This new **coverage ratio** should be a slightly better measure of DNA abundance.

```r
targets[,coverage_ratio:=coverage/2^predict(loess_tumor,gc)]
targets

```

To compare the new metric with the previous, we can plot them side-by-side. With [patchwork](https://patchwork.data-imaginist.com/) we can use arithmetic operators to combine and align multiple plots.

```r
library(patchwork)

p1 <- ggplot(data = targets) + 
    geom_point(mapping = aes(x = gc, y = coverage),alpha=.05)

p2 <- ggplot(data = targets) + 
    geom_point(mapping = aes(x = gc, y = coverage_ratio),alpha=.05)

p1+p2

p1 <- ggplot(data = targets) + ylim(c(0,3500)) +
    geom_point(mapping = aes(x = target, y = coverage, fill=label),shape=21)

p2 <- ggplot(data = targets) + ylim(c(0,2.5)) +
    geom_point(mapping = aes(x = target, y = coverage_ratio, fill=label),shape=21)

p1/p2 + plot_layout(guides = 'collect')

p1 <- ggplot(data = targets) + ylim(c(0,3500)) +
    geom_point(mapping = aes(x = start, y = coverage, fill=label),shape=21) +
    facet_wrap(facets = vars(chromosome),ncol = 2) +
    theme(panel.spacing = unit(0, "lines"),strip.text.x = element_text(size = 4))

p2 <- ggplot(data = targets) + ylim(c(0,2.5)) +
      geom_point(mapping = aes(x = start, y = coverage_ratio, fill=label),shape=21) +
      facet_wrap(facets = vars(chromosome),ncol = 2) +
      theme(panel.spacing = unit(0, "lines"),strip.text.x = element_text(size = 4))

p1+p2 + plot_layout(guides = 'collect')

```

The GC content correction seems to have had a small but positive effect. It is common to also use either a matched normal or a pool of normal reference samples (sequenced similarly) to remove some assay-specific noise. We have omitted that step here.

> One important reason to perform copy number analysis of cancer samples is to find homozygous deletions of tumor suppressor genes. If a patient has a hemizygous germline deletion (loss of one copy) affecting a tumor suppressor gene, and somatic deletion of the healthy copy in their tumor, would the resulting homozygous deletion still be visible if the coverage ratio of the tumor sample was divided by that of the normal sample?

### Segment the data
Although we can already spot some apparent gains and losses, statistical tools can help us better estimate copy number segments, for which we can then calculate the most likely copy number given the observation. **Circular binary segmentation** (CBS) is probably the most commonly used segmentation method. Here we use the **PSCBS** ("Parent specific" CBS, as it can also use SNP allele data) R package as a wrapper to perform basic CBS on the GC content-adjusted coverage ratio. 

CBS requires the DNA abundance ratio to be log-transformed, making its distribution more normal-like. Conveniently we have no zeroes in the data (that would become -Inf). 

```r
library(PSCBS)
targets[,log_ratio:=log2(coverage_ratio)]
targets

segments <- segmentByCBS(y=targets$log_ratio)
segments

```

After segmentation we can transform the segment mean values back and plot the segments with the targets.

```r
segments <- as.data.table(segments)[,coverage_ratio:=2^mean]
segments

ggplot() + ylim(c(0,2.5)) +
  geom_point(data = targets, mapping = aes(x = target, y = coverage_ratio,fill=chromosome),shape=21) + 
  geom_segment(data=segments,col='green',size=2,
               mapping = aes(x=start,xend=end,y=coverage_ratio,yend=coverage_ratio))

```

Note that the segmentation is based only on the target log-ratio and order, and that segment start and end refer to target number. To avoid segments spanning across chromosomes, and breakpoints just missing the chromosome boundary, we can add chromosomes to the segmentation call. The chromosome vector is required to be numeric.

```r
segments <- segmentByCBS(y=targets$log_ratio,
                      chromosome=as.numeric(str_replace(targets$chromosome,'X','23')))
segments

segments <- as.data.table(segments)[,coverage_ratio:=2^mean][!is.na(chromosome),-1]
segments

ggplot() + ylim(c(0,2.5)) +
  geom_point(data = targets, mapping = aes(x = target, y = coverage_ratio,fill=chromosome),shape=21) + 
  geom_segment(data=segments,col='green',size=2,
               mapping = aes(x=start,xend=end,y=coverage_ratio,yend=coverage_ratio))
               
```

CBS can also take a vector of chromosomal positions as input, in which case the resulting segment start and end positions are based on them. We use the target midpoints:

```r
segments_pos <- segmentByCBS(y=targets$log_ratio,
                         chromosome=as.numeric(str_replace(targets$chromosome,'X','23')),
                         x=targets$start+60)
                         
segments_pos <- as.data.table(segments_pos)[,coverage_ratio:=2^mean][!is.na(chromosome),-1]
segments_pos


# convert chromosomes back to NCBI
segments_pos[,chromosome:=str_replace(as.character(chromosome),'23','X')]

ggplot(data = targets) + ylim(c(0,2.5)) +
  geom_point(mapping = aes(x = start, y = coverage_ratio, fill=label),shape=21) +
  geom_segment(data=segments_pos,col='green',
               mapping = aes(x=start,xend=end,y=coverage_ratio,yend=coverage_ratio)) +
  facet_wrap(facets = vars(chromosome),ncol = 2) +
  theme(panel.spacing = unit(0, "lines"),strip.text.x = element_text(size = 6))

```

> Target density is only about one per megabase, with much higher density in some cancer genes. How many targets do you think a deletion would have to cover for us to be able to find it? What else might influence sensitivity?

Let's now take another look at segmented targets plotted in order, as they are very unevenly distributed over the genome. We can easily spot some deletions affecting 2-3 genes in our selection, as well as amplification of another. Take a closer look at BRCA2 and part of PTEN. 

```r             
p1 <- ggplot() + ylim(c(0,2.5)) + 
  geom_point(data = targets, mapping = aes(x = target, y = coverage_ratio,fill=chromosome),shape=21) + 
  geom_segment(data=segments,col='green',size=2,
               mapping = aes(x=start,xend=end,y=coverage_ratio,yend=coverage_ratio))

p2 <- ggplot() + ylim(c(0,2.5)) + 
  geom_point(data = targets, mapping = aes(x = target, y = coverage_ratio,fill=label),shape=21) + 
  geom_segment(data=segments,col='green',size=2,
               mapping = aes(x=start,xend=end,y=coverage_ratio,yend=coverage_ratio))

p1/p2

```


Coverage ratio (observed sequence coverage relative to the expected, given target GC content), is near 0.5, or 50%. For a largely diploid genome this is exactly what to expect if these segments were hemizygously deleted, i.e., one out of the two homologous copies would have been lost. This would also fit reasonably well with gain of 1 copy of 8q (resulting in a coverage ratio of 1.5) and 3 extra copies of AR (being on the X chromosome, tumor cells would have quadrupled its single original AR copy, resulting in the coverage ratio going from about 0.5 to near 2.).

But this solution would require the tumor cell content ("**purity**") in the sample to be near 100%, which is rare unless the sample comes from a cell line. Note also that several other segments have a mean coverage ratio that would not fit this solution well. There are actually multiple other combinations of copy number and purity that would fit the observation reasonably well. 

Our measure of relative DNA abundance - coverage ratio - is automatically centered around 1, regardless of the average copy number ("**ploidy**") of the sequenced genome. Let's assume the cancer cells had twice the number of copies of all chromosomes. This would not double the average sequence coverage of the sample, and we would also not observe twice the sequence coverage relative to our GC content model for most of the genome. 

> We have analyzed a certain amount of DNA, not a certain number of cells. The average DNA abundance measured throughout the genome, whether it is raw sequence coverage, coverage ratio, or log ratio, does not increase in samples with a higher average copy number (more DNA per cell). Instead the average measured DNA abundance will be observed at the average copy number, whatever that is. Does that make sense?

We are also likely to have some normal DNA content in the sample. If both the tumor genome and normal genome are near diploid on average, and the tumor cell content is near 50%, homozygous deletion (loss of both homologous copies) would be observed at a coverge ratio near 0.5, as for every tumor cell that would contribute zero copies to the DNA extraction, one normal cell would contribute its normal two. As loss of one copy in the tumor cell fraction would then appear at about 0.75 and gains would appear near 1.25, 1.5, ..., this solution also appears to fit our data reasonably well.

Fortunately there is another tool we can use to inform our copy number estimates. We can use SNPs from a small variant caller and investigate their variant allele ratios.

### Parse SNP allele ratio

Let's use the VariantAnnotation package to parse a VCF file containing SNPs.

```r
library(VariantAnnotation)

vcf <- readVcf('Sample1.vcf')
vcf

g <- geno(vcf)
g

g$AD

as.data.table(g$AD)

as.data.table(g$DP)

```

The AD column of the VCF file contains the number of reads supporting the reference and alternative allele. The DP column contains the total read depth. Let's make a table of SNPs containing the alt-allele coverage ratio.

```r
snp_table <- data.table(id=names(vcf),
                        AD=sapply(g$AD, "[[", 2),
                        DP=unlist(g$DP[,1]))
snp_table

snp_table[,allele_ratio:=round(AD/DP,3)]
snp_table

```

A GenomicRanges object defining the SNPs' positions on the reference genome can be accessed with the rowRanges() function. We can now overlap that with the GenomicRanges representing our targets and use the resulting Hits object to assign allele ratio to targets.

```r
overlaps <- findOverlaps(target_ranges,rowRanges(vcf))
overlaps

targets$allele_ratio <- NA
targets$allele_ratio[queryHits(overlaps)] <- snp_table$allele_ratio[subjectHits(overlaps)]
targets

```

Let's plot the SNP allele ratio with the coverage ratio and see if the copy number status becomes more clear. Note that many targets contain no SNP, and that many SNPs are homozygous. The allele ratio of a homozygous SNP is 0 or 1 and unaffected by copy number alteration.

```r
p1 <- ggplot() + ylim(c(0,2.5)) +
  geom_point(data = targets, mapping = aes(x = target, y = coverage_ratio,fill=label),shape=21) + 
  geom_segment(data=segments,col='green',size=2,
               mapping = aes(x=start,xend=end,y=coverage_ratio,yend=coverage_ratio))
               
p2 <- ggplot() +
  geom_point(data = targets, mapping = aes(x = target, y = allele_ratio,fill=label),shape=21)

p1 / p2 + plot_layout(guides = 'collect')

```

We can see that for most segments, particularly segments with a coverage ratio near 1, heterozygous SNPs have allele ratios near 0.5. This indicates that the average copy number may be 2. 

> What do you think is the tumor cell content, and what has most likely happened to PTEN and BRCA2? 

### Estimating genome-wide copy number

Although the ploidy and purity of this example are relatively clear, that is not always the case. There are several [methods](https://academic.oup.com/nar/article/44/16/e131/2460163) available that fits your data to potential ploidies and purities, and assigns the copy numbers corresponding to the best fit. Beware, this is somewhat prone to error, especially when the purity is below 50%, if there are few copy number alterations, or if there is some tumor cell heterogeneity. Also, reporting the total copy number of every segment does not always make sense. You should always review the result and curate if necessary. 

### Investigate further

> Most copy number analysis tools use log-transformed coverage ratio (log-ratio). What may be the advantages of this?

> Homozygous deletions are rarely larger than one, occationally a few, megabases. What appears to be the sizes of the deletions affecting PTEN and BRCA2?

> Let's assume you suspect that this patient may have a TMPRSS2-ERG fusion. Is that supported by copy number data?

> Is there a MYC amplification?

> How many copies of AR would you say there are in each tumor cell?

> Let's assume we have a segment with somatic copy-neutral loss of heterozogosity, i.e. one of the two homologous copies was lost and the other duplicated. What coverage ratio would we observe in this sample, and at what allele ratio(s) would the (germline-heterozygous) SNPs likely appear?