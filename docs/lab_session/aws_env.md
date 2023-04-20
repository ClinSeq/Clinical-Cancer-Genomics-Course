
# The workspace on AWS and installations

For the practical sessions each student will work partly using his/her own laptop and partly using a virtual machine on Amazon Web Services (AWS). This machine is an "instance" of a linux computer with the correct configurations. 

To be able to connect to AWS on windows you need to be able to download the SSH client [PuTTY](https://www.putty.org){:target="_blank"}. If you have a mac or a linux/unix computer nothing is needed to preinstall to work on the command line. 

To perform the exercises of the course you also need to install the following software on your local machine: [IGV](https://software.broadinstitute.org/software/igv/download){:target="_blank"}.
Choose IGV that comes bundled with Java for your operating system.
 
This workshop requires a large number of different bioinformatics tools. These have been pre-installed on AWS. For installation instructions of the differnet tools, see the list at the end of this document which links to the website of each software tool. 

## Explore the workspace on AWS
```bash
#Make sure the credentials are set for the "pem" (key) file.
# If needed run:
# chmod 400 course-setup-student-key.pem

# Log onto AWS cloud (make sure in are in the same directory as the pem key).
# You might get a question "Are you sure you want to continue connecting (yes/no/[fingerprint])?". Answer "yes".
ssh -i ./PEM_KEY_ID.pem ubuntu@AWS_ADDRESS_HERE

# Check where you are
pwd

# Note "~"" is a shortcut for the user directory. 
# Pre-load shortcut and settings 
source .bashrc

# the workspace folder contains the course files and pre-installed bioinformatic tools
cd workspace/

# The workspace folder contains the course files and pre-installed bioinformatic tools. Navigate to the folders, list the files and reflecton what is there
cd  ~/workspace/inputs
ls -alF
# Look around in the folders

cd ~/workspace/bin
ls -alF
# Look around in the folders

# Have a look at the raw exone sequencing data for the germline DNA
cd ~/workspace/inputs/data/fastq/Exome_Norm

# What are the files in this folder?

# Explore the file format (gunzip unpacks zipped files)
gunzip -c Exome_Norm_R1.fastq.gz | less -SN

# The format is nicely described here:
# https://en.wikipedia.org/wiki/FASTQ_format
# Read through the fastq format to get an understanding of
# base qualities from illumina DNA sequencing data.
# The format will be discussed tomorrow during the lab intro.
```

## R
[R](https://www.r-project.org){:target="_blank"} is a feature rich interpretive programming language originally released in 1995. It is heavily used in the bioinformatics community largely due to numerous R libraries available on [bioconductor](https://www.bioconductor.org){:target="_blank"}. 

If you have not installed R yet:

  1. Download R to your local machine
    * [R-project](https://www.r-project.org){:target="_blank"}
  2. Double-click on the icon and start the installation process. 

Install the following R libraries on your local machine:
```r
#Open R and run the following code:
packages <- c("devtools", "BiocManager", "tidyverse", "data.table", "patchwork", "jsonlite", "reshape2","cowplot")
install.packages(packages, dependencies = TRUE)
# hat-tip: The tidyverse is an opinionated collection of R packages designed for data science. 
# All packages share an underlying design philosophy, grammar, and data structures. 
# Installing Tidyverse will install several packages like tidyr, dplyr, stringr, ggplot2, lubridate, tibble, readr, purr and forcats.

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(pkgs=c("biomaRt", "Biostrings", "ensembldb", "IRanges", "EnsDb.Hsapiens.v86", "bamsignals", "BSgenome", "BSgenome.Hsapiens.NCBI.GRCh38", "csaw", "DNAcopy", "GenomicRanges", "org.Hs.eg.db", "Rsamtools", "Repitools", "TxDb.Hsapiens.UCSC.hg19.knownGene", "VariantAnnotation", "chimeraviz", "aroma.light"), dependencies = TRUE)

install.packages("PSCBS", dependencies = TRUE)
```
**Note:** 
  During the course you will be visualising data using R. If some of the libraries do not work on your local machine, you can run R and plot on AWS instead. However, then the plots will have to be dowloaded as pdf-files instead of being viewed interactively in R-studio.

The above mentioned libraries are installed in R on AWS, available in the following path, `/home/ubuntu/miniconda3/envs/r-rnaseq/bin/R`, to access this R env follow the below instructions.

- Using the command line do:
```bash
# Start R terminal without GUI
/home/ubuntu/miniconda3/envs/r-rnaseq/bin/R
```
- Using R do:
```r
# Load libraries. Note, not all are needed every time,
# just check here that they load without error. 
# Load each library and wait for it to finish before continuing 
library(BiocManager)
library(dplyr)
library(tidyr)
library(ggplot2)
library(IRanges)
library(data.table)
library(patchwork)
library(PSCBS)
library(stringr)
library(jsonlite)
library(reshape2)
library(cowplot)
library(biomaRt)
library(Biostrings)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(bamsignals)
library(BSgenome)
library(BSgenome.Hsapiens.NCBI.GRCh38)
#library(csaw)
library(DNAcopy)
library(GenomicRanges)
#library(org.Hs.eg.db)
library(Rsamtools)
#library(Repitools)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantAnnotation)
library(chimeraviz)

# Check your current directory
getwd()

# Set to your home directory or some other directory
setwd(~)
# here we can run the R scripts 

# Let us test to plot and save as pdf using pre-loaded data in R

# Check the data
mtcars

# Assign a plot to p
p <- ggplot(mtcars, aes(mpg, wt)) +
  geom_point()

# Save the p plot to the mtcars.pdf
ggsave(filename = "mtcars.pdf", plot = p)

# How to download the pdf-file to your local machine,
# see the section below, however first quit R. "no"
# means that the R-workspace will not be saved.
q("no")
```
## How to use the terminal download files to your local machine
This can be used for both downloading files plotted in R on AWS and other files and file types, such as bam-files that will be visualised in IGV.

```bash
cd TO_PREFERRED_DIRECTORY
scp -ri ~/PEM_KEY_ID.pem ubuntu@AWS_ADDRESS_HERE:~/PATH_TO_FILE/FILENAME .

# Before proceeding - test that you can get the "mtcars.pdf"
# file available here: /home/ubuntu/workspace/
# downloaded to your local computer. If not, ask for assistance.
```

## List of preinstalled bioinformatic tools

### bam-readcount

[bam-readcount](https://github.com/genome/bam-readcount){:target="_blank"} is a program for determing read support for individual variants (SNVs and Indels only). 

### BCFtools

[BCFtools](https://samtools.github.io/bcftools/){:target="_blank"} is an open source program for variant calling and manipulating files in Variant Call Format (VCF) or Binary Variant Call Format (BCF). 

### BWA

[BWA](http://bio-bwa.sourceforge.net){:target="_blank"} is a popular DNA alignment tool used for mapping sequences to a reference genome. It is available under an open source GPLv3 license.

### CNVkit

[CNVkit](https://cnvkit.readthedocs.io/en/stable/){:target="_blank"} is a python based copy number caller designed for use with hybrid capture. It will not be applied during the course but is a frequently used tool for CNV-analysis of copy-number data.

### Delly

[Delly](https://github.com/dellytools/delly){:target="_blank"} is a structural variant caller developed at EMBL. It uses paired-ends, split-reads and read-depth to sensitively and accurately delineate genomic rearrangements throughout the genome. 

### FastQC

[FastQC](https://github.com/s-andrews/FastQC){:target="_blank"} is a quality control program for raw sequencing data. 

### GATK 4

[GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4){:target="_blank"} is a toolkit developed by the broad institute focused primarily on variant discovery and genotyping. It is open source, hosted on github, and available under a BSD 3-clause license. 

### Gffcompare

[Gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml){:target="_blank"} is a program that is used to perform operations on general feature format (GFF) and general transfer format (GTF) files. 

### HISAT2

[HISAT2](http://daehwankimlab.github.io/hisat2/){:target="_blank"} is a graph based alignment algorithm devoloped at Johns Hopkins University. It is heavily used in the bioinformatics community for RNAseq based alignments. 

### Kallisto

[Kallisto](https://pachterlab.github.io/kallisto/){:target="_blank"} is a kmer-based alignment algorithm used for quantifying transcripts in RNAseq data.

### mosdepth

[mosdepth](https://github.com/brentp/mosdepth){:target="_blank"} is a program for determining depth in sequencing data. 

### MultiQC

[MultiQC](https://github.com/ewels/MultiQC){:target="_blank"} is a tool to create a single report with interactive plots for multiple bioinformatics analyses across many samples.

### PICARD

[PICARD](http://broadinstitute.github.io/picard/){:target="_blank"} is a set of java based tools developed by the Broad institute. It is very useful for manipulating next generation sequencing data and is available under an open source MIT license. 

- Try to answer the following:
    - Who is Picard?
        - Clue: "TNG"
### Pizzly

[Pizzly](https://github.com/pmelsted/pizzly){:target="_blank"} is a fusion detection algorithm which uses output from Kallisto. 

### Sambamba

[Sambamba](https://lomereiter.github.io/sambamba/){:target="_blank"} is a high performance alternative to samtools and provides a subset of samtools functionality. 

### Samtools

[Samtools](https://samtools.github.io) is a software package based in C which provies utilities for manipulating alignment files (SAM/BAM/CRAM). It is open source, available on github, and is under an MIT license. 

### seqtk

[Seqtk](https://github.com/lh3/seqtk){:target="_blank"} is a lighweight tool for processing FASTQ and FASTA files. We will use seqtk to subset RNA-seq fastq files to more quickly run the fusion alignment module.

### Strelka

[Strelka](https://github.com/Illumina/strelka){:target="_blank"} is a germline and somatic variant caller developed by illumina and available under an open source GPLv3 license. 

### StringTie

[StringTie](https://ccb.jhu.edu/software/stringtie/){:target="_blank"} is a software program to perform transcript assembly and quantification of RNAseq data. 

### Varscan

[Varscan](http://dkoboldt.github.io/varscan/){:target="_blank"} is a java program designed to call variants in sequencing data. It was developed at the Genome Institute at Washington University and is hosted 
on github. 

### vcf-annotation-tools

[VCF Annotation Tools](https://github.com/griffithlab/VAtools){:target="_blank"} is a python package that includes several tools to annotate VCF files with data from other tools. We will be using this for the purpose of adding bam readcounts to the vcf files.

### VEP

[VEP](https://www.ensembl.org/info/docs/tools/vep/index.html){:target="_blank"} is a variant annotation tool developed by ensembl and written in perl. By default VEP will perform annotations by making web-based API queries however it is much faster to have a local copy of cache and fasta files. The AWS AMI image we’re using already has these files for hg38 in the directory ~/workspace/vep_cache as they can take a bit of time to download. 

### vt

[vt](https://github.com/atks/vt){:target="_blank"} is a variant tool set that discovers short variants from Next Generation Sequencing data. We will use this for the purpose of splitting multi-allelic variants.
