
# The workspace on KI and installations

Due to the computationally intensive nature of this course, we have transitioned the practical sessions to the KI VM server. Each student will be provided access to a virtual machine (VM) on the KI server to perform the lab work. This setup ensures that all required bioinformatics tools are pre-installed and ready to use, eliminating the need for local installations and providing a consistent environment for all participants.

To connect to the KI VM server, you will need an SSH client. If you are using Windows, you can download [PuTTY](https://www.putty.org). For macOS or Linux/Unix users, SSH is available by default in the terminal.

During the lab sessions, our tutors and teaching staffs will be available for discussions and support. To perform the exercises of the course, you will need to install [IGV](https://software.broadinstitute.org/software/igv/download) on your local machine. Choose the IGV version bundled with Java for your operating system.

The KI VM server provides a robust computational environment with all necessary bioinformatics tools pre-installed. For more details about the tools and their usage, refer to the list at the end of this document.

## How to access Bioinformatics Tools
To access the bioinformatics tools on the server, you need to configure your environment to use Conda. Conda is an open-source package and environment management system that simplifies the installation and management of software dependencies. A faster alternative to Conda is Mamba, which is also available on the server.

When you log in to the server for the first time, execute the following commands to set up your environment:

```bash
echo 'export PATH=/nfs/course/miniforge3/bin:$PATH' >> ~/.bashrc
echo 'source activate base' >> ~/.bashrc
source ~/.bashrc
```

There are three Conda environments available on the server:
1. `base`: The default Conda environment.
2. `bioinfo-tools`: Contains most of the bioinformatics tools required for the course.
3. `ensembl-vep-113`: Specifically configured for Ensembl VEP (Variant Effect Predictor) version 113.

To activate the `bioinfo-tools` environment, use the following command:

```bash
source activate bioinfo-tools
```

Once activated, all the tools in the `bioinfo-tools` environment will be available in your current session. You can verify the active environment by running:

```bash
# try to access any tools - bwa
bwa 
samtools --help
```

### About Conda and Mamba

- **Conda**: A powerful package manager that allows you to create isolated environments for different projects, ensuring compatibility and avoiding dependency conflicts.
- **Mamba**: A faster implementation of Conda, written in C++. It significantly speeds up environment creation and package installation.

> **Important Note**: Do not attempt to install or modify any software on the server. The environments are pre-configured to ensure consistency and stability for all users. If you encounter any issues or require additional tools, please contact the course administrators for assistance.

This setup ensures a consistent and efficient environment for running bioinformatics tools during the course.


## R
[R](https://www.r-project.org) is a feature rich interpretive programming language originally released in 1995. It is heavily used in the bioinformatics community largely due to numerous R libraries available on [bioconductor](https://www.bioconductor.org). 

If you have not installed R yet:
- Download R to your local machine:
    - https://www.r-project.org
- Double-click on the icon and start the installation process. 

Install the following R libraries on your local machine:
```r
#Open R and run the following code:
packages <- c("devtools", "BiocManager", "dplyr", "tidyr", "ggplot2", "data.table", "patchwork", "stringr", "jsonlite", "reshape2","cowplot")
install.packages(packages, dependencies = TRUE)

BiocManager::install("biomaRt", "Biostrings", "ensembldb", "IRanges", "EnsDb.Hsapiens.v86", "bamsignals", "BSgenome", "BSgenome.Hsapiens.NCBI.hg19", "csaw", "DNAcopy", "GenomicRanges", "org.Hs.eg.db", "Rsamtools", "Repitools", "TxDb.Hsapiens.UCSC.hg19.knownGene", "VariantAnnotation", "chimeraviz" dependencies = TRUE)

# PSCBS requires DNAcopy 
install.packages("PSCBS")
```
Note during the course you will be visualising data using R. If some of the libraries do not work on your local machine, you can run R and plot on KI instead. However, then the plots will have to be dowloaded as pdf-files instead of being viewed interactively in R-studio.


## How to use the terminal download files to your local machine

This will can be used for both downloading files plotted in R on AWS and other files and file types, such as bam-files that will be visualised in IGV.

```bash
cd TO_PREFERRED_DIRECTORY
scp ki-username@c8cancergen02.ki.se:/PATH_TO_PDF/FILENAME .

# Try to download some test file to your local computer. If not, ask for assistance.
```

## Explore the workspace on KI VM
```bash

ssh ki-username@c8cancergen02.ki.se
 
# Check where you are
pwd

# Note "~"" is a shortcut for the user directory. 
# Pre-load shortcut and settings 
echo 'export PATH=/nfs/course/miniforge3/bin:$PATH' >> ~/.bashrc
echo 'source activate base' >> ~/.bashrc
source .bashrc

# the KI /nfs/course/ contains the course files and pre-installed bioinformatic tools
cd /nfs/course/

# The /nfs/course/ folder contains the course files and pre-installed bioinformatic tools. Navigate to the folders, list the files and reflecton what is there
cd  /nfs/course/inputs/
ls -alF
# Look around in the folders

# Have a look at the raw exone sequencing data for the germline DNA
cd /nfs/course/inputs/data/fastq/Exome_Norm

# What are the files in this folder?

# Explore the file format (gunzip unpacks zipped files)
gunzip -c Exome_Norm_R1.fastq.gz | less -SN

# The format is nicely described here:
# https://en.wikipedia.org/wiki/FASTQ_format
# Read through the fastq format to get an understanding of
# base qualities from illumina DNA sequencing data.
# The format will be discussed tomorrow during the lab intro.

## Each user has their own working directory
# NOTE: Your results, outputs, everything should be written inside your
# working directory
cd /nfs/course/students/KI-USERNAME
```

## List of preinstalled bioinformatic tools

### BCFtools

[BCFtools](https://samtools.github.io/bcftools/) is an open source program for variant calling and manipulating files in Variant Call Format (VCF) or Binary Variant Call Format (BCF). 

### BWA

[BWA](http://bio-bwa.sourceforge.net) is a popular DNA alignment tool used for mapping sequences to a reference genome. It is available under an open source GPLv3 license.

### Delly

[Delly](https://github.com/dellytools/delly) is a structural variant caller developed at EMBL. It uses paired-ends, split-reads and read-depth to sensitively and accurately delineate genomic rearrangements throughout the genome. 

### FastQC

[FastQC](https://github.com/s-andrews/FastQC) is a quality control program for raw sequencing data. 

### GATK 4

[GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4) is a toolkit developed by the broad institute focused primarily on variant discovery and genotyping. It is open source, hosted on github, and available under a BSD 3-clause license. 

### Gffcompare

[Gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) is a program that is used to perform operations on general feature format (GFF) and general transfer format (GTF) files. 

### HISAT2

[HISAT2](http://daehwankimlab.github.io/hisat2/) is a graph based alignment algorithm devoloped at Johns Hopkins University. It is heavily used in the bioinformatics community for RNAseq based alignments. 

### Kallisto

[Kallisto](https://pachterlab.github.io/kallisto/) is a kmer-based alignment algorithm used for quantifying transcripts in RNAseq data.


### MultiQC

[MultiQC](https://github.com/ewels/MultiQC) is a tool to create a single report with interactive plots for multiple bioinformatics analyses across many samples.

### PICARD

[PICARD](http://broadinstitute.github.io/picard/) is a set of java based tools developed by the Broad institute. It is very useful for manipulating next generation sequencing data and is available under an open source MIT license. 

- Try to answer the following:
    - Who is Picard?
        - Clue: "TNG"

### Pizzly

[Pizzly](https://github.com/pmelsted/pizzly) is a fusion detection algorithm which uses output from Kallisto. 

### Sambamba

[Sambamba](https://lomereiter.github.io/sambamba/) is a high performance alternative to samtools and provides a subset of samtools functionality. 

### Samtools

[Samtools](https://samtools.github.io) is a software package based in C which provies utilities for manipulating alignment files (SAM/BAM/CRAM). It is open source, available on github, and is under an MIT license. 

### seqtk

[Seqtk](https://github.com/lh3/seqtk) is a lighweight tool for processing FASTQ and FASTA files. We will use seqtk to subset RNA-seq fastq files to more quickly run the fusion alignment module. 

### StringTie

[StringTie](https://ccb.jhu.edu/software/stringtie/) is a software program to perform transcript assembly and quantification of RNAseq data. 


### vcf-annotation-tools

[VCF Annotation Tools](https://github.com/griffithlab/VAtools) is a python package that includes several tools to annotate VCF files with data from other tools. We will be using this for the purpose of adding bam readcounts to the vcf files.

### VEP

[VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) is a variant annotation tool developed by ensembl and written in perl. By default VEP will perform annotations by making web-based API queries however it is much faster to have a local copy of cache and fasta files. The AWS AMI image we’re using already has these files for hg38 in the directory ~/workspace/vep_cache as they can take a bit of time to download. 

### vt

[vt](https://github.com/atks/vt) is a variant tool set that discovers short variants from Next Generation Sequencing data. We will use this for the purpose of splitting multi-allelic variants.

### SAGE

[SAGE](https://github.com/hartwigmedical/hmftools/tree/master/sage) is a precise and highly sensitive somatic SNV, MNV and small INDEL caller. It has dynamically scaling sensitivity based on the depth of the provided tumor and germline BAMs, but performs best if both BAMs have at least 30x typical depth.