
# THE WORKSPACE DIRECTORIES ON AWS

This workshop requires a large number of different bioinformatics tools. The instructions for installing these tools will vary per operating system. It is beyond the scope of this course to provide installation instructions for all types of operating systems used in this course. Therefore, we have pre-installed the tools on Amazon Web Services (AWS) and will provide access to virtual machines for the lab sessions. 

The only software that needs to be installed on your local machine is "R" to be able to visualize data during the labwork.

## R

[R](https://www.r-project.org) is a feature rich interpretive programming language originally released in 1995. It is heavily used in the bioinformatics community largely due to numerous R libraries available on [bioconductor](https://www.bioconductor.org).

- Download R to your local machine:
    - https://www.r-project.org
- Double-click on the icon and start the installation process. 
- Open R. 
    - In R run othe following lines of code:
        - packages <- c("devtools", "BiocManager", "dplyr", "tidyr", "ggplot2", "IRanges")
        - install.packages(packages, dependencies = TRUE)

## Explore the workspace on AWS

```bash
# Log onto AWS cloud
ssh -i course_EC2_01.pem ubuntu@AWS_ADDRESS_HERE

# Check where you are
pwd

# Note "~"" is a shortcut for the user directory. 
# Pre-load shortcut and settings 
source .bashrc

# the workspace folder contains the course files and pre-installed bioinformatic tools
cd workspace/

# the workspace folder contains the course files and pre-installed bioinformatic tools. Navigate to the folders, list the files and reflecton what is there
cd  ~/workspace/inputs
ls -alF
#Look around in the folders

cd ~/workspace/bin
ls -alF
#Look around in the folders

#Have a look at the raw exone sequencing data for the germline DNA
cd ~/workspace/inputs/data/fastq/Exome_Norm

#What are the files in this folder?

#Explore the file format (gunzip unpacks zipped files)
gunzip -c Exome_Norm_R1.fastq.gz | less -SN
#The format is nicely described here:
#https://en.wikipedia.org/wiki/FASTQ_format
```

## List of preinstalled bioinformatic tools

### bam-readcount

[bam-readcount](https://github.com/genome/bam-readcount) is a program for determing read support for individual variants (SNVs and Indels only). 

### BCFtools

[BCFtools](https://samtools.github.io/bcftools/) is an open source program for variant calling and manipulating files in Variant Call Format (VCF) or Binary Variant Call Format (BCF). 

### BWA

[BWA](http://bio-bwa.sourceforge.net) is a popular DNA alignment tool used for mapping sequences to a reference genome. It is available under an open source GPLv3 license.

### CNVkit

[CNVkit](https://cnvkit.readthedocs.io/en/stable/) is a python based copy number caller designed for use with hybrid capture. 

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

### Manta

[Manta](https://github.com/Illumina/manta) is a structural variant caller developed by Illumina and available on gitub under the GPL_v3 license. 

### mosdepth

[mosdepth](https://github.com/brentp/mosdepth) is a program for determining depth in sequencing data. 

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

### Strelka

[Strelka](https://github.com/Illumina/strelka) is a germline and somatic variant caller developed by illumina and available under an open source GPLv3 license. 

### StringTie

[StringTie](https://ccb.jhu.edu/software/stringtie/) is a software program to perform transcript assembly and quantification of RNAseq data. 

### Varscan

[Varscan](http://dkoboldt.github.io/varscan/) is a java program designed to call variants in sequencing data. It was developed at the Genome Institute at Washington University and is hosted 
on github. 

### vcf-annotation-tools

[VCF Annotation Tools](https://github.com/griffithlab/VAtools) is a python package that includes several tools to annotate VCF files with data from other tools. We will be using this for the purpose of adding bam readcounts to the vcf files.

### VEP

[VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) is a variant annotation tool developed by ensembl and written in perl. By default VEP will perform annotations by making web-based API queries however it is much faster to have a local copy of cache and fasta files. The AWS AMI image we’re using already has these files for hg38 in the directory ~/workspace/vep_cache as they can take a bit of time to download. 

### vt

[vt](https://github.com/atks/vt) is a variant tool set that discovers short variants from Next Generation Sequencing data. We will use this for the purpose of splitting multi-allelic variants.





# FIX BELOW LATER IN CASE IT IS NEEDED - DEPENDS ON WHAT MARKUS M DECIDES


## Install copyCat

copyCat is an R library for detecting copy number aberrations in sequencing data. The library is only available on github so we will have to use the BiocManager library to install a few of the underlying package dependencies. 

- Open R. 
    - In R run othe following lines of code:
        - devtools::install_github("chrisamiller/copycat")

## Install CNVnator

CNVnator is a depth based copy number caller. It is open source and available on github under a creative common public license (CCPL). To install we first download and extract the source code. CNVnator relies on a specific version of samtools which is distributed with CNVnator, so our first step is to run make on that samtools. To finish the installation process we can then run make in CNVnator’s main source directory.

```bash
# download and decompress
cd ~/workspace/bin

#download and install dependency package "root" from Cern (https://root.cern/install/). 

curl -OL https://root.cern/download/root_v6.20.08.Linux-ubuntu20-x86_64-gcc9.3.tar.gz
tar -xvzf root_v6.20.08.Linux-ubuntu20-x86_64-gcc9.3.tar.gz
source root/bin/thisroot.sh

wget https://github.com/abyzovlab/CNVnator/releases/download/v0.3.3/CNVnator_v0.3.3.zip
unzip CNVnator_v0.3.3.zip

# make the samtools dependency distributed with CNVnator
cd CNVnator_v0.3.3/src/samtools
make

# make CNVnator
cd ../
make

# make a symlink
ln -s ~/workspace/bin/CNVnator_v0.3.3/src/cnvnator ~/workspace/bin/cnvnator

# test installation
~/workspace/bin/cnvnator
```

MultiQC is a tool to create a single report with interactive plots for multiple bioinformatics analyses across many samples.
