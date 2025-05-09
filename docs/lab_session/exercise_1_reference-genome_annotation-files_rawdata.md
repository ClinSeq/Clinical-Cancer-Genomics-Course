---
title: Clinical Cancer Genomics Course Exercise Set 1 Part 2

---

# Parts of a bioinformatics analysis: reference genome, annotation files, and raw data files

In this part of Exercise Set 1, we will be briefly exploring some input file types that we will later use to analyze DNA-sequencing data in Exercise Set 2 of the course. These input files include the human reference genome, annotation files for variation in the human genome, and raw sequencing data files.

## Human Reference Genome
As we learned during the lecture, the human reference genome is a work-in-progress. Different versions (referred to as *assemblies*) exist, and these assemblies have been updated in smaller ways with new information in the form of patches. If you are interested, you can find a great overview of the different human reference genome versions [here](https://genome.ucsc.edu/FAQ/FAQreleases.html).

Importantly, in genomics, the reference genome offers a scaffold upon which new data can be mapped, which is a much more efficient way rather than building a genome from scratch.

During this course, we will mostly use the hg38 human reference genome assembly. We will be using a smaller version of the reference genome containing only chromosomes 6 and 17 to shorten our analysis runtimes. The files are already available in our computing environment.

### Exploring the human reference genome file
Log into our course computing environment. Then, run the commands below and answer the accompanying questions to explore the human reference genome file.

```bash
# Go to your folder in the computing environment
# Replace YOUR_KI_ID below with your specific ID
cd /nfs/course/students/YOUR_KI_ID

# Create a directory for the reference genome files
mkdir references references/genome/
cd references/genome/

# Download the reference genome and unpack it
wget http://genomedata.org/pmbio-workshop/references/genome/chr6_and_chr17/ref_genome.tar
tar -xvf ref_genome.tar
gunzip ref_genome.fa.gz
rm -f ref_genome.fa.fai ref_genome.dict
rm ref_genome.tar

# View the downloaded files
ls -halt

# Check what the reference genome looks like. When you have the reference genome open just type "1000" and press enter. You will then jump 1000 lines.
less -SN ref_genome.fa

# The reference genome is in the FASTA format.
# Check the chromosome headers in the FASTA file.
cat ref_genome.fa | grep "^>"

# View the first 10 lines of this file. Note the header line starting with `>`. 
# Answer Question 1 below.
head -n 10 ref_genome.fa

# Determine the count of each base in the entire reference genome file (skipping the header lines for each sequence). 
# Note that this takes several seconds to run.
# Answer Question 2 below.
cat ref_genome.fa | grep -v ">" | perl -ne 'chomp $_; $bases{$_}++ for split //; if (eof){print "$_ $bases{$_}\n" for sort keys %bases}'

```
**Question 1:** Above you viewed the first 10 lines of the FASTA file and should have noticed that the chromosome 6 sequence started with many Ns. Why do you think this is? (*Hint: N is used to denote an unknown base, i.e. the base could be A, C, T, or G, but we are not sure which it is.*)

**Question 2:** What is the count of each base in the entire reference genome file (skipping the header lines for each sequence)? What does each of these bases refer to? What are the "unexpected bases"? (*Hint: Google "IUPAC nucleic acid codes"*)

### Splitting the human reference genome file by chromosome
Occasionally, some bioinformatics tools might expect our genome FASTA file to be split by chromosome. We can achieve this with the `faSplit` utility.
```bash
# Make a directory for the results
mkdir ref_genome_split

# Split by chromosome
faSplit byname ref_genome.fa ref_genome_split/

# View folder contents
tree
```
### Indexing the human reference genome file
Indexing is widely used in bioinformatics workflows to improve performance. Typically, it is applied to large data files with many records to improve the ability of tools to rapidly access specific locations of the file. For example, if we have alignments against the entire genome and we want to visualize the alignments for a single gene on chr17 in a genome viewer such as IGV, we don’t want the viewer to have to scan through the entire file. Indexing allows us to jump to the correct place in the file and pull out just the information we need without reading much of the file.

Below, we create some index files that we will need for future Exercise Sets.

```bash
# Use samtools to create a FASTA index file
samtools faidx ref_genome.fa

# View the contents of the index file
# Answer Question 3 below.
head ref_genome.fa.fai

# Use picard to create a dictionary file.
picard CreateSequenceDictionary -R ref_genome.fa -O ref_genome.dict

# View the content of the dictionary file.
head ref_genome.dict
# less can also be applied

#Also index the split chromosomes.
samtools faidx ./ref_genome_split/chr6.fa
samtools faidx ./ref_genome_split/chr17.fa
```
**Question 3:** Based on viewing the contents on the index file, what information does the index file store? (*Hint: see https://www.htslib.org/doc/faidx.html for an explanation of index file contents*)

## Annotation files
In addition to the reference genome, we usually need one or more **annotation files** to perform our bioinformatics analyses. Broadly speaking, these are files that provide information and coordinates for features in the genome, with features being for example genes, transcripts, exons, variants, or similar.

Below we will take a look at some of the annotation files we will be using for the analysis of DNA-sequencing data. There are more we will encounter in Exercise Set 2.

To start, one of the bioinformatics analyses that can require several annotation files is **variant calling** (our topic on Wednesday). We download these below for variant calling with the Genome Analysis Toolkit (GATK). There are annotation files of different types:
* **Known sites**: lists of variants that have been previously identified and reported, such as dbSNP below.
* **Training sets**: lists of variants that are used by machine learning algoriths to model the properties of true variation vs artifacts.
* **Truth sets**: lists of variants that are used to evaluate the quality of a variant callset (eg. sensitivity and specificity, or recall).

```bash
# Remember to replace YOUR_KI_ID below with your ID
cd /nfs/course/students/YOUR_KI_ID/references
mkdir -p gatk
cd gatk

# SNP calibration call sets - dbSNP, hapmap, omni, and 1000 Genomes
# Make symbolic links of these files
ln -s /nfs/course/inputs/references/gatk/Homo_sapiens_assembly38.dbsnp138.vcf.gz* .
ln -s /nfs/course/inputs/references/gatk/hapmap_3.3.hg38.vcf.gz* .
ln -s /nfs/course/inputs/references/gatk/1000G_omni2.5.hg38.vcf.gz* .
ln -s /nfs/course/inputs/references/gatk/1000G_phase1.snps.high_confidence.hg38.vcf.gz* .

# Indel calibration call sets - dbSNP, Mills
ln -s /nfs/course/inputs/references/gatk/Homo_sapiens_assembly38.known_indels.vcf.gz* .
ln -s /nfs/course/inputs/references/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz* .

# Interval lists that can be used to parallelize certain GATK tasks
ln -s /nfs/course/inputs/references/gatk/wgs_calling_regions.hg38.interval_list .
ln -s /nfs/course/inputs/references/gatk/scattered_calling_intervals .

# List the files we just made symbolic links of
ls -lh

# Look into the Mills and 1000 Genomes gold standard indels set and answer Question 4 below
gunzip -c Mills_and_1000G_gold_standard.indels.hg38.vcf.gz | grep -v "##" | head -n 2
```
**Question 4**: What information do the CHROM, POS, REF, and ALT fields contain in the Mills and 1000 Genomes gold standard indel set? For the indel we viewed above, is it an insertion or a deletion?

For many of our analyses, we will also need a file that specifies the coordinates of the regions we have sequenced from our samples. This allows us to focus our analyses on the regions we are interested in and have reliable information from. As we will see in the section below, the data we will analyze is generated via exome sequencing (meaning only the exonic regions of the genome have been sequenced). 

The reagent used to produce the exome data was the SeqCapEZ_Exome_v3.0 from Roche Nimblegen. In our computing environment we have two files, `SeqCap_EZ_Exome_v3_hg19_capture_targets.bed` and `SeqCap_EZ_Exome_v3_hg19_primary_targets.bed`, from Roche that contain the chromosome, start, stop, and gene annotation for each probe used in the reagent. The only problem is that the files use the hg19 genome assembly, whereas our analysis will use the newer hg38 assembly. In the section below we therefore convert the hg19 coordinates to hg38 coordinates using the [liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) tool.
```bash
# Make a new directory for exome files and copy them
# Remember to replace YOUR_KI_ID below with your ID
cd /nfs/course/students/YOUR_KI_ID/references
mkdir exome
cd exome
cp /nfs/course/inputs/references/gatk/SeqCap_EZ_Exome_v3_hg19* .

# Download a chain file that provides a mapping between the hg19 and hg38 assemblies for conversion
wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# Use liftOver to convert hg19 to hg38
# Answer Question 5 below
/nfs/course/inputs/references/exome/liftOver SeqCap_EZ_Exome_v3_hg19_primary_targets.bed  hg19ToHg38.over.chain.gz SeqCap_EZ_Exome_v3_hg38_primary_targets.bed unMapped.bed

/nfs/course/inputs/references/exome/liftOver SeqCap_EZ_Exome_v3_hg19_capture_targets.bed  hg19ToHg38.over.chain.gz SeqCap_EZ_Exome_v3_hg38_capture_targets.bed unMapped1.bed

# Create a version in standard BED format (chr, start, stop- removing gene information)
cut -f 1-3 SeqCap_EZ_Exome_v3_hg38_primary_targets.bed > SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed

cut -f 1-3 SeqCap_EZ_Exome_v3_hg38_capture_targets.bed > SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.bed

# Now create a subset of these BED files for the chromosomes we are using in our analysis (chr6 and chr17)
grep -w -P "^chr6|^chr17" SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed > exome_regions.bed

grep -w -P "^chr6|^chr17" SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.bed > probe_regions.bed
```
**Question 5:** Based on the `Unmapped` files created during the liftOver process, how many primary and capture targets could not be mapped from hg19 to hg38 because they were partially or completely deleted in the new assembly?

## Investigating raw DNA-sequencing data files
In this last section of the Exercise Set, we will familiarize ourselves with the DNA-sequencing data we will be further analyzing in Exercise Set 2 on Wednesday. We will also perform some quality control analyses to check the quality of our data.

### About the data
We will be analyzing two samples: `Exome_Tumor` and `Exome_Norm`. These samples are from a well-described breast cancer cell line (HCC1395, the tumor sample) and its matched lymphoblastoid cell line (HCC1395BL, the normal sample). These are meant to represent a tumor and a normal tissue samples from a hypothetical patient.

To generate the data we are about to analyze, exome sequencing was performed on the cell line samples. This is a form of targeted DNA-sequencing where only the exonic regions are sequenced.

### Overview of the raw data
The raw data is in the form of FASTQ files, which we learned about during the lecture. Let's inspect these files now in the computing environment.

```bash
# Make a new folder to hold the FASTQ files
# Remember to replace YOUR_KI_ID with you ID below
cd /nfs/course/students/YOUR_KI_ID/
mkdir data data/fastq
cd data/fastq

# Copy the FASTQ files into your folder
cp -r /nfs/course/inputs/data/fastq/* .

# List all files
tree
```

## Investigate the FASTQ files
```bash
cd /nfs/course/students/YOUR_KI_ID/data/fastq/

# Show the first ten lines of the Exome Tumor FASTQ files
zcat Exome_Tumor/Exome_Tumor_R1.fastq.gz | head
zcat Exome_Tumor/Exome_Tumor_R2.fastq.gz | head

# This wikipedia file gives a very nice overview of FASTQ files and Illumina base qualities
# https://en.wikipedia.org/wiki/FASTQ_format
# Let's look at the FASTQ fields for the first read and then answer Question 6 below
zcat Exome_Tumor/Exome_Tumor_R1.fastq.gz | head -n 4

# What is the length of each read?
zcat Exome_Tumor/Exome_Tumor_R1.fastq.gz | head -n 2 | tail -n 1 | wc
# The length of each read is 101 bases

# How many lines are there in the Exome_Tumor file
zcat Exome_Tumor/Exome_Tumor_R1.fastq.gz | wc -l 
# There are: 33,326,620

# How many paired reads or fragments are there?
expr 33326620 / 4 
# There are: 8,331,655 paired end reads

# How many total bases of data are in the Exome Tumor dataset?
echo "8331655 * (101 * 2)" | bc 
# There are: 1,682,994,310 bases of data

# How many total bases are there when expressed as "gigabases" (specify 2 decimal points using `scale`)
echo "scale=2; (8331655 * (101 * 2))/1000000000" | bc 
# There are: 1.68 Gbp of data

# What is the average coverage we expect to achieve with this much data for the targeted exome region?
# First determine the size of our exome regions (answer = 6683920). 
cat /nfs/course/students/YOUR_KI_ID/references/exome/exome_regions.bed | perl -ne 'chomp; @l=split("\t", $_); $size += $l[2]-$l[1]; if (eof){print "size = $size\n"}' 

# Now determine the average coverage of these positions by our bases of data
echo "scale=2; (8331655 * (101 * 2))/6683920" | bc 
# Average coverage expected = 251.79x

# Answer Question 8 below.
```
**Question 6:** Why are there two FASTQ files per sample? What do the R1 and R2 in the file names refer to?

**Question 7:** What is the symbol and numeric quality value (Q-score) of the lowest quality base(s) in the first read in Exome_Tumor_R1.fastq.gz? (*Hint: see for example https://help.basespace.illumina.com/files-used-by-basespace/quality-scores for an explanation of the quality symbols*)

**Question 8:** Above we calculated the average coverage of our data. What is the fundamental assumption of this calculation that is at least partially not true? What effect will this have on the observed coverage?

## Run FastQC to check the quality of the FASTQ-files
Next, we will run some quality checks on the FASTQ data. DNA-sequencing quality control will be discussed more in depth during Wednesday's lectures, but we will already get started today.

Have a look [here](https://www.youtube.com/watch?v=lUk5Ju3vCDM) for a short tutorial on the FastQC tool output.
```bash
cd /nfs/course/students/YOUR_KI_ID/data/fastq/

fastqc -t 4 Exome_Norm/Exome_Norm*.fastq.gz
fastqc -t 4 Exome_Tumor/Exome_Tumor*.fastq.gz
tree
```

## Run MultiQC 
We will run MultiQC to provide a combined report of the FastQC data. 

More info on MultiQC is available [here](https://multiqc.info).

```bash
cd /nfs/course/students/YOUR_KI_ID/
mkdir qc
cd qc
multiqc /nfs/course/students/YOUR_KI_ID/data/fastq/
tree

# Download MultiQC HTML output to the local computer, open the .html in you favourite browser.
scp -r YOUR_KI_ID@c8cancergen02.ki.se:/nfs/course/students/YOUR_KI_ID/qc/multiqc* .
```
**Question 9:** In your MultiQC report, select the Sequence Quality Histograms section from the menu on the left hand side. Here we can see the mean quality value across each base position in our sequencing reads. Does our data seem to be of good quality, or are there issues with low quality scores somewhere in the reads (and if so, where?)? 

*Hint: you can read more about Phred quality scores and what they mean [here](mean quality value across each base position in the read).*
