---
title: BIOINFORMATIC PROCESSING OF RNA DATA

---

# Bioinformatic processing of RNA data

In this exercise set we will learn how to process data produced from an RNA-sequencing experiment. This will involve FASTQ file processing and quality control, alignment, gene expression abundance estimation, and detection of fusions from the data.

## The data
We will analyze RNA-sequencing data from the same samples we analyzed during the DNA-sequencing exercises. As a reminder, these were a sample from the breast cancer cell line HCC1395 (the tumor sample) and its matched lymphoblastoid cell line HCC1395BL (the normal sample). These are meant to represent a tumor and a normal tissue samples from a hypothetical patient. **As before, we are focusing our analyses on chromosomes 6 and 17 to reduce computation times.**

You will notice that our tumor and normal samples have different tissues of origin. While we will compare their RNA-seq samples to each other in this exercise set for practice, note that doing so would not make sense in a real analysis to e.g., uncover features of tumor vs normal.

These are the names of the samples we will analyze (**note that we have data from two sequencing lanes**- this was required to achieve the target total depth in the RNA-seq experiment):
`RNAseq_Tumor_Lane1`
`RNAseq_Tumor_Lane2`
`RNAseq_Norm_Lane1`
`RNAseq_Norm_Lane2`

The files are several gigabytes in size and can take a while to process, so we have run some of the preliminary analyses (FastQC/MultiQC and alignment) for you to make the exercises more manageable. Let's take a look at these analyses below.
```bash!
# FastQC and MultiQC have been run using commands like what you saw in Exercise Set 1.
# Copy the MultiQC HTML file to your local computer for inspection and then answer Question 1 below.
scp YOUR_KI_ID@c8cancergen02.ki.se:/nfs/course/inputs/data/rna/multiqc_report.html .
```
**Question 1:** Based on the MultiQC report, is there adapter contamination in the samples? Are there issues we should address with sequence duplication levels? (Hint: recall we are working with data only from chr6 and chr17. Also, see a brief discussion of sequence duplication in RNA-seq under `Common reasons for warnings` [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html)).

We have used the aligner [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) to perform spliced alignments to the reference genome. You can see the commands used below (**you do not need to run them yourself!**). For efficiency, the output of HISAT2 (SAM format) was piped directly to another program called [sambamba](http://lomereiter.github.io/sambamba/) to first convert to BAM format and then sort the BAM file. Before each command we also created and assigned a path for temporary directories.

```bash!
# All the commands below have already been run, they are only for your reference! 

# Align tumor data using HISAT2
# First for Lane 1
TUMOR_DATA_1_TEMP=`mktemp -d /nfs/course/inputs/data/rna/alignments/2895626107.XXXXXXXXXXXX`
hisat2 -p 7 --dta -x /nfs/course/inputs/references/transcriptome/ref_genome --rg-id 2895626107 --rg PL:ILLUMINA --rg PU:H3MYFBBXX-GCCAAT.4 --rg LB:rna_tumor_lib1 --rg SM:HCC1395_RNA --rna-strandness RF -1 /nfs/course/inputs/data/rna/RNAseq_Tumor/RNAseq_Tumor_Lane1_R1.fastq.gz -2  /nfs/course/inputs/data/rna/RNAseq_Tumor/RNAseq_Tumor_Lane1_R2.fastq.gz | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t 8 -m 24G --tmpdir $TUMOR_DATA_1_TEMP -o /nfs/course/inputs/data/rna/alignments/RNAseq_Tumor_Lane1.bam /dev/stdin

# Then for Lane 2
TUMOR_DATA_2_TEMP=`mktemp -d /nfs/course/inputs/data/rna/alignments/2895626112.XXXXXXXXXXXX`
hisat2 -p 7 --dta -x /nfs/course/inputs/references/transcriptome/ref_genome --rg-id 2895626112 --rg PL:ILLUMINA --rg PU:H3MYFBBXX-GCCAAT.5 --rg LB:rna_tumor_lib1 --rg SM:HCC1395_RNA --rna-strandness RF -1 /nfs/course/inputs/data/rna/RNAseq_Tumor/RNAseq_Tumor_Lane2_R1.fastq.gz -2 /nfs/course/inputs/data/rna/RNAseq_Tumor/RNAseq_Tumor_Lane2_R2.fastq.gz | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t 8 -m 24G --tmpdir $TUMOR_DATA_2_TEMP -o /nfs/course/inputs/data/rna/alignments/RNAseq_Tumor_Lane2.bam /dev/stdin

# Align normal data using HISAT2
# First for Lane 1
NORMAL_DATA_1_TEMP=`mktemp -d /nfs/course/inputs/data/rna/alignments/2895625992.XXXXXXXXXXXX`
hisat2 -p 7 --dta -x /nfs/course/inputs/references/transcriptome/ref_genome --rg-id 2895625992 --rg PL:ILLUMINA --rg PU:H3MYFBBXX-CTTGTA.4 --rg LB:rna_norm_lib1 --rg SM:HCC1395BL_RNA --rna-strandness RF -1 /nfs/course/inputs/data/rna/RNAseq_Norm/RNAseq_Norm_Lane1_R1.fastq.gz -2 /nfs/course/inputs/data/rna/RNAseq_Norm/RNAseq_Norm_Lane1_R2.fastq.gz | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t 8 -m 24G --tmpdir $NORMAL_DATA_1_TEMP -o /nfs/course/inputs/data/rna/alignments/RNAseq_Norm_Lane1.bam /dev/stdin

# Then for Lane 2
NORMAL_DATA_2_TEMP=`mktemp -d /nfs/course/inputs/data/rna/alignments/2895626097.XXXXXXXXXXXX`
hisat2 -p 7 --dta -x /nfs/course/inputs/references/transcriptome/ref_genome --rg-id 2895626097 --rg PL:ILLUMINA --rg PU:H3MYFBBXX-CTTGTA.5 --rg LB:rna_norm_lib1 --rg SM:HCC1395BL_RNA --rna-strandness RF -1 /nfs/course/inputs/data/rna/RNAseq_Norm/RNAseq_Norm_Lane2_R1.fastq.gz -2 /nfs/course/inputs/data/rna/RNAseq_Norm/RNAseq_Norm_Lane2_R2.fastq.gz | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t 8 -m 24G --tmpdir $NORMAL_DATA_2_TEMP -o /nfs/course/inputs/data/rna/alignments/RNAseq_Norm_Lane2.bam /dev/stdin

# Since we have multiple BAMs of each sample that just represent additional data for the same sequence library, we combine them into a single BAM for convenience:
sambamba merge -t 8 /nfs/course/inputs/data/rna/alignments/RNAseq_Norm.bam /nfs/course/inputs/data/rna/alignments/RNAseq_Norm_Lane1.bam /nfs/course/inputs/data/rna/alignments/RNAseq_Norm_Lane2.bam

sambamba merge -t 8 /nfs/course/inputs/data/rna/alignments/RNAseq_Tumor.bam /nfs/course/inputs/data/rna/alignments/RNAseq_Tumor_Lane1.bam /nfs/course/inputs/data/rna/alignments/RNAseq_Tumor_Lane2.bam

# Index the BAM files
samtools index RNAseq_Norm.bam
samtools index RNAseq_Tumor.bam
```

## GTF (General Transfer Format) and RNA-seq post-alignment quality control
We will encounter the GTF file format during the exercises below, described in detail [here](https://www.ensembl.org/info/website/upload/gff.html#fields). As you might recall from the lectures, the GTF format is used to describe genes and other features of DNA, RNA, and proteins.

In the section below we will inspect the files needed for RNA-seq post-alignment quality control.

```bash
#Go to this directory:
cd /nfs/course/inputs/references/transcriptome

#Have a peek into the transcriptome GTF file and try to make sense of it using the description behind the link above
less -SN ref_transcriptome.gtf

#Check which chromosomes we have transcripts from
cut -f1 ref_transcriptome.gtf | sort | uniq -c

#We will also use a BED file with ribosomal RNA coordinates in our analyses
#Check out the format of the ribosome BED file
#Answer Question 2 below
less -SN ref_ribosome.bed
```
**Question 2:** How could we have generated the ribosome file from our `ref_transcriptome.gtf` file? In other words, what information/fields/labels in the GTF file can be used to identify ribosomal RNAs?

Next, we will perform post-alignment quality control. For this, we will be using `samtools flagstat`, `fastqc`, and `picard CollectRnaSeqMetrics`.

```bash
#Start by making a folder for your RNA analyses
cd /nfs/course/students/YOUR_KI_ID
mkdir rna rna/qc
cd rna/qc

#Run samtools flagstat
#Runtime: ~2min, start and send to background using nohup
nohup samtools flagstat /nfs/course/inputs/data/rna/alignments/RNAseq_Norm.bam > /nfs/course/students/YOUR_KI_ID/rna/qc/RNAseq_Norm_flagstat.txt 2> flagstat_RNA_N_bam.out &
nohup samtools flagstat /nfs/course/inputs/data/rna/alignments/RNAseq_Tumor.bam > /nfs/course/students/YOUR_KI_ID/rna/qc/RNAseq_Tumor_flagstat.txt 2> flagstat_RNA_T_bam.out &

#Run FastQC
#Runtime: ~12 min
nohup fastqc -t 4 -o /nfs/course/students/YOUR_KI_ID/rna/qc /nfs/course/inputs/data/rna/alignments/RNAseq_Tumor.bam > fastqc_RNA_T_bam.out 2>&1 &
nohup fastqc -t 4 -o /nfs/course/students/YOUR_KI_ID/rna/qc /nfs/course/inputs/data/rna/alignments/RNAseq_Norm.bam > fastqc_RNA_N_bam.out 2>&1 &

#You can follow the progress of the program by looking in the nohup output files
less -SN fastqc_RNA_T_bam.out
less -SN fastqc_RNA_N_bam.out

#Run Picard CollectRnaSeqMetrics
#Notice that we run this for versions of the BAM files that contain only properly paired reads, as required by the tool
nohup picard CollectRnaSeqMetrics I=/nfs/course/inputs/data/rna/alignments/RNAseq_Norm_paired.bam O=/nfs/course/students/YOUR_KI_ID/rna/qc/RNAseq_Norm.RNA_Metrics REF_FLAT=/nfs/course/inputs/references/transcriptome/ref_flat.txt STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=/nfs/course/inputs/references/transcriptome/ref_ribosome.interval_list > collectRnaSeqMetrics_N.out 2>&1 &
nohup picard CollectRnaSeqMetrics I=/nfs/course/inputs/data/rna/alignments/RNAseq_Tumor_paired.bam O=/nfs/course/students/YOUR_KI_ID/rna/qc/RNAseq_Tumor.RNA_Metrics REF_FLAT=/nfs/course/inputs/references/transcriptome/ref_flat.txt STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=/nfs/course/inputs/references/transcriptome/ref_ribosome.interval_list > qcollectRnaSeqMetrics_T.out 2>&1 &

cd /nfs/course/students/YOUR_KI_ID/rna/qc/
mkdir post_align_qc
cd post_align_qc
multiqc /nfs/course/students/YOUR_KI_ID/rna/qc/

#Finally, download the MultiQC HTML report files to your local computer and view them in your favorite browser
#Familiarize yourself with the RNA-seq MultiQC data
scp YOUR_KI_ID@c8cancergen02.ki.se:/nfs/course/students/YOUR_KI_ID/rna/qc/post_align_qc/multiqc_report.html .
```
**Question 3:** In the MultiQC report, look at the metrics and plots from Picard and answer the following:
* For good quality RNA-seq data and its mapping, we would expect a large portion (some say at least 60%) of the reads to map to coding or UTR regions of the genome. Look at the RnaSeqMetrics Assignment plot. Is this the case for our samples?
* In the Gene Coverage plot, do you see evidence for 5'-3' bias?

You will have noticed we used the `nohup` command above. `nohup` is very useful for sending processes to the background so that they will continue if connection breaks or if you need to log out of the server.

For a review on how to send stdout and stderr to a file while running `nohup`, have a look [here](https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file).


## Run a simplified "reference only" StringTie expression abundance estimation approach

[`StringTie`](https://ccb.jhu.edu/software/stringtie/) is a fast and highly efficient assembler of RNA-seq alignments into potential transcripts. 

The `StringTie` developers recommend to: 
  - perform reference guided transcript compilation (aka transcript assembly) on each individual sample.
  - merge transcript predictions from all samples into a single model of the transcriptome.
  - annotate this predicted transcriptome with known transcriptome information.
  - estimate abundance for each of the transcripts in this final transcriptome model in each sample.
    
The final result csn have abundances for all transcripts across all samples. This includes a combination of predicted and known transcripts. 

It is sometimes convenient to have a more simplified workflow where we only have values for known transcripts. This is particularly true in species where we already have comprehensive high quality transcriptome annotations and there is less of a focus on de novo transcript discovery.

The following workflow produces such a "reference-only" transcriptome result in which we will perform abundance calculations on each lane of data individually.

```bash
cd /nfs/course/students/YOUR_KI_ID/rna/
mkdir ref-only-expression
cd ref-only-expression

nohup stringtie -p 4 -e -G /nfs/course/inputs/references/transcriptome/ref_transcriptome.gtf -o /nfs/course/students/YOUR_KI_ID/rna/ref-only-expression/RNAseq_Tumor_Lane1/transcripts.gtf -A /nfs/course/students/YOUR_KI_ID/rna/ref-only-expression/RNAseq_Tumor_Lane1/gene_abundances.tsv /nfs/course/inputs/data/rna/alignments/RNAseq_Tumor_Lane1.bam > knowntranscripts_T1_L1.out 2>&1 &
nohup stringtie -p 4 -e -G /nfs/course/inputs/references/transcriptome/ref_transcriptome.gtf -o /nfs/course/students/YOUR_KI_ID/rna/ref-only-expression/RNAseq_Tumor_Lane2/transcripts.gtf -A /nfs/course/students/YOUR_KI_ID/rna/ref-only-expression/RNAseq_Tumor_Lane2/gene_abundances.tsv /nfs/course/inputs/data/rna/alignments/RNAseq_Tumor_Lane2.bam > knowntranscripts_T1_L2.out 2>&1 &
nohup stringtie -p 4 -e -G /nfs/course/inputs/references/transcriptome/ref_transcriptome.gtf -o /nfs/course/students/YOUR_KI_ID/rna/ref-only-expression/RNAseq_Norm_Lane1/transcripts.gtf -A /nfs/course/students/YOUR_KI_ID/rna/ref-only-expression/RNAseq_Norm_Lane1/gene_abundances.tsv /nfs/course/inputs/data/rna/alignments/RNAseq_Norm_Lane1.bam > knowntranscripts_N1_L1.out 2>&1 &
nohup stringtie -p 4 -e -G /nfs/course/inputs/references/transcriptome/ref_transcriptome.gtf -o /nfs/course/students/YOUR_KI_ID/rna/ref-only-expression/RNAseq_Norm_Lane2/transcripts.gtf -A /nfs/course/students/YOUR_KI_ID/rna/ref-only-expression/RNAseq_Norm_Lane2/gene_abundances.tsv /nfs/course/inputs/data/rna/alignments/RNAseq_Norm_Lane2.bam > knowntranscripts_T2_L2.out 2>&1 &

#Wait until the scripts above have finished before proceeding.
```

Next we will create tidy expression matrix files from the `StringTie` results. This will be done at both the gene- and transcript-level and also will take into account the various expression measures produced: coverage, fragments per kilobase of transcript per million mapped reads (FPKM), and transcripts per million (TPM).

```bash
cd /nfs/course/students/YOUR_KI_ID/rna/ref-only-expression

#Download script for generating tidy matrix
wget https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/scripts/stringtie_expression_matrix.pl
chmod +x stringtie_expression_matrix.pl

#Then generate the expression matrices
#TPM
./stringtie_expression_matrix.pl --expression_metric=TPM --result_dirs='RNAseq_Norm_Lane1,RNAseq_Norm_Lane2,RNAseq_Tumor_Lane1,RNAseq_Tumor_Lane2' --transcript_matrix_file=transcript_tpm_all_samples.tsv --gene_matrix_file=gene_tpm_all_samples.tsv

#FPKM
./stringtie_expression_matrix.pl --expression_metric=FPKM --result_dirs='RNAseq_Norm_Lane1,RNAseq_Norm_Lane2,RNAseq_Tumor_Lane1,RNAseq_Tumor_Lane2' --transcript_matrix_file=transcript_fpkm_all_samples.tsv --gene_matrix_file=gene_fpkm_all_samples.tsv

#Coverage
./stringtie_expression_matrix.pl --expression_metric=Coverage --result_dirs='RNAseq_Norm_Lane1,RNAseq_Norm_Lane2,RNAseq_Tumor_Lane1,RNAseq_Tumor_Lane2' --transcript_matrix_file=transcript_coverage_all_samples.tsv --gene_matrix_file=gene_coverage_all_samples.tsv

#Have a look at the output files
head gene_coverage_all_samples.tsv transcript_coverage_all_samples.tsv gene_fpkm_all_samples.tsv transcript_fpkm_all_samples.tsv gene_tpm_all_samples.tsv transcript_tpm_all_samples.tsv
```

## Reference-free expression analysis with Kallisto

Remember that in previous sections we have been using reference genome FASTA sequences for the reference for alignment and subsequent steps. However, [`Kallisto`](https://pachterlab.github.io/kallisto/about) is a tool that works directly on target cDNA/transcript sequences. 

For `StringTie`, we have used transcript annotations for genes on our subset of chromosomes (i.e. chr6 and chr17). The transcript models were downloaded from Ensembl in GTF format. This GTF contains a description of the coordinates of exons that make up each transcript but it does not contain the transcript sequences themselves, which is something `Kallisto` needs. There are many places we could obtain such transcript sequences. For example, we could have download them directly in FASTA format from the Ensembl FTP site (or from UCSC or NCBI). Here, we have downloaded them for you.

```bash
cd /nfs/course/students/YOUR_KI_ID/rna/
mkdir kallisto
cd kallisto

#First check that the GTF and FASTA files are present for Kallisto
head /nfs/course/inputs/references/transcriptome/ref_transcriptome.gtf
head /nfs/course/inputs/references/kallisto/ref_transcriptome_clean.fa

#Now check for the kallisto index is there
ls -halt /nfs/course/inputs/references/kallisto/ref_transcriptome_kallisto_index

#Create a list of all transcript IDs for later use:
cat /nfs/course/inputs/references/kallisto/ref_transcriptome_clean.fa | grep ">" | perl -ne '$_ =~ s/\>//; print $_' | sort | uniq > transcript_id_list.txt
head -n 10 transcript_id_list.txt
```

## Generate abundance estimates for all samples using Kallisto
As we did with `StringTie`, we will generate transcript abundances for each of our demonstration samples using `Kallisto`. Here we are treating the two lanes for each sample as if they were independent samples.

```bash
cd /nfs/course/students/YOUR_KI_ID/rna/kallisto/
mkdir quants
cd quants

nohup /nfs/course/bin/kallisto/kallisto quant --index=/nfs/course/inputs/references/kallisto/ref_transcriptome_kallisto_index --output-dir=RNAseq_Norm_Lane1 --threads=2 --plaintext /nfs/course/inputs/data/rna/RNAseq_Norm/RNAseq_Norm_Lane1_R1.fastq.gz /nfs/course/inputs/data/rna/RNAseq_Norm/RNAseq_Norm_Lane1_R2.fastq.gz > kallisto_quant_N_L1.out 2>&1 &
nohup /nfs/course/bin/kallisto/kallisto quant --index=/nfs/course/inputs/references/kallisto/ref_transcriptome_kallisto_index --output-dir=RNAseq_Norm_Lane2 --threads=2 --plaintext /nfs/course/inputs/data/rna/RNAseq_Norm/RNAseq_Norm_Lane2_R1.fastq.gz /nfs/course/inputs/data/rna/RNAseq_Norm/RNAseq_Norm_Lane2_R2.fastq.gz > kallisto_quant_N_L2.out 2>&1 &
nohup /nfs/course/bin/kallisto/kallisto quant --index=/nfs/course/inputs/references/kallisto/ref_transcriptome_kallisto_index --output-dir=RNAseq_Tumor_Lane1 --threads=2 --plaintext /nfs/course/inputs/data/rna/RNAseq_Tumor/RNAseq_Tumor_Lane1_R1.fastq.gz /nfs/course/inputs/data/rna/RNAseq_Tumor/RNAseq_Tumor_Lane1_R2.fastq.gz > kallisto_quant_T_L1.out 2>&1 &
nohup /nfs/course/bin/kallisto/kallisto quant --index=/nfs/course/inputs/references/kallisto/ref_transcriptome_kallisto_index --output-dir=RNAseq_Tumor_Lane2 --threads=2 --plaintext /nfs/course/inputs/data/rna/RNAseq_Tumor/RNAseq_Tumor_Lane2_R1.fastq.gz /nfs/course/inputs/data/rna/RNAseq_Tumor/RNAseq_Tumor_Lane2_R2.fastq.gz > kallisto_quant_T_L2.out 2>&1 &

#Check nohup output
less -SN kallisto_quant_T_L2.out

#Create a single TSV file that has the TPM abundance estimates for all samples.
#First check the contets of an abudance output file from Kallisto
less -SN ./RNAseq_Norm_Lane1/abundance.tsv

#Merge all files
paste */abundance.tsv | cut -f 1,2,5,10,15,20 > transcript_tpms_all_samples.tsv
ls -1 */abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\tlength$_\n"' > header.tsv
cat header.tsv transcript_tpms_all_samples.tsv | grep -v "tpm" > transcript_tpms_all_samples.tsv2
mv transcript_tpms_all_samples.tsv2 transcript_tpms_all_samples.tsv
rm -f header.tsv
head transcript_tpms_all_samples.tsv

#To create a gene version of the Kallisto TPM matrix, we can sum the TPM values for transcripts of the same gene.

#Remember, one gene can have many transcripts. For example, the Ensembl gene ID for TP53 is ENSG00000141510. Check which Ensembl transcripts exist in the reference GTF for TP53 and how many there are.
#We will use this information later in one of our final visualizations
grep ENST /nfs/course/inputs/references/transcriptome/ref_transcriptome.gtf | grep ENSG00000141510 | perl -nle '@a = split /\t/; $_ =~ /transcript_id\s"(\S*)";/g; print $1;' | sort | uniq

#Download a script for generating the gene matrix
wget https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/scripts/kallisto_gene_matrix.pl
chmod +x kallisto_gene_matrix.pl

#Generate a gene-level matrix
./kallisto_gene_matrix.pl --gtf_file=/nfs/course/inputs/references/transcriptome/ref_transcriptome.gtf --kallisto_transcript_matrix_in=transcript_tpms_all_samples.tsv --kallisto_transcript_matrix_out=gene_tpms_all_samples.tsv
less -SN gene_tpms_all_samples.tsv
```

# Compare expression values between Kallisto and StringTie

To compare the `StringTie` and `Kallisto` approaches, we can use the expression value for each Ensembl transcript.

To do this comparison, we need to gather the expression estimates for each of our replicates from each approach. We will perform this comparison in R/RStudio on our local machines.

First download the needed files to your local machine:

```bash
#Run these commands on your local machine

#Download the Kallisto expression estimates
scp -r YOUR_KI_ID@c8cancergen02.ki.se:/nfs/course/students/YOUR_KI_ID/rna/kallisto/quants/transcript_tpms_all_samples.tsv .
mv transcript_tpms_all_samples.tsv kallisto_transcript_tpms_all_samples.tsv

#Download the StringTie expression estimates
scp -r YOUR_KI_ID@c8cancergen02.ki.se:/nfs/course/students/YOUR_KI_ID/rna/ref-only-expression/transcript_tpm_all_samples.tsv .
mv transcript_tpm_all_samples.tsv stringtie_transcript_tpms_all_samples.tsv

#Start R or RStudio on your local machine

#Set the working directory to where the Kallisto and StringTie files were downloaded above
setwd("DIR GOES HERE")

# load libraries
library(ggplot2)
library(reshape2)

# read in data
kallisto_transcript_tpm <- read.delim("kallisto_transcript_tpms_all_samples.tsv")
head(kallisto_transcript_tpm)

stringtie_transcript_tpm <- read.delim("stringtie_transcript_tpms_all_samples.tsv")
head(stringtie_transcript_tpm)

# minor reformatting
kallisto_transcript_tpm <- kallisto_transcript_tpm[,-2]
kallisto_transcript_tpm <- melt(kallisto_transcript_tpm, id.vars=c("target_id"))
head(kallisto_transcript_tpm)

stringtie_transcript_tpm <- melt(stringtie_transcript_tpm, id.vars=c("Transcript_ID"))
head(stringtie_transcript_tpm)

# merge the data
kallisto_stringtie_tpm <- merge(kallisto_transcript_tpm, stringtie_transcript_tpm, by.x=c("target_id", "variable"), by.y=c("Transcript_ID", "variable"), suffixes=c(".kallisto", ".stringtie"))

#The numeric vector is read as a character vector, correcting ..
kallisto_stringtie_tpm$value.stringtie = as.numeric(kallisto_stringtie_tpm$value.stringtie)
str(kallisto_stringtie_tpm) #Now should be numeric

#TP53 transcript ID vector from the transcript GTF (we determined these earlier)
TP53 = c("ENST00000269305", "ENST00000359597", "ENST00000413465", "ENST00000420246", "ENST00000445888", "ENST00000455263", "ENST00000503591", "ENST00000504290", "ENST00000504937", "ENST00000505014", "ENST00000508793", "ENST00000509690", "ENST00000510385", "ENST00000514944", "ENST00000574684", "ENST00000576024", "ENST00000604348", "ENST00000610292", "ENST00000610538", "ENST00000610623", "ENST00000615910", "ENST00000617185", "ENST00000618944", "ENST00000619186", "ENST00000619485", "ENST00000620739", "ENST00000622645", "ENST00000635293")

idx = match(kallisto_stringtie_tpm$target_id, TP53)
idx.1 = which(is.na(idx) == FALSE)
idx.2 = idx[idx.1]
TP53[idx.2] == kallisto_stringtie_tpm$target_id[idx.1]
#Check TP53 expression values
kallisto_stringtie_tpm[idx.2,]
#To plot TP53, make a new column
kallisto_stringtie_tpm$TP53 = "Other transcript"
kallisto_stringtie_tpm$TP53[idx.2] = "TP53"

kallisto_stringtie_tpm$TP53 = factor(kallisto_stringtie_tpm$TP53, levels = c("Other transcript", "TP53"))
str(kallisto_stringtie_tpm)

#Remove 0/NA values from the expression estimates
#This is needed for visualization on the log2 scale below
kallisto_stringtie_tpm$value.stringtie[which(is.na(kallisto_stringtie_tpm$value.stringtie)==TRUE)] = 0.001
kallisto_stringtie_tpm$value.stringtie = kallisto_stringtie_tpm$value.stringtie + 0.1
kallisto_stringtie_tpm$value.kallisto = kallisto_stringtie_tpm$value.kallisto + 0.1

########### Plot the results ######################

#Plot a dot plot of expression estimates using the two methods
#Expression values are plotted on the log2 scale
ggplot(data=kallisto_stringtie_tpm, aes(x=log2(value.kallisto), y=log2(value.stringtie), colour=TP53)) +
geom_point(alpha=1) +
scale_colour_manual(values = c("forestgreen", "firebrick", "dodgerblue"), name = "Transcripts", drop=TRUE) +
geom_point(data = subset(kallisto_stringtie_tpm, TP53 == "TP53"), aes(x=log2(value.kallisto), y=log2(value.stringtie), colour=TP53), inherit.aes=FALSE)+
facet_wrap(~variable)+
theme_bw()

#Plot density disbributions of the expression estimates
df = rbind(data.frame(variable = kallisto_stringtie_tpm$variable, value =  kallisto_stringtie_tpm$value.kallisto, type = "kallisto"), data.frame(variable = kallisto_stringtie_tpm$variable, value =  kallisto_stringtie_tpm$value.stringtie, type = "stringtie"))
df$type = factor(df$type, levels =c("kallisto","stringtie"))
head(df)

ggplot(data=df, aes(x=log2(value))) +
geom_density() +
facet_wrap(type~variable, ncol = 4)+
ggtitle("Density distributions")+
xlab("Log2 expression value")+
theme(legend.position = "none")
#dev.off()
```
**Question 4:** In the plots of expression values from Kallisto vs StringTie, would you say that the expression values are similar between the two tools? If you observe differences, why do you think there are differences? 
**Question 5:** If you were interested in quantifying the expression levels of **known transcripts**, would you choose StringTie or Kallisto? Why?

# Gene fusion detection

## Introduction

In addition to providing information about gene expression, RNA-seq data can be used to discover transcripts which result from chromosomal translocations. Translocations and their resultant chimeric (AKA fusion) transcripts are important driver mutations in many cancers. A [variety of specialized alignment and filtering strategies](https://www.ncbi.nlm.nih.gov/pubmed/27485475) have been developed to identify fusion transcripts from RNA, but these programs can suffer from low specificity (i.e. many false-positives) and poor correlation across methods.

This exercise uses the [kallisto](https://pachterlab.github.io/kallisto/about) and [pizzly](https://github.com/pmelsted/pizzly) tools for fusion detection from RNA-seq data. Kallisto quantifies transcript abundance through [pseudoalignment](http://tinyheero.github.io/2015/09/02/pseudoalignments-kallisto.html). Pizzly aligns reads which kallisto has flagged as potentially spanning fusion junctions.

Some of the other tools used to detection gene fusions from RNA-seq are Arriba, EricScript, FusionCatcher, Fusion-Inspector, fusion-report, Squid, and Star-Fusion. All these tools have their own merits and demerits. They are often used in combinations or all together. Hence, in a real analysis, we need to assess which one is suitable for our particular data.

## Setup information

Pizzly does not work with the most recent Ensembl human GTF annotation file. We have therefore downloaded the version 87 GTF and subset the reference transcriptome FASTA file accordingly. These new reference data are located under `/nfs/course/inputs/references/rna_fusion`.

We have also merged the FASTQ files from separate sequencing lanes into single files to generate one R1 and one R2 file for the Tumor and Normal samples. We the further subsampled these (using `seqtk sample`) to only contain 30% of the reads for each sample, to allow for fusion alignment to complete in a reasonable time. The files are under `/nfs/course/inputs/data/rna/fusion`.
 
## Calling fusions with Kallisto and pizzly

```bash
#Create a folder for fusion results
cd /nfs/course/students/YOUR_KI_ID/rna
mkdir fusion
cd fusion

#Quantify potential fusions using Kallisto:
nohup /nfs/course/bin/kallisto/kallisto quant -i /nfs/course/inputs/references/rna_fusion/index.617.idx --fusion -o kquant-norm617 /nfs/course/inputs/data/rna/fusion/subRNAseq_Norm_R1.fastq.gz /nfs/course/inputs/data/rna/fusion/subRNAseq_Norm_R2.fastq.gz 2> kallisto_quant_N_bam.out &
nohup /nfs/course/bin/kallisto/kallisto quant -i /nfs/course/inputs/references/rna_fusion/index.617.idx --fusion -o kquant-tumor617 /nfs/course/inputs/data/rna/fusion/subRNAseq_Tumor_R1.fastq.gz /nfs/course/inputs/data/rna/fusion/subRNAseq_Tumor_R2.fastq.gz 2> kallisto_quant_T_bam.out &

#Call fusions with pizzly
#Note that this step can take 
nohup pizzly -k 31 --gtf /nfs/course/inputs/references/rna_fusion/chr617.gtf --cache index-norm617.cache.txt --align-score 2 --insert-size 400 --fasta /nfs/course/inputs/references/rna_fusion/chr617.fa --output norm-fuse617 kquant-norm617/fusion.txt 2> pizzly_N_bam.out &
nohup pizzly -k 31 --gtf /nfs/course/inputs/references/rna_fusion/chr617.gtf --cache index-tumor617.cache.txt --align-score 2 --insert-size 400 --fasta /nfs/course/inputs/references/rna_fusion/chr617.fa --output tumor-fuse617 kquant-tumor617/fusion.txt 2> pizzly_T_bam.out &
```
**Question 6:** Based on the `pizzly_N_bam.out` and `pizzly_T_bam.out` files, how many records (transcripts) were kept by pizzly for each sample from those that were suggested by kallisto? 

## Investigating the output of pizzly fusion calling
Investigation of the pizzly output will be done on your local computer using R/RStudio. See below to first download the required files to your local machine, followed by code in R to visualize the fusion calling output.

```bash
#To go a directory on your local machine where you would like to download the fusion calling files
scp YOUR_KI_ID@c8cancergen02.ki.se:/nfs/course/students/YOUR_KI_ID/rna/fusion/*fuse617* .

# Also download R scripts for later use
scp YOUR_KI_ID@c8cancergen02.ki.se:/nfs/course/inputs/data/rna/scripts/* .

# Open R or RStudio on your computer and run the following code
#Change working directory to location of fusion files
setwd("path_to_fusion_files_local_machine")
#Check the path
getwd()

# Load required packages
library(EnsDb.Hsapiens.v86)
library(chimeraviz)
library(cowplot)
library(jsonlite)
library(dplyr)
library(ggplot2)

# Load data
suffix = "617.json"
JSON_files = list.files(path = ".", pattern = paste0("*",suffix))
print(JSON_files)
Ids = gsub(suffix, "", JSON_files)
print(Ids)

# Loading gene models from the annotation package
edb <- EnsDb.Hsapiens.v86
listColumns(edb)
supportedFilters(edb)

# Load one of the downloaded R scripts for parsing of the pizzly output
# This contains the function GetFusionz_and_namez() that we use below
source("mod.grolar.R")

# This function will flatten out the JSON giving you a list of gene A and gene B sorted by splitcount then paircount
# Each fusion now has a unique identifier, sequence positions, and distance values for genes from the same chromosome
lapply(Ids, function(x) GetFusionz_and_namez(x, suffix = "617.json") )

# Read the flattened files produced by the function
normal=read.table("norm-fuse_fusions_filt_sorted.txt", header=T)
tumor=read.table("tumor-fuse_fusions_filt_sorted.txt", header=T)

#Investigate the first five rows
head(normal, 5)
head(tumor, 5)

# Common filtering tasks for fusion output include removing fusions from the tumor sample which are present in the normal, and removing fusions for which the is little support by pair and split read counts:
normal$genepairs=paste0(normal$geneA.name,".",normal$geneB.name)
tumor$genepairs=paste0(tumor$geneA.name,".",tumor$geneB.name)
uniqueTumor=subset(tumor, !(tumor$genepairs %in% normal$genepairs))
nrow(uniqueTumor)!=nrow(tumor) # Did we filter out any fusions from the tumor file?
nrow(tumor)-nrow(uniqueTumor) # How many fusions did we filter out?

# There are two fusions (or at least fusion gene pairs) from the normal sample which are also present in the tumor. 
# Examine these shared fusion gene pairs:
shared_src_tumor=subset(tumor, (tumor$genepairs %in% normal$genepairs))
shared_src_normal=subset(normal, (normal$genepairs %in% tumor$genepairs))
shared_src_tumor
shared_src_normal

# Merge normal and tumor data
normal$sample="normal"
tumor$sample="tumor"
allfusions=rbind(normal,tumor)

# Compare counts of paired and split reads
tapply(allfusions$paircount, allfusions$sample, summary)
tapply(allfusions$splitcount, allfusions$sample, summary)

# As a density plot
p1=ggplot(allfusions, aes(paircount, fill=sample))+
  geom_density(alpha=.4)+
  geom_vline(xintercept=2)+
  coord_fixed(ratio=15)

p2=ggplot(allfusions, aes(splitcount, fill=sample))+
  geom_density(alpha=.4)+coord_cartesian(ylim= c(0,.2))+
  geom_vline(xintercept=5)+
  coord_fixed(ratio=200)

plot_grid(p1,p2, ncol = 2, rel_heights = c(1,1))

# Filter by requiring each fusion to be supported by two or more paired fusion reads and five or more split fusion reads
allfusions_filtered=allfusions[which(allfusions$paircount >= 2 & allfusions$splitcount >= 5),]
nrow(allfusions)
nrow(allfusions_filtered) #How many fusions are left after filtering?

# Write out a table of the filtered fusions
write.table(allfusions_filtered, file="allfusions.txt",row.names=F, col.names=T, quote=F, sep="\t")

#Next, we will do a visualization using chimeraviz
#Chimeraviz is an R package for visualizing fusions from RNA-seq data. The chimeraviz package has import functions built in for a variety of fusion-finder programs, but not for pizzly. We will have to load our own import function that you downloaded above:
source("import_Pizzly.R")

#The code above loads the importPizzly() function
#We use it below to generate a circos plot

fusions = importPizzly("allfusions.txt","hg38")
plot_circle(fusions)
```
**Question 7:** Above we performed fusion filtering by requiring each fusion to be supported by two or more paired fusion reads and five or more split fusion reads. How many fusions were left after this filtering? How many fusions were present prior to filtering?

**Question 8:** In the filtered fusion table`allfusions.txt`, what are the two genes involved in the fusion with the highest paired and split fusion read support?
