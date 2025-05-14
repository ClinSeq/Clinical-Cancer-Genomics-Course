---
title: Bioinformatic processing of DNA data

---

# Bioinformatic processing of DNA data

In this exercise set, we will start by aligning the DNA data we inspected on Monday as FASTQ files. We will also inspect some additional quality metrics from the aligned data. After alignment, we will proceed to germline and somatic variant calling in the second part of the exercise set, in ´variant_calling.md´.

## Running BWA

To align our DNA data to a reference genome, we will be using the [Burrows-Wheeler Aligner (BWA)](https://bio-bwa.sourceforge.net/). BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome.

```bash
# Go to your student folder and make a folder for alignments
cd /nfs/course/students/YOUR_KI_ID/data
mkdir dna_alignments

# First, align the normal sample
# Remember to replace YOUR_KI_ID below with your ID
# Runtime: ~ 10 min
bwa mem -t 2 -Y -R "@RG\tID:Exome_Norm\tPL:ILLUMINA\tPU:C1TD1ACXX-CGATGT.7\tLB:exome_norm_lib1\tSM:HCC1395BL_DNA" -o /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Norm.sam /nfs/course/students/YOUR_KI_ID/references/genome/ref_genome.fa /nfs/course/students/YOUR_KI_ID/data/fastq/Exome_Norm/Exome_Norm_R1.fastq.gz /nfs/course/students/YOUR_KI_ID/data/fastq/Exome_Norm/Exome_Norm_R2.fastq.gz

# Second, align the tumor sample
# Remember to replace YOUR_KI_ID below with your ID
# Runtime: ~ 15 min
bwa mem -t 2 -Y -R "@RG\tID:Exome_Tumor\tPL:ILLUMINA\tPU:C1TD1ACXX-ATCACG.7\tLB:exome_tumor_lib1\tSM:HCC1395_DNA" -o /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Tumor.sam /nfs/course/students/YOUR_KI_ID/references/genome/ref_genome.fa /nfs/course/students/YOUR_KI_ID/data/fastq/Exome_Tumor/Exome_Tumor_R1.fastq.gz /nfs/course/students/YOUR_KI_ID/data/fastq/Exome_Tumor/Exome_Tumor_R2.fastq.gz
```

## Converting SAM to BAM format
The output of the BWA alignment is a SAM file, which is a human-readable alignment file. Detailed information on the SAM format can be found [HERE](https://samtools.github.io/hts-specs/SAMv1.pdf). Generally, we want to convert our alignment files into the BAM format, which contains the same information as the SAM file but is binary, therefore taking up less space on our servers and computers. We perform this conversion from SAM to BAM below.

```bash
# Make sure you are in the folder with the SAM files
cd /nfs/course/students/YOUR_KI_ID/data/dna_alignments

samtools view -@ 2 -h -b -o Exome_Norm.bam Exome_Norm.sam
samtools view -@ 2 -h -b -o Exome_Tumor.bam Exome_Tumor.sam

# Have a look at the BAM file header and try to figure out what type of information it contains:
samtools view -H Exome_Tumor.bam | less -SN

# Have a look at the BAM file data and try to figure out what type of information is in the BAM file:
samtools view Exome_Tumor.bam | less -SN 
```

## Sort the BAM files
Below, we sort our BAM files by query or read name. Reads that have the same query name are derived from the same template, so this ensures that read pairs are grouped together in the sorted BAM file. This is needed for duplicate marking in the next section.
```bash
cd /nfs/course/students/YOUR_KI_ID/data/dna_alignments

# Runtime: ~ 4 min
picard SortSam I=Exome_Norm.bam O=Exome_Norm_namesorted.bam SO=queryname
picard SortSam I=Exome_Tumor.bam O=Exome_Tumor_namesorted.bam SO=queryname
```

## Mark duplicates
During library preparation, PCR is performed. Removing the PCR duplicates is necessary for substantial removal of PCR artefacts. Picard MarkDuplicates uses the names of the sequencing reads and alignment positions to identify PCR duplicates. 

```bash
cd /nfs/course/students/YOUR_KI_ID/data/dna_alignments

picard MarkDuplicates -I Exome_Norm_namesorted.bam -O Exome_Norm_namesorted_mrkdup.bam -ASSUME_SORT_ORDER queryname -METRICS_FILE Exome_Norm_mrkdup_metrics.txt -QUIET true -COMPRESSION_LEVEL 0 -VALIDATION_STRINGENCY LENIENT
picard MarkDuplicates -I Exome_Tumor_namesorted.bam -O Exome_Tumor_namesorted_mrkdup.bam -ASSUME_SORT_ORDER queryname -METRICS_FILE Exome_Tumor_mrkdup_metrics.txt -QUIET true -COMPRESSION_LEVEL 0 -VALIDATION_STRINGENCY LENIENT
```
## Position sort the BAM files

For the subsequent steps of our DNA analysis, a position-sorted BAM file is required. Therefore, we will sort the BAM files by coordinate.

```bash
cd /nfs/course/students/YOUR_KI_ID/data/dna_alignments

picard SortSam -I Exome_Norm_namesorted_mrkdup.bam -O Exome_Norm_sorted_mrkdup.bam -SO coordinate
picard SortSam -I Exome_Tumor_namesorted_mrkdup.bam -O Exome_Tumor_sorted_mrkdup.bam -SO coordinate
```
## Create a BAM index for use with GATK, IGV, and other tools

In order to efficiently load and search a BAM file, downstream applications typically require an index file. We will create this file next.

```bash
cd /nfs/course/students/YOUR_KI_ID/data/dna_alignments

picard BuildBamIndex -I Exome_Norm_sorted_mrkdup.bam
picard BuildBamIndex -I Exome_Tumor_sorted_mrkdup.bam
```

## Indel realignment

Depending on the downstream approach, a process called indel realignment may be needed. During alignment of the DNA data, each individual read is aligned to the reference genome individually. This may lead to misaligned reads in e.g. repetitive regions or due to somatic/germline variation. Documentation for the indel realignment step is [here](https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tutorials/(howto)_Perform_local_realignment_around_indels.md). In this course, HaplotypeCaller/MuTect2 will be applied to identify germline/somatic variation. These two variant callers do realignment internally and therefore we do not need to perform it here.

## Perform Base Quality Score Recalibration (BQSR)

In a nutshell, Base Quality Score Recalibration (BQSR) is a data pre-processing step that detects systematic errors made by the sequencing machine when it estimates the accuracy of each base call. This is the next step of our processing of the DNA data.

### Calculate BQSR Table

First, we calculate the BQSR table.

NOTE: For exome data, we should usually use the `--intervals` (`-L`) option in the commands below. This option excludes off-target sequences and sequences that may be poorly mapped, which have a higher error rate. Including them could lead to a skewed model and bad recalibration. **We do not use this option in this exercise for simplicity, but in a real analysis it should be included.**

```bash
cd /nfs/course/students/YOUR_KI_ID/data/dna_alignments

gatk --java-options '-Xmx12g' BaseRecalibrator -R /nfs/course/students/YOUR_KI_ID/references/genome/ref_genome.fa -I /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Norm_sorted_mrkdup.bam -O /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Norm_sorted_mrkdup_bqsr.table --known-sites /nfs/course/students/YOUR_KI_ID/references/gatk/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites /nfs/course/students/YOUR_KI_ID/references/gatk/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites /nfs/course/students/YOUR_KI_ID/references/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --preserve-qscores-less-than 6 --disable-bam-index-caching

gatk --java-options '-Xmx12g' BaseRecalibrator -R /nfs/course/students/YOUR_KI_ID/references/genome/ref_genome.fa -I /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Tumor_sorted_mrkdup.bam -O /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Tumor_sorted_mrkdup_bqsr.table --known-sites /nfs/course/students/YOUR_KI_ID/references/gatk/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites /nfs/course/students/YOUR_KI_ID/references/gatk/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites /nfs/course/students/YOUR_KI_ID/references/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --preserve-qscores-less-than 6 --disable-bam-index-caching
```

### Apply BQSR
Now, apply the BQSR table to the BAM files.
```bash
cd /nfs/course/students/YOUR_KI_ID/data/dna_alignments

gatk --java-options '-Xmx12g' ApplyBQSR -R /nfs/course/students/YOUR_KI_ID/references/genome/ref_genome.fa -I /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Norm_sorted_mrkdup.bam -O /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Norm_sorted_mrkdup_bqsr.bam --bqsr-recal-file /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Norm_sorted_mrkdup_bqsr.table --preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30

gatk --java-options '-Xmx12g' ApplyBQSR -R /nfs/course/students/YOUR_KI_ID/references/genome/ref_genome.fa -I /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Tumor_sorted_mrkdup.bam -O /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Tumor_sorted_mrkdup_bqsr.bam --bqsr-recal-file /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Tumor_sorted_mrkdup_bqsr.table --preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30
```

Create an index for these new BAM files.

```bash
cd /nfs/course/students/YOUR_KI_ID/data/dna_alignments

picard BuildBamIndex I=Exome_Norm_sorted_mrkdup_bqsr.bam
picard BuildBamIndex I=Exome_Tumor_sorted_mrkdup_bqsr.bam
```

## Clean up the SAM/BAM files we no longer need

Keep the final sorted, duplicate-marked, BQSR BAM/bai/table files and the mrkdup.txt files. Delete everything else.

```bash
cd /nfs/course/students/YOUR_KI_ID/data/dna_alignments
mkdir final
mv *_sorted_mrkdup_bqsr.* final/
mv *.txt final/

rm *.sam
rm *.bam
rm *.bai

mv final/* .
rmdir final/
```

Next, we will look at some quality control metrics from the BAM files.

## Run the samtools flagstat tool

`samtools flagstat` provides QC-metrics after aligning the reads, e.g. the fraction of read-pairs that map to different chromosomes (which should be limited for a good quality alignment).

```bash
cd /nfs/course/students/YOUR_KI_ID/data/dna_alignments

samtools flagstat /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Norm_sorted_mrkdup_bqsr.bam > /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Norm_flagstat.txt
samtools flagstat /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Tumor_sorted_mrkdup_bqsr.bam > /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Tumor_flagstat.txt
```
**Question 1:** Look inside the `Exome_Tumor_flagstat.txt` file and refer to the `samtools flagstat` documentation [here](https://www.htslib.org/doc/samtools-flagstat.html). 
* How many total reads are in the file?
* How many reads are marked as duplicates?
* What percentage of reads had a mate mapped to a different chromosome?

## Run various picard CollectMetrics tools

Picard is a widely used quality control tool in genomics and it is named after a famous [Star Trek Captain](https://en.wikipedia.org/wiki/Jean-Luc_Picard).

Descriptions of the outputs of picard tools can be found [here](http://broadinstitute.github.io/picard/). 

```bash
# Create interval files for our exome BED files we created during Exercise Set 1
cd /nfs/course/students/YOUR_KI_ID/references/exome
picard BedToIntervalList I=exome_regions.bed O=exome_regions.bed.interval_list SD=/nfs/course/students/YOUR_KI_ID/references/genome/ref_genome.dict
picard BedToIntervalList I=probe_regions.bed O=probe_regions.bed.interval_list SD=/nfs/course/students/YOUR_KI_ID/references/genome/ref_genome.dict

# Then proceed to run Picard tools
cd /nfs/course/students/YOUR_KI_ID/data/dna_alignments

# Picard CollectInsertSizeMetrics
picard CollectInsertSizeMetrics -I /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Norm_sorted_mrkdup_bqsr.bam -O /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Norm_insert_size_metrics.txt -H /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Norm_insert_size_metrics.pdf
picard CollectInsertSizeMetrics -I /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Tumor_sorted_mrkdup_bqsr.bam -O /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Tumor_insert_size_metrics.txt -H /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Tumor_insert_size_metrics.pdf

# Picard CollectAlignmentSummaryMetrics
picard CollectAlignmentSummaryMetrics -I /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Norm_sorted_mrkdup_bqsr.bam -O /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Norm_alignment_metrics.txt -R /nfs/course/students/YOUR_KI_ID/references/genome/ref_genome.fa
picard CollectAlignmentSummaryMetrics -I /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Tumor_sorted_mrkdup_bqsr.bam -O /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Tumor_exome_alignment_metrics.txt -R /nfs/course/students/YOUR_KI_ID/references/genome/ref_genome.fa

#Picard CollectHsMetrics
#Exome only
picard CollectHsMetrics -I /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Norm_sorted_mrkdup_bqsr.bam -O /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Norm_hs_metrics.txt -R /nfs/course/students/YOUR_KI_ID/references/genome/ref_genome.fa -BI /nfs/course/students/YOUR_KI_ID/references/exome/probe_regions.bed.interval_list -TI /nfs/course/students/YOUR_KI_ID/references/exome/exome_regions.bed.interval_list
picard CollectHsMetrics -I /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Tumor_sorted_mrkdup_bqsr.bam -O /nfs/course/students/YOUR_KI_ID/data/dna_alignments/Exome_Tumor_hs_metrics.txt -R /nfs/course/students/YOUR_KI_ID/references/genome/ref_genome.fa -BI /nfs/course/students/YOUR_KI_ID/references/exome/probe_regions.bed.interval_list -TI /nfs/course/students/YOUR_KI_ID/references/exome/exome_regions.bed.interval_list
```

## Run FastQC

```bash
cd /nfs/course/students/YOUR_KI_ID/data/dna_alignments
fastqc -t 8 Exome_Norm_sorted_mrkdup_bqsr.bam
fastqc -t 8 Exome_Tumor_sorted_mrkdup_bqsr.bam
tree
```

## Run MultiQC to produce a final QC report
```bash
cd /nfs/course/students/YOUR_KI_ID/data/dna_alignments
mkdir post_align_qc
cd post_align_qc
multiqc /nfs/course/students/YOUR_KI_ID/data/dna_alignments
tree

# Download the MultiQC output to your local computer and open the .html in your favorite browser
scp YOUR_KI_ID@c8cancergen02.ki.se:/nfs/course/students/YOUR_KI_ID/data/dna_alignments/post_align_qc/multiqc_report.html .
```

**Question 2:** Look through the MultiQC report and answer the following:
* Look at the Picard Alignment Summary. What percentage of reads in the Tumor and Normal were unaligned? Do you think this is a good quality alignment result?
* Look at the Picard Mark Duplicates plot. What was the most common type of duplicate that was marked? (*FYI- optical duplicates are reads that result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument.*)