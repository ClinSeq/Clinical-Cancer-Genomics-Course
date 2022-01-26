# Bioinformatic processing of DNA data

```bash
cd ~/workspace/inputs/data/fastq

# list all files
tree

# view the exome normal sample data files
tree Exome_Norm

# view the RNA-seq tumor sample data files.
tree RNAseq_Tumor
```
## Investigate the raw data files (fastq-files)
```bash

cd ~/workspace/inputs/data/fastq/Exome_Tumor

# show the first ten lines of the Exome Tumor fastq files
# note this is how paird-end data is stored. 
# Read 1
zcat Exome_Tumor/Exome_Tumor_R1.fastq.gz | head
# Read 2
zcat Exome_Tumor/Exome_Tumor_R2.fastq.gz | head

# This wikipedia file gives a very nice overview of fastq-files and Illumina base qualities
# https://en.wikipedia.org/wiki/FASTQ_format
zcat Exome_Tumor/Exome_Tumor_R1.fastq.gz | head -n 8

# what do R1 and R2 refer to? What is the length of each read?
zcat Exome_Tumor/Exome_Tumor_R1.fastq.gz | head -n 2 | tail -n 1 | wc

# how many lines are there in the Exome_Tumor file
zcat Exome_Tumor/Exome_Tumor_R1.fastq.gz | wc -l # There are: 33,326,620

# how many paired reads or fragments are there then?
expr 33326620 / 4 # There are: 8,331,655 paired end reads

# how many total bases of data are in the Exome Tumor data set?
echo "8331655 * (101 * 2)" | bc # There are: 1,682,994,310 bases of data

# how many total bases when expressed as "gigabases" (specify 2 decimal points using `scale`)
echo "scale=2; (8331655 * (101 * 2))/1000000000" | bc # There are: 1.68 Gbp of data

# what is the average coverage we expect to achieve with this much data for the exome region targeted?

# first determine the size of our exome regions (answer = 6683920). 
cat /workspace/inputs/references/exome/exome_regions.bed | perl -ne 'chomp; @l=split("\t", $_); $size += $l[2]-$l[1]; if (eof){print "size = $size\n"}' 

# now determine the average coverage of these positions by our bases of data
echo "scale=2; (8331655 * (101 * 2))/6683920" | bc # Average covered expected = 251.79x

# what is the fundamental assumption of this calculation that is at least partially not true? What effect will this have on the observed coverage?
```
## Run fastqc to check the quality of the fastq-files
Have a look [here](https://www.youtube.com/watch?v=lUk5Ju3vCDM) for a short tutorial on the fastqc output.
```bash
cd ~/workspace/inputs/data/fastq

fastqc -t 8 Exome_Norm/Exome_Norm*.fastq.gz
fastqc -t 8 Exome_Tumor/Exome_Tumor*.fastq.gz
tree

fastqc -t 8 RNAseq_Norm/RNAseq_Norm*.fastq.gz
fastqc -t 8 RNAseq_Tumor/RNAseq_Tumor*.fastq.gz
tree
```

## Run MultiQC 
We will run mutliQC to provide a combined report of the fastQC data. 

More info on MultiQC is available [here](https://multiqc.info).

```bash
cd ~/workspace/inputs
#Dont do if directory exist
#mkdir qc
cd qc
multiqc ~/workspace/inputs/data/fastq/
tree
# Download MultiQC output to the local computer, open the .html in you favourite browser.
scp -ri ~/PEM_KEY_ID.pem ubuntu@AWS_ADDRESS_HERE:~/workspace/inputs/qc/multiqc* .
```

## Run BWA

```bash
# Run bwa mem using the above information

cd ~/workspace/align
# Runtime: ~ 4min
bwa mem -t 2 -Y -R "@RG\tID:Exome_Norm\tPL:ILLUMINA\tPU:C1TD1ACXX-CGATGT.7\tLB:exome_norm_lib1\tSM:HCC1395BL_DNA" -o ~/workspace/align/Exome_Norm.sam ~/workspace/inputs/references/genome/ref_genome.fa ~/workspace/inputs/data/fastq/Exome_Norm/Exome_Norm_R1.fastq.gz ~/workspace/inputs/data/fastq/Exome_Norm/Exome_Norm_R2.fastq.gz
# Runtime: ~ 4min
bwa mem -t 2 -Y -R "@RG\tID:Exome_Tumor\tPL:ILLUMINA\tPU:C1TD1ACXX-ATCACG.7\tLB:exome_tumor_lib1\tSM:HCC1395_DNA" -o ~/workspace/align/Exome_Tumor.sam ~/workspace/inputs/references/genome/ref_genome.fa ~/workspace/inputs/data/fastq/Exome_Tumor/Exome_Tumor_R1.fastq.gz ~/workspace/inputs/data/fastq/Exome_Tumor/Exome_Tumor_R2.fastq.gz
```

## Convert sam to bam format
```bash
cd ~/workspace/align
# Runtime: ~30s
samtools view -@ 2 -h -b -o Exome_Norm.bam Exome_Norm.sam
# Runtime: ~45s
samtools view -@ 2 -h -b -o Exome_Tumor.bam Exome_Tumor.sam

# Have a look at the bam header and try to figure out what type of information that is in the bam file:
samtools view -H Exome_Tumor.bam | less -SN 
# Have a look at the bam file data and try to figure out what type of information that is in the bam file:
samtools view Exome_Tumor.bam | less -SN 
```
Detailed information on the sam format can be found [HERE](https://samtools.github.io/hts-specs/SAMv1.pdf) .

## Sort bam files
```bash
cd ~/workspace/align
# Runtime: ~ 4min
java -Xmx12g -jar $PICARD SortSam I=Exome_Norm.bam O=Exome_Norm_namesorted.bam SO=queryname
java -Xmx12g -jar $PICARD SortSam I=Exome_Tumor.bam O=Exome_Tumor_namesorted.bam SO=queryname
```

## Mark duplicates
During library preparation, PCR is performed. Removing the PCR duplicates is necessary for substantial removal of PCR artefacts. Picard MarkDuplicates use the read names and alignment positions to identify PCR duplicates. 

```bash
cd ~/workspace/align
java -Xmx12g -jar $PICARD MarkDuplicates -I Exome_Norm_namesorted.bam -O Exome_Norm_namesorted_mrkdup.bam -ASSUME_SORT_ORDER queryname -METRICS_FILE Exome_Norm_mrkdup_metrics.txt -QUIET true -COMPRESSION_LEVEL 0 -VALIDATION_STRINGENCY LENIENT
java -Xmx12g -jar $PICARD MarkDuplicates -I Exome_Tumor_namesorted.bam -O Exome_Tumor_namesorted_mrkdup.bam -ASSUME_SORT_ORDER queryname -METRICS_FILE Exome_Tumor_mrkdup_metrics.txt -QUIET true -COMPRESSION_LEVEL 0 -VALIDATION_STRINGENCY LENIENT
```
## Position sort bam file

For indexing and subsequent steps a position-sorted bam is required. Therefore, we will sort bam file by coordinate.

```bash
 cd ~/workspace/align
java -Xmx12g -jar $PICARD SortSam -I Exome_Norm_namesorted_mrkdup.bam -O Exome_Norm_sorted_mrkdup.bam -SO coordinate
java -Xmx12g -jar $PICARD SortSam -I Exome_Tumor_namesorted_mrkdup.bam -O Exome_Tumor_sorted_mrkdup.bam -SO coordinate
```
## Create bam index for use with GATK, IGV, etc

In order to efficiently load and search a bam file, downstream applications typically require an index

```bash
cd ~/workspace/align
java -Xmx12g -jar $PICARD BuildBamIndex -I Exome_Norm_sorted_mrkdup.bam
java -Xmx12g -jar $PICARD BuildBamIndex -I Exome_Tumor_sorted_mrkdup.bam
```

## Indel Realignment

Depending on the downstream appoarch something called Indel realignment may be needed. During alignment of the DNA data, each individual read is aligned to the reference sequence individually. This may lead to misaligned reads in e.g. repetitive regions or due to somatic/germline variation. If needed, add the realignment-step here, see docs [here](https://software.broadinstitute.org/gatk/documentation/article?id=7156). In this course the HaplotypeCaller/MuTect2 will be applied to identify germline/somatic variantion. These two variant callers do realignment internally and therefore it is not needed. See the following [blog](https://software.broadinstitute.org/gatk/blog?id=7847) for a detailed discussion of this issue. See [here](https://drive.google.com/drive/folders/1U6Zm_tYn_3yeEgrD1bdxye4SXf5OseIt) for latest versions of all GATK tutorials.

## Perform Base Quality Score Recalibration

BQSR stands for Base Quality Score Recalibration. In a nutshell, it is a data pre-processing step that detects systematic errors made by the sequencing machine when it estimates the accuracy of each base call.

### Calculate BQSR Table

First, calculate the BQSR table.

NOTE: For exome data, we should modify the below to use the `--intervals` (`-L`) option. This excludes off-target sequences and sequences that may be poorly mapped, which have a higher error rate. Including them could lead to a skewed model and bad recalibration.

```bash
cd ~/workspace/align
gatk --java-options '-Xmx12g' BaseRecalibrator -R ~/workspace/inputs/references/genome/ref_genome.fa -I ~/workspace/align/Exome_Norm_sorted_mrkdup.bam -O ~/workspace/align/Exome_Norm_sorted_mrkdup_bqsr.table --known-sites ~/workspace/inputs/references/gatk/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites ~/workspace/inputs/references/gatk/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites ~/workspace/inputs/references/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --preserve-qscores-less-than 6 --disable-bam-index-caching $GATK_REGIONS

gatk --java-options '-Xmx12g' BaseRecalibrator -R ~/workspace/inputs/references/genome/ref_genome.fa -I ~/workspace/align/Exome_Tumor_sorted_mrkdup.bam -O ~/workspace/align/Exome_Tumor_sorted_mrkdup_bqsr.table --known-sites ~/workspace/inputs/references/gatk/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites ~/workspace/inputs/references/gatk/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites ~/workspace/inputs/references/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --preserve-qscores-less-than 6 --disable-bam-index-caching $GATK_REGIONS

```

### Apply BQSR
Now, apply the BQSR table to bam files.
```bash
cd ~/workspace/align
gatk --java-options '-Xmx12g' ApplyBQSR -R ~/workspace/inputs/references/genome/ref_genome.fa -I ~/workspace/align/Exome_Norm_sorted_mrkdup.bam -O ~/workspace/align/Exome_Norm_sorted_mrkdup_bqsr.bam --bqsr-recal-file ~/workspace/align/Exome_Norm_sorted_mrkdup_bqsr.table --preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30

gatk --java-options '-Xmx12g' ApplyBQSR -R ~/workspace/inputs/references/genome/ref_genome.fa -I ~/workspace/align/Exome_Tumor_sorted_mrkdup.bam -O ~/workspace/align/Exome_Tumor_sorted_mrkdup_bqsr.bam --bqsr-recal-file ~/workspace/align/Exome_Tumor_sorted_mrkdup_bqsr.table --preserve-qscores-less-than 6 --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30
```

Create an index for these new bams

```bash
cd ~/workspace/align
java -Xmx12g -jar $PICARD BuildBamIndex I=Exome_Norm_sorted_mrkdup_bqsr.bam
java -Xmx12g -jar $PICARD BuildBamIndex I=Exome_Tumor_sorted_mrkdup_bqsr.bam
```

## Clean up un-needed sam/bam files

Keep final sorted, duplicated marked, bqsr bam/bai/table files and mrkdup.txt files. Delete everything else.

```bash
cd ~/workspace/align
mkdir final
mv *_sorted_mrkdup_bqsr.* final/
mv *.txt final/

rm *.sam
rm *.bam
rm *.bai

mv final/* .
rmdir final/
```

## Run Samtools flagstat

Samtools flagstat provides QC-metrics after aligning the reads, e.g. the fraction of read-pairs that map to different chromosomes (which should be limited).

```bash
cd ~/workspace/align/
samtools flagstat ~/workspace/align/Exome_Norm_sorted_mrkdup_bqsr.bam > ~/workspace/align/Exome_Norm_flagstat.txt
samtools flagstat ~/workspace/align/Exome_Tumor_sorted_mrkdup_bqsr.bam > ~/workspace/align/Exome_Tumor_flagstat.txt
```
## Run various picard CollectMetrics tools

Picard is a widely used QC tool in genomics and it is named after a famous [Star Trek Captain](https://en.wikipedia.org/wiki/Jean-Luc_Picard).

Description of the output of picard can be found [here](http://broadinstitute.github.io/picard/). 

```bash
cd ~/workspace/align/
java -Xmx12g -jar $PICARD CollectInsertSizeMetrics -I ~/workspace/align/Exome_Norm_sorted_mrkdup_bqsr.bam -O ~/workspace/align/Exome_Norm_insert_size_metrics.txt -H ~/workspace/align/Exome_Norm_insert_size_metrics.pdf
java -Xmx12g -jar $PICARD CollectInsertSizeMetrics -I ~/workspace/align/Exome_Tumor_sorted_mrkdup_bqsr.bam -O ~/workspace/align/Exome_Tumor_insert_size_metrics.txt -H ~/workspace/align/Exome_Tumor_insert_size_metrics.pdf

java -Xmx12g -jar $PICARD CollectAlignmentSummaryMetrics -I ~/workspace/align/Exome_Norm_sorted_mrkdup_bqsr.bam -O ~/workspace/align/Exome_Norm_alignment_metrics.txt -R ~/workspace/inputs/references/genome/ref_genome.fa
java -Xmx12g -jar $PICARD CollectAlignmentSummaryMetrics -I ~/workspace/align/Exome_Tumor_sorted_mrkdup_bqsr.bam -O ~/workspace/align/Exome_Tumor_exome_alignment_metrics.txt -R ~/workspace/inputs/references/genome/ref_genome.fa

#Picard CollectHsMetrics
#Exome Only
java -Xmx12g -jar $PICARD CollectHsMetrics -I ~/workspace/align/Exome_Norm_sorted_mrkdup_bqsr.bam -O ~/workspace/align/Exome_Norm_hs_metrics.txt -R ~/workspace/inputs/references/genome/ref_genome.fa -BI ~/workspace/inputs/references/exome/probe_regions.bed.interval_list -TI ~/workspace/inputs/references/exome/exome_regions.bed.interval_list
java -Xmx12g -jar $PICARD CollectHsMetrics -I ~/workspace/align/Exome_Tumor_sorted_mrkdup_bqsr.bam -O ~/workspace/align/Exome_Tumor_hs_metrics.txt -R ~/workspace/inputs/references/genome/ref_genome.fa -BI ~/workspace/inputs/references/exome/probe_regions.bed.interval_list -TI ~/workspace/inputs/references/exome/exome_regions.bed.interval_list
```

## Run FastQC

```bash
cd ~/workspace/align
fastqc -t 8 Exome_Norm_sorted_mrkdup_bqsr.bam
fastqc -t 8 Exome_Tumor_sorted_mrkdup_bqsr.bam
tree
```
## Run MultiQC to produce a final report
```bash
cd ~/workspace/align
#Dont do if directory exist
#mkdir post_align_qc
cd post_align_qc
multiqc ~/workspace/align/
tree

# Download MultiQC output to the local computer, open the .html in you favourite browser.

# Discuss the MultiQC output with a fellow student
scp -ri ~/PEM_KEY_ID.pem ubuntu@AWS_ADDRESS_HERE:~/workspace/align/post_align_qc/multiqc* .
```

