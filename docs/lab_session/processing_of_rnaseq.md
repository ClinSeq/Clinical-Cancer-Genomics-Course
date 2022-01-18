# BIOINFORMATIC PROCESSING OF RNA DATA 

## Adapter Trimming of the FASTQ files

The purpose of adapter trimming is to remove sequences in our data that correspond to the Illumina sequence adapters.  The most common adapter trimming scenario is the removal of adapter sequences that occur at the end of read sequences. This happens when a DNA (or cDNA) fragment is shorter than the read length.  For example if we sequence our RNA-seq fragments to 150 base length and a fragment is only 140 bases long the read will end with 10 bases of adapter sequence. Since this adapter sequence does not correspond to the genome, it will not align. Too much adapter sequence can actually prevent reads from aligning at all. Adapter trimming may therefore sometime improve the overall alignment success rate for an RNA-seq data set.  Adapter trimming involves a simplistic alignment itself and therefore can be computationally expensive.

```bash
cd ~/workspace/inputs/references

#Command already run
#wget -c http://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa

#Have a look at the adapter sequences that we want to remove from the RNA data
less -SN ~/workspace/inputs/references/illumina_multiplex.fa

#Start with tumor RNAseq data
cd ~/workspace/inputs/data/fastq/RNAseq_Tumor

#Have a look to see what is in the folder
ls -halt

#pPerform trimming using flexbar

#Before using a new tool, make a habit of looking the tool up,
# google "flexbar trimming adapters" or something similar

#As this step takes approxmately 25 min, we will use nohup. 
#Google "nohup" ...
#If running with nohup the command will continue even if the connection breaks to the aws server

#Also, notice that we are using the flexbar argument --threads 7 which means that we are using 7 cpu threads. To check how many you have at your disposal use the command:
lscpu
#We have 2 thread left which is good, not to hog all resources
#To check the available RAM use the command:
cat /proc/meminfo

#When trimming the adapters run two simultaneous jobs using 7 cores/job and send both to the background if connection breaks. Redirect the standard output to a custom log file
#nohup cmd > custom-out.log &

nohup flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters ~/workspace/inputs/references/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 7 --zip-output GZ --reads RNAseq_Tumor_Lane1_R1.fastq.gz --reads2 RNAseq_Tumor_Lane1_R2.fastq.gz --target ~/workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane1 > RNAseq_Tumor_Lane1.log &

#Start the second process
nohup flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters ~/workspace/inputs/references/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 7 --zip-output GZ --reads RNAseq_Tumor_Lane2_R1.fastq.gz --reads2 RNAseq_Tumor_Lane2_R2.fastq.gz --target ~/workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane2 > RNAseq_Tumor_Lane2.log &

#Use htop to monitor the processor load
htop
#quit htop by pressing "q"

#have a look at one of the log files
less -SN RNAseq_Tumor_Lane1.log

####### Wait for the processes to stop before continuing #######

#Start the 
#Normal RNAseq data
cd ~/workspace/inputs/data/fastq/RNAseq_Norm

nohup flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters ~/workspace/inputs/references/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 7 --zip-output GZ --reads RNAseq_Norm_Lane1_R1.fastq.gz --reads2 RNAseq_Norm_Lane1_R2.fastq.gz --target ~/workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane1 > RNAseq_Normal_Lane1.log &

flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters ~/workspace/inputs/references/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 7 --zip-output GZ --reads RNAseq_Norm_Lane2_R1.fastq.gz --reads2 RNAseq_Norm_Lane2_R2.fastq.gz --target ~/workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane2 > RNAseq_Normal_Lane2.log &

####### Wait for the processes to stop #######

```
As you can see, nohup is very useful for sending processes to the backgroud that will continue if connection breaks or if you need to log out.
For repetition how to send stdout and stderr to a file while running nohup, have a look [here](https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file).

## Run fastqc to check the quality of the fastq-files (as for DNA data)

```bash
#Run fastqc and compre the trimmed and non-trimmed file
cd ~/workspace/inputs/data/fastq/RNAseq_Tumor
fastqc RNAseq_Tumor_Lane1_R1.fastq.gz
fastqc RNAseq_Tumor_Lane1_1.fastq.gz

# Download the output to your computer
scp -ri ~/course_EC2_01.pem ubuntu@AWS_ADDRESS_HERE:~/workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane1_1_fastqc.html .
scp -ri ~/course_EC2_01.pem ubuntu@AWS_ADDRESS_HERE:~/workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane1_R1_fastqc.html .
```
Compare the two html files side by side:
- Any obivous difference?

```bash
#Move the trimmed data
mkdir -p ~/workspace/inputs/data/fastq/RNAseq_Tumor/trimmed
mkdir -p ~/workspace/inputs/data/fastq/RNAseq_Norm/trimmed
mv ~/workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane1_1.fastq.gz ~/workspace/inputs/data/fastq/RNAseq_Tumor/trimmed/RNAseq_Tumor_Lane1_R1.fastq.gz
mv ~/workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane1_2.fastq.gz ~/workspace/inputs/data/fastq/RNAseq_Tumor/trimmed/RNAseq_Tumor_Lane1_R2.fastq.gz
mv ~/workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane2_1.fastq.gz ~/workspace/inputs/data/fastq/RNAseq_Tumor/trimmed/RNAseq_Tumor_Lane2_R1.fastq.gz
mv ~/workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane2_2.fastq.gz ~/workspace/inputs/data/fastq/RNAseq_Tumor/trimmed/RNAseq_Tumor_Lane2_R2.fastq.gz

mv ~/workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane1_1.fastq.gz ~/workspace/inputs/data/fastq/RNAseq_Norm/trimmed/RNAseq_Norm_Lane1_R1.fastq.gz
mv ~/workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane1_2.fastq.gz ~/workspace/inputs/data/fastq/RNAseq_Norm/trimmed/RNAseq_Norm_Lane1_R2.fastq.gz
mv ~/workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane2_1.fastq.gz ~/workspace/inputs/data/fastq/RNAseq_Norm/trimmed/RNAseq_Norm_Lane2_R1.fastq.gz
mv ~/workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane2_2.fastq.gz ~/workspace/inputs/data/fastq/RNAseq_Norm/trimmed/RNAseq_Norm_Lane2_R2.fastq.gz

#make a habit of always checking the result of any command
ls -halt ~/workspace/inputs/data/fastq/RNAseq_Tumor/trimmed
ls -halt ~/workspace/inputs/data/fastq/RNAseq_Norm/trimmed
```

## Alignment

We will use the aligner [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) to perform spliced alignments to our reference genome. For efficiency, the output of HISAT (SAM format) will be piped directly to another program called [sambamba](http://lomereiter.github.io/sambamba/) to first convert to BAM format and then sort that BAM file. Before each command we will also create and assign a path for temporary directories.

```bash
mkdir -p ~/workspace/rnaseq/alignments
cd ~/workspace/rnaseq/alignments

# Align tumor data
# Runtime: ~15 min each run
TUMOR_DATA_1_TEMP=`mktemp -d ~/workspace/rnaseq/alignments/2895626107.XXXXXXXXXXXX`
nohup hisat2 -p 7 --dta -x ~/workspace/inputs/references/transcriptome/ref_genome --rg-id 2895626107 --rg PL:ILLUMINA --rg PU:H3MYFBBXX-GCCAAT.4 --rg LB:rna_tumor_lib1 --rg SM:HCC1395_RNA --rna-strandness RF -1 ~/workspace/inputs/data/fastq/RNAseq_Tumor/trimmed/RNAseq_Tumor_Lane1_R1.fastq.gz -2  ~/workspace/inputs/data/fastq/RNAseq_Tumor/trimmed/RNAseq_Tumor_Lane1_R2.fastq.gz | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t 8 -m 24G --tmpdir $TUMOR_DATA_1_TEMP -o ~/workspace/rnaseq/alignments/RNAseq_Tumor_Lane1.bam /dev/stdin &

TUMOR_DATA_2_TEMP=`mktemp -d ~/workspace/rnaseq/alignments/2895626112.XXXXXXXXXXXX`
nohup hisat2 -p 7 --dta -x ~/workspace/inputs/references/transcriptome/ref_genome --rg-id 2895626112 --rg PL:ILLUMINA --rg PU:H3MYFBBXX-GCCAAT.5 --rg LB:rna_tumor_lib1 --rg SM:HCC1395_RNA --rna-strandness RF -1 ~/workspace/inputs/data/fastq/RNAseq_Tumor/trimmed/RNAseq_Tumor_Lane2_R1.fastq.gz -2  ~/workspace/inputs/data/fastq/RNAseq_Tumor/trimmed/RNAseq_Tumor_Lane2_R2.fastq.gz | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t 8 -m 24G --tmpdir $TUMOR_DATA_2_TEMP -o ~/workspace/rnaseq/alignments/RNAseq_Tumor_Lane2.bam /dev/stdin &

#Wait for the processes to finish, monitor using htop
rmdir $TUMOR_DATA_2_TEMP $TUMOR_DATA_1_TEMP

# Align normal data

NORMAL_DATA_1_TEMP=`mktemp -d ~/workspace/rnaseq/alignments/2895625992.XXXXXXXXXXXX`
nohup hisat2 -p 7 --dta -x ~/workspace/inputs/references/transcriptome/ref_genome --rg-id 2895625992 --rg PL:ILLUMINA --rg PU:H3MYFBBXX-CTTGTA.4 --rg LB:rna_norm_lib1 --rg SM:HCC1395BL_RNA --rna-strandness RF -1 ~/workspace/inputs/data/fastq/RNAseq_Norm/trimmed/RNAseq_Norm_Lane1_R1.fastq.gz -2  ~/workspace/inputs/data/fastq/RNAseq_Norm/trimmed/RNAseq_Norm_Lane1_R2.fastq.gz | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t 8 -m 24G --tmpdir $NORMAL_DATA_1_TEMP -o ~/workspace/rnaseq/alignments/RNAseq_Norm_Lane1.bam /dev/stdin &

NORMAL_DATA_2_TEMP=`mktemp -d ~/workspace/rnaseq/alignments/2895626097.XXXXXXXXXXXX`

nohup hisat2 -p 7 --dta -x ~/workspace/inputs/references/transcriptome/ref_genome --rg-id 2895626097 --rg PL:ILLUMINA --rg PU:H3MYFBBXX-CTTGTA.5 --rg LB:rna_norm_lib1 --rg SM:HCC1395BL_RNA --rna-strandness RF -1 ~/workspace/inputs/data/fastq/RNAseq_Norm/trimmed/RNAseq_Norm_Lane2_R1.fastq.gz -2  ~/workspace/inputs/data/fastq/RNAseq_Norm/trimmed/RNAseq_Norm_Lane2_R2.fastq.gz | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort -t 8 -m 24G --tmpdir $NORMAL_DATA_2_TEMP -o ~/workspace/rnaseq/alignments/RNAseq_Norm_Lane2.bam /dev/stdin &

#Wait for the processes to finish, monitor using htop
rmdir $NORMAL_DATA_2_TEMP $NORMAL_DATA_1_TEMP
```
## Merging BAMs
Since we have multiple BAMs of each sample that just represent additional data for the same sequence library, we should combine them into a single BAM for convenience before proceeding.

```bash
cd ~/workspace/rnaseq/alignments
#Runtime: ~ 8m each merging command
sambamba merge -t 8 ~/workspace/rnaseq/alignments/RNAseq_Norm.bam ~/workspace/rnaseq/alignments/RNAseq_Norm_Lane1.bam ~/workspace/rnaseq/alignments/RNAseq_Norm_Lane2.bam

sambamba merge -t 8 ~/workspace/rnaseq/alignments/RNAseq_Tumor.bam ~/workspace/rnaseq/alignments/RNAseq_Tumor_Lane1.bam ~/workspace/rnaseq/alignments/RNAseq_Tumor_Lane2.bam
```

## GTF (General Transfer Format) and creating files for downstream processing
We will encounter the GTF file format during the exercises below, described in detail [here](https://www.ensembl.org/info/website/upload/gff.html#fields).

The GTF format is use to describe genes and other features of DNA, RNA and proteins.

In the section below we will prepare files needed for RNAseq QC

**NOTE: for the sake of time we are only investigating genes on chr6 and chr17**

```bash
#Go to this directory:
cd ~/workspace/inputs/references/transcriptome

#Have a peak into the transcriptome gtf file and try to make sense of it using the ensembl description above
less -SN ref_transcriptome.gtf

#check on which chromsomes from which we have transcripts
cut -f1 ref_transcriptome.gtf | sort | uniq -c

# Generating the necessary input files for picard CollectRnaSeqMetrics
cd ~/workspace/inputs/references/transcriptome

grep -i rrna ref_transcriptome.gtf > ref_ribosome.gtf

#A bed file with Ribosomal RNA
gff2bed < ~/workspace/inputs/references/transcriptome/ref_ribosome.gtf > ref_ribosome.bed

#Check the ribosome bed file
less -SN ref_ribosome.bed

java -jar $PICARD BedToIntervalList I=~/workspace/inputs/references/transcriptome/ref_ribosome.bed O=~/workspace/inputs/references/transcriptome/ref_ribosome.interval_list SD=~/workspace/inputs/references/genome/ref_genome.dict

gtfToGenePred -genePredExt ~/workspace/inputs/references/transcriptome/ref_transcriptome.gtf ~/workspace/inputs/references/transcriptome/ref_flat.txt

cat ref_flat.txt | awk '{print $12"\t"$0}' | cut -d$'\t' -f1-11 > ref_flat_final.txt

mv ref_flat_final.txt ref_flat.txt
```

## Post-alignment QC
```bash
cd ~/workspace/rnaseq/alignments
#Run samtools flagstat
#Runtime: ~2min, start and send to background
nohup samtools flagstat ~/workspace/rnaseq/alignments/RNAseq_Norm.bam > ~/workspace/rnaseq/alignments/RNAseq_Norm_flagstat.txt 2> flagstat_RNA_N_bam.out &
nohup samtools flagstat ~/workspace/rnaseq/alignments/RNAseq_Tumor.bam > ~/workspace/rnaseq/alignments/RNAseq_Tumor_flagstat.txt 2> flagstat_RNA_T_bam.out &

#Runtime: ~12 min, start and send to background
nohup fastqc -t 8 ~/workspace/rnaseq/alignments/RNAseq_Tumor.bam > fastqc_RNA_T_bam.out 2>&1 &
nohup fastqc -t 8 ~/workspace/rnaseq/alignments/RNAseq_Norm.bam > fastqc_RNA_N_bam.out 2>&1 &

# Runtime: 26min
nohup java -jar $PICARD CollectRnaSeqMetrics I=~/workspace/rnaseq/alignments/RNAseq_Norm.bam O=~/workspace/rnaseq/alignments/RNAseq_Norm.RNA_Metrics REF_FLAT=~/workspace/inputs/references/transcriptome/ref_flat.txt STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=~/workspace/inputs/references/transcriptome/ref_ribosome.interval_list > collectRnaSeqMetrics_N.out 2>&1 &
#Follow the progress of the program by looking in the nohup output file
less -SN collectRnaSeqMetrics_N.out

nohup java -jar $PICARD CollectRnaSeqMetrics I=~/workspace/rnaseq/alignments/RNAseq_Tumor.bam O=~/workspace/rnaseq/alignments/RNAseq_Tumor.RNA_Metrics REF_FLAT=~/workspace/inputs/references/transcriptome/ref_flat.txt STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=~/workspace/inputs/references/transcriptome/ref_ribosome.interval_list > collectRnaSeqMetrics_T.out 2>&1 &

cd ~/workspace/rnaseq
#mkdir post_align_qc
cd post_align_qc
multiqc ~/workspace/rnaseq/alignments/

#Finally, download multiqc files to your local computer
# Download MultiQC output to the local computer, open the .html in you favourite browser.
scp -ri ~/course_EC2_01.pem ubuntu@AWS_ADDRESS_HERE:~/workspace/rnaseq/post_align_qc/multiqc* .
```
Familiarize yourself with the RNAseq MultiQC data

## Indexing BAMs
In order to be able to view our BAM files in IGV, as usual we need to index them
```bash
cd  ~/workspace/rnaseq/alignments/
samtools index RNAseq_Norm.bam
samtools index RNAseq_Tumor.bam

```

## Run a simplified "reference only" StringTie expression approach

- The StringTie developer's recommend to: 
    - perform reference guided transcript compilation (aka transcript assembly) on each individual sample.
    - merge transcript predictions from all samples into a single model of the transcriptome.
    - annotate this predicted transcriptome with known transcriptome information.
    - estimate abundance for each of the transcripts in this final transcriptome model in each sample.
    
In the final result we would have abundance for all the same transcripts across all samples. This includes a combination of predicted and known transcripts. 

It is sometimes convenient to have a more simplified workflow where we only have values for known transcripts. This is particularly true in species where we already have comprehensive high quality transcriptome annotations and there is less of a focus on de novo transcript discovery.

The following workflow produces a "reference-only" transcriptome result in which we will perform abundance calculations on each lane of data individually.

```bash
cd ~/workspace/rnaseq
#mkdir ref-only-expression
cd ref-only-expression

nohup stringtie -p 4 -e -G ~/workspace/inputs/references/transcriptome/ref_transcriptome.gtf -o ~/workspace/rnaseq/ref-only-expression/RNAseq_Tumor_Lane1/transcripts.gtf -A ~/workspace/rnaseq/ref-only-expression/RNAseq_Tumor_Lane1/gene_abundances.tsv ~/workspace/rnaseq/alignments/RNAseq_Tumor_Lane1.bam > knowntranscripts_T1_L1.out 2>&1 &
nohup stringtie -p 4 -e -G ~/workspace/inputs/references/transcriptome/ref_transcriptome.gtf -o ~/workspace/rnaseq/ref-only-expression/RNAseq_Tumor_Lane2/transcripts.gtf -A ~/workspace/rnaseq/ref-only-expression/RNAseq_Tumor_Lane2/gene_abundances.tsv ~/workspace/rnaseq/alignments/RNAseq_Tumor_Lane2.bam > knowntranscripts_T1_L2.out 2>&1 &
nohup stringtie -p 4 -e -G ~/workspace/inputs/references/transcriptome/ref_transcriptome.gtf -o ~/workspace/rnaseq/ref-only-expression/RNAseq_Norm_Lane1/transcripts.gtf -A ~/workspace/rnaseq/ref-only-expression/RNAseq_Norm_Lane1/gene_abundances.tsv ~/workspace/rnaseq/alignments/RNAseq_Norm_Lane1.bam > knowntranscripts_N1_L1.out 2>&1 &
nohup stringtie -p 4 -e -G ~/workspace/inputs/references/transcriptome/ref_transcriptome.gtf -o ~/workspace/rnaseq/ref-only-expression/RNAseq_Norm_Lane2/transcripts.gtf -A ~/workspace/rnaseq/ref-only-expression/RNAseq_Norm_Lane2/gene_abundances.tsv ~/workspace/rnaseq/alignments/RNAseq_Norm_Lane2.bam > knowntranscripts_T2_L2.out 2>&1 &

#Wait to continue submitting jobs until the scripts have finished

```
Create a tidy expression matrix files for the StringTie results. This will be done at both the gene and transcript level and also will take into account the various expression measures produced: coverage, FPKM, and TPM.

```bash
cd ~/workspace/rnaseq/ref-only-expression

#Already run
#wget https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/scripts/stringtie_expression_matrix.pl
#chmod +x stringtie_expression_matrix.pl

./stringtie_expression_matrix.pl --expression_metric=TPM --result_dirs='RNAseq_Norm_Lane1,RNAseq_Norm_Lane2,RNAseq_Tumor_Lane1,RNAseq_Tumor_Lane2' --transcript_matrix_file=transcript_tpm_all_samples.tsv --gene_matrix_file=gene_tpm_all_samples.tsv

./stringtie_expression_matrix.pl --expression_metric=FPKM --result_dirs='RNAseq_Norm_Lane1,RNAseq_Norm_Lane2,RNAseq_Tumor_Lane1,RNAseq_Tumor_Lane2' --transcript_matrix_file=transcript_fpkm_all_samples.tsv --gene_matrix_file=gene_fpkm_all_samples.tsv

./stringtie_expression_matrix.pl --expression_metric=Coverage --result_dirs='RNAseq_Norm_Lane1,RNAseq_Norm_Lane2,RNAseq_Tumor_Lane1,RNAseq_Tumor_Lane2' --transcript_matrix_file=transcript_coverage_all_samples.tsv --gene_matrix_file=gene_coverage_all_samples.tsv

#Have a look at the output files
head gene_coverage_all_samples.tsv transcript_coverage_all_samples.tsv gene_fpkm_all_samples.tsv transcript_fpkm_all_samples.tsv gene_tpm_all_samples.tsv transcript_tpm_all_samples.tsv

```

## Reference Free Expression Analysis with Kallisto

Remember that in previous sections we have been using reference genome fasta sequences for the reference for alignment and subsequent steps. However, Kallisto works directly on target cDNA/transcript sequences. 

We have for stringtie used transcript annotations for genes on our subset of chromosomes (i.e. chr6 and chr17). The transcript models were downloaded from Ensembl in GTF format. This GTF contains a description of the coordinates of exons that make up each transcript but it does not contain the transcript sequences themselves which Kallisto is using. There are many places we could obtain such transcript sequences. For example, we could have download them directly in Fasta format from the Ensembl FTP site (or from UCSC or NCBI).

```bash
cd ~/workspace/rnaseq/
#mkdir kallisto
cd kallisto

# first check that the GTF and Fasta file are present for Kallisto
head ~/workspace/inputs/references/transcriptome/ref_transcriptome.gtf
head ~/workspace/inputs/references/transcriptome/kallisto/ref_transcriptome_clean.fa

#Each transcript and its sequence
head ~/workspace/inputs/references/transcriptome/kallisto/ref_transcriptome_clean.fa

# now check for the kallisto index is there
ls -halt ~/workspace/inputs/references/transcriptome/kallisto/ref_transcriptome_kallisto_index

# create a list of all transcript IDs for later use:
cd ~/workspace/rnaseq/kallisto/
cat ~/workspace/inputs/references/transcriptome/kallisto/ref_transcriptome_clean.fa | grep ">" | perl -ne '$_ =~ s/\>//; print $_' | sort | uniq > transcript_id_list.txt

head -n 10 transcript_id_list.txt
```

## Generate abundance estimates for all samples using Kallisto
As we did with StringTie we will generate transcript abundances for each of our demonstration samples using Kallisto. Here we are treating the two lanes for each sample as if they were independent samples.

```bash

cd ~/workspace/rnaseq/kallisto/
#mkdir quants
cd quants

nohup kallisto quant --index=/home/ubuntu/workspace/inputs/references/transcriptome/kallisto/ref_transcriptome_kallisto_index --output-dir=RNAseq_Norm_Lane1 --threads=8 --plaintext /home/ubuntu/workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane1_R1.fastq.gz /home/ubuntu/workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane1_R2.fastq.gz > kallisto_quant_N_L1.out 2>&1 &

nohup kallisto quant --index=/home/ubuntu/workspace/inputs/references/transcriptome/kallisto/ref_transcriptome_kallisto_index --output-dir=RNAseq_Norm_Lane2 --threads=8 --plaintext /home/ubuntu/workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane2_R1.fastq.gz /home/ubuntu/workspace/inputs/data/fastq/RNAseq_Norm/RNAseq_Norm_Lane2_R2.fastq.gz > kallisto_quant_N_L2.out 2>&1 &

nohup kallisto quant --index=/home/ubuntu/workspace/inputs/references/transcriptome/kallisto/ref_transcriptome_kallisto_index --output-dir=RNAseq_Tumor_Lane1 --threads=8 --plaintext /home/ubuntu/workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane1_R1.fastq.gz /home/ubuntu/workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane1_R2.fastq.gz > kallisto_quant_T_L1.out 2>&1 &

nohup kallisto quant --index=/home/ubuntu/workspace/inputs/references/transcriptome/kallisto/ref_transcriptome_kallisto_index --output-dir=RNAseq_Tumor_Lane2 --threads=8 --plaintext /home/ubuntu/workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane2_R1.fastq.gz /home/ubuntu/workspace/inputs/data/fastq/RNAseq_Tumor/RNAseq_Tumor_Lane2_R2.fastq.gz > kallisto_quant_T_L2.out 2>&1 &

#check nohup output
less -SN kallisto_quant_T_L2.out

#Create a single TSV file that has the TPM abundance estimates for all six samples.
#First check the contets of an abudance output file from Kallisto
less -SN ./RNAseq_Norm_Lane1/abundance.tsv

#Merge all files
paste */abundance.tsv | cut -f 1,2,5,10,15,20 > transcript_tpms_all_samples.tsv
ls -1 */abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\tlength$_\n"' > header.tsv

cat header.tsv transcript_tpms_all_samples.tsv | grep -v "tpm" > transcript_tpms_all_samples.tsv2

mv transcript_tpms_all_samples.tsv2 transcript_tpms_all_samples.tsv

rm -f header.tsv

head transcript_tpms_all_samples.tsv

#First create a gene version of the Kallisto TPM matrix, we will simply sum the TPM values for transcripts of the same gene.

# Command already run
#wget https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/scripts/kallisto_gene_matrix.pl
#chmod +x kallisto_gene_matrix.pl

./kallisto_gene_matrix.pl --gtf_file=/home/ubuntu/workspace/inputs/references/transcriptome/ref_transcriptome.gtf --kallisto_transcript_matrix_in=transcript_tpms_all_samples.tsv --kallisto_transcript_matrix_out=gene_tpms_all_samples.tsv

less -SN gene_tpms_all_samples.tsv
```

# Compare expression values between Kallisto and StringTie

To compare the two approaches we can use the expression value for each Ensembl transcript.

To do this comparison, we need to gather the expression estimates for each of our replicates from each approach. 

First - download the files to your local machine

```bash

#Run these commands on your local machine

scp -ri ~/course_EC2_01.pem ubuntu@AWS_ADDRESS_HERE:/home/ubuntu/workspace/rnaseq/kallisto/transcript_tpms_all_samples.tsv .

mv transcript_tpms_all_samples.tsv kallisto_transcript_tpms_all_samples.tsv

scp -ri ~/course_EC2_01.pem ubuntu@AWS_ADDRESS_HERE:/home/ubuntu/workspace/rnaseq/ref-only-expression/transcript_tpm_all_samples.tsv .
mv transcript_tpm_all_samples.tsv stringtie_transcript_tpms_all_samples.tsv

#Remember, one gene can have many transcripts. The ensembl gene ID for TP53 is ENSG00000141510. Check which ensembl transcripts that exist in the reference gtf.

#RUN ON THE AWS SERVER
grep ENST ~/workspace/inputs/references/transcriptome/ref_transcriptome.gtf | grep ENSG00000141510 | perl -nle '@a = split /\t/; $_ =~ /transcript_id\s"(\S*)";/g; print $1;' | sort | uniq

# start R

# set the working directory to where the files were downloaded
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
str(kallisto_stringtie_tpm)

#TP53 transcript ID vector from the transcript GTF
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

#To remove 0/NA values
kallisto_stringtie_tpm$value.stringtie[which(is.na(kallisto_stringtie_tpm$value.stringtie)==TRUE)] = 0.001
#kallisto_stringtie_tpm$value.kallisto[which(is.na(kallisto_stringtie_tpm$value.kallisto)==TRUE)] = 0.001
kallisto_stringtie_tpm$value.stringtie = kallisto_stringtie_tpm$value.stringtie + 0.1
kallisto_stringtie_tpm$value.kallisto = kallisto_stringtie_tpm$value.kallisto + 0.1
########### plot the result ######################

#pdf(file="transcript_stringtie_v_kallisto.pdf", height=8, width=8)
ggplot(data=kallisto_stringtie_tpm, aes(x=log2(value.kallisto), y=log2(value.stringtie), colour=TP53)) +
geom_point(alpha=1) +
scale_colour_manual(values = c("forestgreen", "firebrick", "dodgerblue"), name = "Transcripts", drop=TRUE) +
geom_point(data = subset(kallisto_stringtie_tpm, TP53 == "TP53"), aes(x=log2(value.kallisto), y=log2(value.stringtie), colour=TP53), inherit.aes=FALSE)+
facet_wrap(~variable)+
theme_bw()

#To plot density disbributions
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

# Gene Fusion detection

# Introduction

In addition to providing information about gene expression, RNA-seq data can be used to discover transcripts which result from chromosomal translocations. Translocations and their resultant chimeric (AKA fusion) transcripts are important driver mutations in many cancers. A [variety of specialized alignment and filtering strategies](https://www.ncbi.nlm.nih.gov/pubmed/27485475) have been developed to identify fusion transcripts from RNA, but these programs suffer from low specificity (i.e. many false-positives) and poor correlation across methods.

This tutorial uses the [kallisto](https://pachterlab.github.io/kallisto/about) and [pizzly](https://github.com/pmelsted/pizzly) tools for fusion detection from RNA-seq data. Kallisto quantifies transcript abundance through [pseudoalignment](http://tinyheero.github.io/2015/09/02/pseudoalignments-kallisto.html). Pizzly aligns reads which kallisto has flagged as potentially spanning fusion junctions. Running the tutorial requires RNA fastq files, a reference transcriptome, and a gene annotation file- see below.

Some of the tools used to detection gene fusions from RNAseq are Arriba, EricScript, FusionCatcher, Fusion-Inspector, fusion-report, Pizzly, Squid, Star-Fusion.

All these tools have their own merits and demerits. They are often used in combinations or all together. Hence, we need to assess which one is suitable for our data.

We will use Pizzly in this exercise.

# Setup

**Prerequisites- This module assumes you have completed the initial setup process for the course including:**
- installed Conda package manager, R, pizzly, and kallisto.
    - have fastq sequence files of normal and tumor RNA at ~/workspace/inputs/data/fastq/chr6_and_chr17/. 

**Additional setup:**

**_Important: pizzly will not work with the most recent Ensembl human GTF annotation file. Download the version 87 GTF as shown in the below code block. We will subset the reference transcriptome fasta file. Fasta splitting programs which do not preserve the full Ensembl header such as gffread or gtf_to_fasta will not work with pizzly._**

- Download Ensembl GTF and fasta and parse to include only chromosomes 6 and 17 (**15 min**): 

```bash
# Get files from source
#Already run
#mkdir -p ~/workspace/inputs/reference/fusion
#cd ~/workspace/inputs/reference/fusion/
#wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
#gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz       
#wget ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz
#gunzip -k Homo_sapiens.GRCh38.87.gtf.gz

cd ~/workspace/inputs/reference/fusion/
ls -halt

# Get annotations for only chromosomes 6 and 17
cat Homo_sapiens.GRCh38.87.gtf | grep --color=never -w "^6" > chr617.gtf
cat Homo_sapiens.GRCh38.87.gtf | grep --color=never -w "^17">> chr617.gtf 
head chr617.gtf
tail chr617.gtf

#### CODE BELOW ALREDY RUN - TAKES >2h
# Parse transcriptome fasta, preserving full ensembl headers
#  Setup directory
#mkdir -p ~/workspace/inputs/reference/fusion/per-feature
#cd ~/workspace/inputs/reference/fusion/per-feature
#  Split fasta at each instance of a sequence header and write to new file
#csplit -s -z ../Homo_sapiens.GRCh38.cdna.all.fa '/>/' '{*}'
#  If from chromosomes 6 or 17, rename files using the columns of the original ensemble header
#####  (This step takes about 120 minutes. You can proceed with the next section in ~/workspace/inputs/data/fastq/chr6_and_chr17)

#Run with Nohup in case connection is lost
#for f in xx*; do awk -F ":" 'NR==1 && $3=="6" || NR==1 && $3=="17"{print $2 "." $3 "." $4 "." $5}' $f | xargs -I{} mv $f {}.fa; done
#  Concatenate features from chromsomes 6 and 17 to a new reference fasta  
#cd ~/workspace/inputs/reference/fusion
#cat ./per-feature/GRCh38.6.*.fa ./per-feature/GRCh38.17.*.fa > chr617.fa
#rm -rf per-feature
```
To get one read pair each for normal and tumor, merge the chr6_and_chr17 only RNA-seq fastqs.

```bash
mkdir -p ~/workspace/inputs/data/fastq/chr6_and_chr17/fusion
cd ~/workspace/inputs/data/fastq/chr6_and_chr17/fusion
cat ../../RNAseq_Norm/RNAseq_Norm_Lane1_R1.fastq.gz ../../RNAseq_Norm/RNAseq_Norm_Lane2_R1.fastq.gz > RNAseq_Norm_R1.fastq.gz
cat ../../RNAseq_Norm/RNAseq_Norm_Lane1_R2.fastq.gz ../../RNAseq_Norm/RNAseq_Norm_Lane2_R2.fastq.gz > RNAseq_Norm_R2.fastq.gz
cat ../../RNAseq_Tumor/RNAseq_Tumor_Lane1_R1.fastq.gz ../../RNAseq_Tumor/RNAseq_Tumor_Lane2_R1.fastq.gz > RNAseq_Tumor_R1.fastq.gz
cat ../../RNAseq_Tumor/RNAseq_Tumor_Lane1_R2.fastq.gz ../../RNAseq_Tumor/RNAseq_Tumor_Lane2_R2.fastq.gz > RNAseq_Tumor_R2.fastq.gz
```

- Subsample fastqs to allow fusion alignment to run quickly (**5 min**):
```bash
# Use seqtk to take subsamples of the 30% of the fastq read pairs  
seqtk sample -s100 RNAseq_Norm_R1.fastq.gz 0.3 > subRNAseq_Norm_R1.fastq.gz
seqtk sample -s100 RNAseq_Norm_R2.fastq.gz 0.3 > subRNAseq_Norm_R2.fastq.gz
seqtk sample -s100 RNAseq_Tumor_R1.fastq.gz 0.3 > subRNAseq_Tumor_R1.fastq.gz
seqtk sample -s100 RNAseq_Tumor_R2.fastq.gz 0.3 > subRNAseq_Tumor_R2.fastq.gz
```

## Run Fusion Alignment 
- Create kallisto index:

```bash
#mkdir -p ~/workspace/rnaseq/fusion
cd ~/workspace/rnaseq/fusion
kallisto index -i index.617.idx -k 31 --make-unique ~/workspace/inputs/reference/fusion/chr617.fa
```
#Call fusions
```bash

cd ~/workspace/rnaseq/fusion
#Quantify potential fusions (**15 min**):
nohup kallisto quant -i index.617.idx --fusion -o kquant-norm617 ~/workspace/inputs/data/fastq/chr6_and_chr17/fusion/subRNAseq_Norm_R1.fastq.gz ~/workspace/inputs/data/fastq/chr6_and_chr17/fusion/subRNAseq_Norm_R2.fastq.gz 2> kallisto_quant_N_bam.out &

nohup kallisto quant -i index.617.idx --fusion -o kquant-tumor617 ~/workspace/inputs/data/fastq/chr6_and_chr17/fusion/subRNAseq_Tumor_R1.fastq.gz ~/workspace/inputs/data/fastq/chr6_and_chr17/fusion/subRNAseq_Tumor_R2.fastq.gz 2> kallisto_quant_T_bam.out &

#Call fusions with pizzly: **120 mins**
cd ~/workspace/rnaseq/fusion
nohup pizzly -k 31 --gtf ~/workspace/inputs/reference/fusion/chr617.gtf --cache index-norm617.cache.txt --align-score 2 --insert-size 400 --fasta ~/workspace/inputs/reference/fusion/chr617.fa --output norm-fuse617 kquant-norm617/fusion.txt 2> pizzly_N_bam.out &

nohup pizzly -k 31 --gtf ~/workspace/inputs/reference/fusion/chr617.gtf --cache index-tumor617.cache.txt --align-score 2 --insert-size 400 --fasta ~/workspace/inputs/reference/fusion/chr617.fa --output tumor-fuse617 kquant-tumor617/fusion.txt 2> pizzly_T_bam.out &

```

- If using 30% of reads with the above process, expect about 13,000 retained transcripts from normal and about 3,000 retained transcripts from tumor. See next section to investigate the output of the pizzly fusion calling.

```bash
#To go the fusion directory
cd ~/workspace/rnaseq/fusion

# ALREADY RUN: Get R scripts for later use
#wget https://raw.githubusercontent.com/griffithlab/pmbio.org/master/assets/course_scripts/mod.grolar.R
#wget https://raw.githubusercontent.com/griffithlab/pmbio.org/master/assets/course_scripts/import_Pizzly.R
scp -ri ~/course_EC2_01.pem ubuntu@ec2-34-203-226-9.compute-1.amazonaws.com:~/workspace/rnaseq/fusion/*.json .
scp -ri ~/course_EC2_01.pem ubuntu@ec2-34-203-226-9.compute-1.amazonaws.com:~/workspace/rnaseq/fusion/mod20220117.grolar.R .

#list files
ls -halt

# Open R on AWS and set working directory 
R
setwd("~/workspace/rnaseq/fusion/")
#Check the path
getwd()

# Load required packages
library(cowplot)
library(jsonlite)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(chimeraviz)

# Change to your own output location
setwd("~")

# Load data
suffix = "617.json"
print(suffix)
JSON_files = list.files(path = "~", pattern = paste0("*",suffix))
print(JSON_files)
Ids = gsub(suffix, "", JSON_files)
print(Ids)

# Assuming you have used GRCh38 gene models
edb <- EnsDb.Hsapiens.v86
listColumns(edb)
supportedFilters(edb)

# Load the functoin https://github.com/MattBashton/grolar/blob/master/grolar.R
source("mod20220117.grolar.R")

# Suffix which apears after sample id in output file name
# Use above funciton on all output files
lapply(Ids, function(x) GetFusionz_and_namez(x, suffix = "617.json") )

# This function will flatten out the JSON giving you a list of gene A and gene B sorted by splitcount then paircount
# Each fusion now has a unique identifier, sequence positions, and distance values for genes from the same chromosome:

# Read the flattened files
normal=read.table("~/norm-fuse_fusions_filt_sorted.txt", header=T)
tumor=read.table("~/tumor-fuse_fusions_filt_sorted.txt", header=T)

#Investigate the first five rows
head(normal, 5)
head(tumor, 5)

# Common filtering tasks for fusion output include removing fusions from the tumor sample which are present in the normal, and removing fusions for which the is little support by pair and split read counts:

normal$genepr=paste0(normal$geneA.name,".",normal$geneB.name)
tumor$genepr=paste0(tumor$geneA.name,".",tumor$geneB.name)
uniqueTumor=subset(tumor, !(tumor$genepr %in% normal$genepr))
nrow(uniqueTumor)==nrow(tumor)
[1] FALSE
nrow(tumor)-nrow(uniqueTumor)
[1] 2
# There are two fusions (or at least fusion gene pairs) from the normal sample which are also present in the tumor. 
# Examine the output of- 
shared_src_tumor=subset(tumor, (tumor$genepr %in% normal$genepr))
shared_src_normal=subset(normal, (normal$genepr %in% tumor$genepr))
shared_src_tumor
shared_src_normal

# Filtering by counts:

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

nrow(allfusions)
allfusions=allfusions[which(allfusions$paircount >= 2 & allfusions$splitcount >= 5),]
nrow(allfusions)
#write.table(allfusions, "allfusions.txt")
write.table(allfusions, file="allfusions.txt",row.names=F, col.names=T, quote=F, sep="\t")

#Chimeraviz is an R package for visualizing fusions from RNA-seq data. The chimeraviz package has import functions built in for a variety of fusion-finder programs, but not for pizzly. We will have to load our own import function that you downloaded above:

# Enter R, install and load chimeraviz 
R
source("https://bioconductor.org/biocLite.R")
biocLite("chimeraviz")
# (if asked to update old packages, you can ignore- Update all/some/none? [a/s/n]:)
library(chimeraviz)

# Use the pizzly importer script to import fusion data
source("./import_Pizzly.R")
#  You can view the function by calling it without variables
#importPizzly
fusions = importPizzly("./allfusions.txt","hg38")
plot_circle(fusions)
```