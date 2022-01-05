
# HUMAN GENOME REFERENCE FILES

```bash
# Make sure CHRS environment variable is set. 
echo $CHRS

# Create a directory for reference genome files and enter this dir.
#Command below already run
#mkdir -p ~/workspace/inputs/references/genome
cd ~/workspace/inputs/references/genome

# Dowload human reference genome files.
#Command below already run
#wget http://genomedata.org/pmbio-workshop/references/genome/$CHRS/ref_genome.tar

# Unpack the archive using `tar -xvf` (`x` for extract, `v` for verbose,
#  `f` for file).
tar -xvf ref_genome.tar

# View contents.
tree

# Remove the archive.
rm -f ref_genome.tar

# Uncompress the reference genome FASTA file.
gunzip ref_genome.fa.gz

# View contents.
tree

# Check the chromosome headers in the fasta file.
cat ref_genome.fa | grep -P "^>"
```

## Split the long fasta by chromosome

```bash
# Make new directory and change directories.
mkdir -p ~/workspace/inputs/references/genome/ref_genome_split/
cd ~/workspace/inputs/references/genome

# Split.
faSplit byname ref_genome.fa ./ref_genome_split/
```

## Explore the contents of the reference genome file

```bash
# View the first 10 lines of this file. Note the header line starting with `>`. 
# Why does the sequence look like this?
cd ~/workspace/inputs/references/genome
head -n 10 ref_genome.fa

# Pull out only the header lines.
grep ">" ref_genome.fa

# How many lines and characters are in this file?
wc ref_genome.fa

# How long are to two chromosomes combined (in bases and Mbp)? Use grep to skip the header lines for each chromosome.
grep -v ">" ref_genome.fa | wc

# How long does that command take to run?
time grep -v ">" ref_genome.fa | wc

# View 10 lines from approximately the middle of this file
head -n 2500000 ref_genome.fa | tail

# What is the count of each base in the entire reference genome file (skipping the header lines for each sequence)?
# Runtime: ~30s
cat ref_genome.fa | grep -v ">" | perl -ne 'chomp $_; $bases{$_}++ for split //; if (eof){print "$_ $bases{$_}\n" for sort keys %bases}'
# What does each of these bases refer to? What are the "unexpected bases"?
```

## Index the fasta files

```bash
# first remove the .fai and .dict files that were downloaded. Do not remove the .fa file though!
cd ~/workspace/inputs/references/genome
rm -f ref_genome.fa.fai ref_genome.dict

# Use samtools to create a fasta index file.
samtools faidx ref_genome.fa

# View the contents of the index file.
head ref_genome.fa.fai

# Use picard to create a dictionary file.
java -jar $PICARD CreateSequenceDictionary -R ref_genome.fa -O ref_genome.dict

# View the content of the dictionary file.
cat ref_genome.dict

#Also index the split chromosomes.
samtools faidx ./ref_genome_split/chr6.fa
samtools faidx ./ref_genome_split/chr17.fa

# Create reference index for the genome to use BWA
bwa index ref_genome.fa
```

# TRANSCRIPTOME REFERENCE FILES

```bash
# Make sure CHRS environment variable is set.
echo $CHRS

# Create a directory for transcriptome input files.
#Command below already run
#mkdir -p ~/workspace/inputs/references/transcriptome
cd ~/workspace/inputs/references/transcriptome

# Download the files.
#Command below already run
#wget http://genomedata.org/pmbio-workshop/references/transcriptome/$CHRS/ref_transcriptome.gtf
#wget http://genomedata.org/pmbio-workshop/references/transcriptome/$CHRS/ref_transcriptome.fa

# Take a look at the contents of the gtf file. Press 'q' to exit the 'less' display.
less -p start_codon -S ref_transcriptome.gtf
```

## Explore the contents of the transcriptome reference files

```bash
#How many chromsomes are represented?
cut -f1 ref_transcriptome.gtf | sort | uniq -c

# How many unique gene IDs are in the .gtf file?
# We can use a perl command-line command to find out:
perl -ne 'if ($_ =~ /(gene_id\s\"ENSG\w+\")/){print "$1\n"}' ref_transcriptome.gtf | sort | uniq | wc -l

# what are all the feature types listed in the third column of the GTF?
# how does the following command (3 commands piped together) answer that question?
cut -f 3 ref_transcriptome.gtf | sort | uniq -c
```
## Create a reference index for transcriptome with HISAT for splice RNA alignments to the genome

```bash
cd ~/workspace/inputs/references/transcriptome

# Create a database of observed splice sites represented in our reference transcriptome GTF
~/workspace/bin/hisat2-2.1.0/hisat2_extract_splice_sites.py ref_transcriptome.gtf > splicesites.tsv
head splicesites.tsv

# Create a database of exon regions in our reference transcriptome GTF
~/workspace/bin/hisat2-2.1.0/hisat2_extract_exons.py ref_transcriptome.gtf > exons.tsv
head exons.tsv

# build the reference genome index for HISAT and supply the exon and splice site information extracted in the previous steps
# specify to use 8 threads with the `-p 8` option
# run time for this index is ~5 minutes
~/workspace/bin/hisat2-2.1.0/hisat2-build -p 8 --ss splicesites.tsv --exon exons.tsv ~/workspace/inputs/references/genome/ref_genome.fa ref_genome
```

## Create a reference transcriptome index for use with Kallisto

```bash
cd ~/workspace/inputs/references/transcriptome
mkdir kallisto
cd kallisto

# tidy up the headers to just include the ensembl transcript ids
cat ../ref_transcriptome.fa | perl -ne 'if ($_ =~ /\d+\s+(ENST\d+)/){print ">$1\n"}else{print $_}' > ref_transcriptome_clean.fa

# run time for this index is ~30 seconds
kallisto index --index=ref_transcriptome_kallisto_index ref_transcriptome_clean.fa

```
