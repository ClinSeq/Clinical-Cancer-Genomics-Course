
# Human genome reference files

The human genome is work in progress. Different versions exist (assemblies) and is often a source of confusion in genomics. Quick use of published data is a bit more challenging if a differnt genome build was applied than what is used in your lab. 

Also, there may even be differences in the *latest* versions, e.g. the GRCh37 and the hg19 assemblies from NCBI and USCS, respectively, had different mitochondrial genome. 

In genomics, the reference genome offers a scaffold upon which new data can be mapped, which is a much more efficient way rather than building a genome from scratch.

Avery nice overview of the human genome versions can be found [here](http://genome.ucsc.edu/FAQ/FAQreleases.html).

During this course version GRCh38, from the Genome Reference Consortium, with modifications fromt he 1000 genomes consortium will be used. It includes extra decoy and HLA sequences. The files are already available on AWS but can be downloaded from [here](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/).

As several of the analysis steps wold take many hours we will use a smaller version of the reference genome which contains chr6 and chr17. These were choosen because the datat that will be analysed contain interesting varaiants on these chromosomes. 
- Connect to your AWS instance


```bash
# Remember, source the settings, paths and variables file
source .bashrc

# Check the CHRS variable.
echo $CHRS

# This variable should return "chr6_and_chr17", if not, reach of for assistance.

# Go to the directory for the reference genome files
cd ~/workspace/inputs/references/genome

# List the files
ls -halt

# Check what the reference genome looks like. When you have the reference genome open just type "1000" and press enter. You will then jump 1000 lines.
less -SN ref_genome.fa

#less -SN will be used during the course. Check the man page. Give it some time to understand how it it structured, it is very helpful to be able to read the help-pages fast when working on the command line.
man less

# Check the chromosome headers in the fasta file.
cat ref_genome.fa | grep -P "^>"

# View the first 10 lines of this file. Note the header line starting with `>`. 
# Why does the sequence look like this?
cd ~/workspace/inputs/references/genome
head -n 10 ref_genome.fa

# How many lines and characters are in this file?
wc ref_genome.fa

# How long are to two chromosomes combined (in bases and Mbp)? Use grep to skip the header lines for each chromosome.
grep -v ">" ref_genome.fa | wc

# How long does that command take to run?
time grep -v ">" ref_genome.fa | wc

# View 10 lines from the end of the filefile
tail -n 10 ref_genome.fa

# What is the count of each base in the entire reference genome file (skipping the header lines for each sequence)?
# Runtime: ~30s
cat ref_genome.fa | grep -v ">" | perl -ne 'chomp $_; $bases{$_}++ for split //; if (eof){print "$_ $bases{$_}\n" for sort keys %bases}'
# What does each of these bases refer to? What are the "unexpected bases"?
# Google "IUPAC nucleic acid codes"
```

## Split the long fasta by chromosome
```bash
cd ~/workspace/inputs/references/genome

# Split.
faSplit byname ref_genome.fa ./ref_genome_split/

# View contents.
tree
```

## Index the fasta files

Indexing is widely used in bioinformatics workflows to improve performance. Typically it is applied to large data files with many records to improve the ability of tools to rapidly access specific locations of the file. For example, if we have alignments against the entire genome and we want to visualize those alignments for a single gene on chr17 in a genome viewer such as IGV, we don’t want the viewer to have to scan through the entire file. Indexing allows us to jump right to the correct place in the file and pull out just the information we need without reading much of the file.

```bash
# first remove the .fai and .dict files that were downloaded. Do not remove the .fa file though!
cd ~/workspace/inputs/references/genome
ls -halt
rm -f ref_genome.fa.fai ref_genome.dict

# Use samtools to create a fasta index file.
samtools faidx ref_genome.fa

# View the contents of the index file.
head ref_genome.fa.fai
# What information does the index fail store?
# google "fasta index fai format" or something similar

# Use picard to create a dictionary file.
java -jar $PICARD CreateSequenceDictionary -R ref_genome.fa -O ref_genome.dict

# View the content of the dictionary file.
head ref_genome.dict
# less can also be applied

#Also index the split chromosomes.
samtools faidx ./ref_genome_split/chr6.fa
samtools faidx ./ref_genome_split/chr17.fa

# Create reference index for the genome to use BWA
bwa index ref_genome.fa
```

# Human genome transcriptome reference files
```bash
# Make sure CHRS environment variable is set.
echo $CHRS
cd ~/workspace/inputs/references/transcriptome

# List files in the directory
ls -halt

# Take a look at the contents of the gtf file. Press 'q' to exit the 'less' display. 
# When the file is displayed usin "less" type "/" and then a something you want to highlight e.g. "chr6". This is very useful when searching for specific infomration
less -SN ref_transcriptome.gtf
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

Below we will prepare files for the bioinformatic tools (HISAT/Kallisto) used to analyse RNA-sequencing data.

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
