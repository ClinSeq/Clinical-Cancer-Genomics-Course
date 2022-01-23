# Annotation files

## Files for annotation of variation in the human genome
```bash

#These are available from the Broad institute in Boston, who also provides the GATK software suite.
cd ~/workspace/inputs/references/gatk

#Download the following files:
# SNP calibration call sets - dbsnp, hapmap, omni, and 1000G
# Runtime: < 2min
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf .
# Runtime: ~ 2min
bgzip --threads 8 Homo_sapiens_assembly38.dbsnp138.vcf
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz .
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz .
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz .

# Have a loot what the files contain. 
# Dont worry, we will look at VCF files later.
gunzip -c 1000G_phase1.snps.high_confidence.hg38.vcf.gz | less -SN

# Indel calibration call sets - dbsnp, Mills
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz .
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz .

# Interval lists that can be used to parallelize certain GATK tasks
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list .
gsutil cp -r gs://genomics-public-data/resources/broad/hg38/v0/scattered_calling_intervals/ .

# list the files we just downloaded
ls -lh

```

## Index the variation files

```bash
cd ~/workspace/inputs/references/gatk/

#SNP calibration call sets - dbsnp, hapmap, omni, and 1000G
# Runtime: ~ 4min
gatk --java-options '-Xmx12g' IndexFeatureFile -I ~/workspace/inputs/references/gatk/Homo_sapiens_assembly38.dbsnp138.vcf.gz
gatk --java-options '-Xmx12g' IndexFeatureFile -I ~/workspace/inputs/references/gatk/hapmap_3.3.hg38.vcf.gz
gatk --java-options '-Xmx12g' IndexFeatureFile -I ~/workspace/inputs/references/gatk/1000G_omni2.5.hg38.vcf.gz
# Runtime: ~ 3min
gatk --java-options '-Xmx12g' IndexFeatureFile -I ~/workspace/inputs/references/gatk/1000G_phase1.snps.high_confidence.hg38.vcf.gz

#Indel calibration call sets - dbsnp, Mills
gatk --java-options '-Xmx12g' IndexFeatureFile -I ~/workspace/inputs/references/gatk/Homo_sapiens_assembly38.known_indels.vcf.gz
gatk --java-options '-Xmx12g' IndexFeatureFile -I ~/workspace/inputs/references/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

```

## Interval files and coordinates for the exome sequencing assay
These files contain the location of the exons in the human genome that were "targeted", that is for which targeting oligos were created (primary_targets). 

The files also contain the location of the actual oligos (capture_targets).

```bash
# change directories
cd ~/workspace/inputs/references/exome

# download the files
wget -c https://sequencing.roche.com/content/dam/rochesequence/worldwide/shared-designs/SeqCapEZ_Exome_v3.0_Design_Annotation_files.zip
unzip SeqCapEZ_Exome_v3.0_Design_Annotation_files.zip

# remove the zip
rm -f SeqCapEZ_Exome_v3.0_Design_Annotation_files.zip

# Lift-over of the Roche coordinates from hg19 to the hg38 assembly.
# the software is availble from USCS, downloadble here:
#wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
#chmod +x liftOver
# the chain file
wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# run liftover
liftOver SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_primary_targets.bed  hg19ToHg38.over.chain.gz SeqCap_EZ_Exome_v3_hg38_primary_targets.bed unMapped.bed

liftOver SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_capture_targets.bed  hg19ToHg38.over.chain.gz SeqCap_EZ_Exome_v3_hg38_capture_targets.bed unMapped1.bed

# create a version in standard bed format (chr, start, stop)
cut -f 1-3 SeqCap_EZ_Exome_v3_hg38_primary_targets.bed > SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed

cut -f 1-3 SeqCap_EZ_Exome_v3_hg38_capture_targets.bed > SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.bed

# take a quick look at the format of these files
head SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed
head SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.bed
```

## Calculate the size of the SeqCap v3 exome
```bash
#This can be done in many ways - give it a try yourself before trying the code below and compare results

# first sort the bed files and store the sorted versions
bedtools sort -i SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed > SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.bed
bedtools sort -i SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.bed > SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.sort.bed

# now merge the bed files to collapse any overlapping regions so they are not double counted.
bedtools merge -i SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.bed > SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.merge.bed

bedtools merge -i SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.sort.bed > SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.sort.merge.bed

# finally use a Perl one liner to determine the size of the files in Mb
FILES=(SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.sort.merge.bed SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.merge.bed)
echo ${FILES[0]}

for FILE in ${FILES[@]}
do
echo "--------------------------------------------------------"
echo $FILE
#With merge
cat $FILE | sortBed -i stdin | mergeBed -i stdin | perl -e '$counter=0; while(<>){chomp; @array = split "\t", $_; $counter = $counter + ($array[2] - $array[1]);}; print $counter/1000000, "\n";'
done
# note that the size of the space targeted by the exome reagent is ~63 Mbp. Does that sound reasonable?

# now create a subset of these bed files
grep -w -P "^chr6|^chr17" SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.merge.bed > exome_regions.bed
#When creating files, make a habit to investigate the output to avoid downstream confusion
head -n 10 exome_regions.bed

grep -w -P "^chr6|^chr17" SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.sort.merge.bed > probe_regions.bed
head -n 10 probe_regions.bed

# clean up intermediate files
#rm -fr SeqCapEZ_Exome_v3.0_Design_Annotation_files/ SeqCap_EZ_Exome_v3_hg38_primary_targets.bed SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.bed SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.bed unMapped.bed
```

## Create an inverval list for the exome bed files

```bash
# first for the complete exome and probe bed file
cd ~/workspace/inputs/references/
#Commands below already run
#mkdir temp
cd temp
#Commands below already run
#wget http://genomedata.org/pmbio-workshop/references/genome/all/ref_genome.dict
cd ~/workspace/inputs/references/exome
java -jar $PICARD BedToIntervalList I=SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.merge.bed O=SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.merge.interval_list SD=~/workspace/inputs/references/temp/ref_genome.dict
java -jar $PICARD BedToIntervalList I=SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.sort.merge.bed O=SeqCap_EZ_Exome_v3_hg38_capture_targets.v2.sort.merge.interval_list SD=~/workspace/inputs/references/temp/ref_genome.dict
#rm -fr ~/workspace/inputs/references/temp/
#Explore the interval-lists and what the files contain, scroll down until the coordinates are shown
less -SN SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.sort.merge.interval_list

# next for our subset exome and probe regions file
cd ~/workspace/inputs/references/exome
java -jar /usr/local/bin/picard.jar BedToIntervalList I=exome_regions.bed O=exome_regions.bed.interval_list SD=~/workspace/inputs/references/genome/ref_genome.dict
java -jar /usr/local/bin/picard.jar BedToIntervalList I=probe_regions.bed O=probe_regions.bed.interval_list SD=~/workspace/inputs/references/genome/ref_genome.dict
```
