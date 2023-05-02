# Calling germline- and somatic mutations

## Germline variants
We will start by calling germline variants. For this purpose we will use the [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/4414586765723-HaplotypeCaller) from the Broad Institute. We will use this tool on the germline DNA bam files that were processed yesterday. 

Active forum and tutorials for the GATK software suit is available [here](https://gatk.broadinstitute.org/hc/en-us/community/topics).

The HaplotypeCaller can be run in single- or joint sample mode. The joint sample mode gives improved accuracy, especially in samples with low coverage. For simplicity, we will only use the single sample mode in this course. 

When haplotype caller detects variation it performes local realignment of the reads which improves variant accuracy, especially in e.g. repetitive regions. 

```bash
#Make sure that $GATK_REGIONS is set correctly
echo $GATK_REGIONS

#Create working dir for germline results if not already created
mkdir -p ~/workspace/germline
cd ~/workspace/germline

#The argument --java-options '-Xmx12g' tells GATK to use 12GB of memory
#To list the haplotype caller arguments run 
gatk --java-options '-Xmx12g' HaplotypeCaller --help

#Call variants for exome data
gatk --java-options '-Xmx12g' HaplotypeCaller \
    -R ~/workspace/inputs/references/genome/ref_genome.fa \
    -I ~/workspace/align/Exome_Norm_sorted_mrkdup_bqsr.bam \
    -O ~/workspace/germline/Exome_Norm_HC_calls.vcf \
    --bam-output ~/workspace/germline/Exome_Norm_HC_out.bam $GATK_REGIONS
```

## Explore the ouput files

### Realigned bam file

Navigate to the folder where you want to place the BAMs
e.g. ~/course_bams

```bash
#Copy the haplotype caller realigned bam-file to your local computer 
scp -ri ~/PEM_KEY_ID.pem ubuntu@AWS_ADDRESS_HERE:~/workspace/germline/Exome_Norm_HC_out.ba* .
```

Launch IGV (remember **hg38**)
- Open the bam files
    - Exome_Norm_sorted_mrkdup_bqsr.bam (intput file to haplotype caller).
    - Exome_Norm_HC_out.bam (realigned bam file by haplotype caller).
- Navigate to chr17:76286974-76287203.
    - Color alignments by: no color.
    - Enable ruler in the menue bar.
    - In the first base of the indel: Sort aligments according to base.
    - Select squished view.
    - Zoom out.
- Try to answer the following:
    - What characterises misaligned variant-supporting reads in non-realigned bam file?
    - Are there any false positive variants take in the non-realigned bam file?
- ***OBS***: Leave IGV open, we will use it in the below section.

![](https://i.imgur.com/yQS7rie.png)


### VCF (variant call format) file

```bash
#Make sure you hare in the same directory as the output files
cd ~/workspace/germline

#List files
ls -halt

#Open the VCF file
less -SN Exome_Norm_HC_calls.vcf
```

The VCF file format is quite complex, the specification can be found [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

Try to do the following:

- Identify the VCF header section.
- Identify the data column header.
- Find the variant on chr17 pos 106742 in the VCF file.
    - Lead, type /106742 and press ENTER
    - What is the read depth in this position?
    - Is it a high quality variant?
    - Go go IGV in this position
        - Select expanded view.
        - Color alignments by read strand.
        - Sort alignments by base.
        - Does the visual impression in IGV support the quality of the variant?
- Identify a way to quickly separate SNPs and INDELs in the VCF.
- Find the variant on chr17 pos 1787008 in the VCF file
    - Is this a high quality variant?
    - Navgiate in IGV to the position.
    - Try to explain the information in the VCF to what is shown in IGV, is it correct?

### Filtering germline variants

As noted above, the raw output of GATK HaplotypeCaller will include many variants with varying degrees of quality. For various reasons we might wish to further filter these to a higher confidence set of variants. The hard-filtering approached applied here is futher described in detail at the [GATK webside](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants). 

### Separate SNPs from Indels in the variant call set

First, we will separate out the SNPs and Indels from the VCF into new separate VCFs. Note that the variant type (SNP, INDEL, MIXED, etc) is not stored explicitly in the vcf but instead inferred from the genotypes. We will use a versatile GATK tool called `SelectVariants`. This command can be used for all kinds of simple filtering or subsetting purposes. We will run it twice to select by variant type, once for SNPs, and then again for Indels, to produce two new VCFs.

```bash
cd ~/workspace/germline/
gatk --java-options '-Xmx12g' SelectVariants \
    -R ~/workspace/inputs/references/genome/ref_genome.fa \
    -V ~/workspace/germline/Exome_Norm_HC_calls.vcf -select-type SNP \
    -O ~/workspace/germline/Exome_Norm_HC_calls.snps.vcf
gatk --java-options '-Xmx12g' SelectVariants \
    -R ~/workspace/inputs/references/genome/ref_genome.fa \
    -V ~/workspace/germline/Exome_Norm_HC_calls.vcf -select-type INDEL \
    -O ~/workspace/germline/Exome_Norm_HC_calls.indels.vcf
```
### Apply filters to the SNP and Indel call sets

Next, we will perform so-called hard-filtering by applying a number of cutoffs. Multiple filters can be combined arbitrarily. Each can be given its own name so that you can later determine which one or more filters a variant fails. Visit the [GATK documentation on hard-filtering](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants) to learn more about the following hard filtering options. Notice that different filters and cutoffs are recommended for SNVs and Indels. This is why we first split them into separate files.

```bash
cd ~/workspace/germline/
gatk --java-options '-Xmx12g' VariantFiltration \
    -R ~/workspace/inputs/references/genome/ref_genome.fa \
    -V ~/workspace/germline/Exome_Norm_HC_calls.snps.vcf \
    --filter-expression "QD < 2.0" --filter-name "QD_lt_2" \
    --filter-expression "FS > 60.0" --filter-name "FS_gt_60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ_lt_40" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRS_lt_n12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS_lt_n8" \
    --filter-expression "SOR > 3.0" --filter-name "SOR_gt_3" \
    -O ~/workspace/germline/Exome_Norm_HC_calls.snps.filtered.vcf

gatk --java-options '-Xmx12g' VariantFiltration \
    -R ~/workspace/inputs/references/genome/ref_genome.fa \
    -V ~/workspace/germline/Exome_Norm_HC_calls.indels.vcf \
    --filter-expression "QD < 2.0" --filter-name "QD_lt_2" \
    --filter-expression "FS > 200.0" --filter-name "FS_gt_200" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "RPRS_lt_n20" \
    --filter-expression "SOR > 10.0" --filter-name "SOR_gt_10" \
    -O ~/workspace/germline/Exome_Norm_HC_calls.indels.filtered.vcf
```

Notice that warnings are given with regards to MQRankSum and ReadPosRankSum, these are only calculated if the site is called as heterozygous. E.g. for this site (inspect in IGV and in the VCF file):

- chr17:1067361

This was discussed at the [GATK forum](https://sites.google.com/a/broadinstitute.org/legacy-gatk-forum-discussions/2017-01-18-2016-08-11/8323-Why-do-MQRankSum-and-ReadPosRankSum-not-appear-in-some-vcf-file-entries).

The variants are now marked as PASS or or *filter-name*. Use `grep -v` and `head` to skip past all the VCF header lines and view the first few records. 

```bash
#inspect the 10 first variants that did not pass
grep -v "##" Exome_Norm_HC_calls.snps.filtered.vcf | grep -vP "\tPASS\t" | head -10
#inspect the 10 first variants that did pass
grep -v "##" Exome_Norm_HC_calls.snps.filtered.vcf | grep -P "\tPASS\t" | head -10
```
Try to do the following:

- Count the number of variants that failed and passed filtering.

#### Merge filtered SNP and INDEL vcfs back together

```bash
gatk --java-options '-Xmx12g' MergeVcfs \
    -I ~/workspace/germline/Exome_Norm_HC_calls.snps.filtered.vcf \
    -I ~/workspace/germline/Exome_Norm_HC_calls.indels.filtered.vcf \
    -O ~/workspace/germline/Exome_Norm_HC_calls.filtered.vcf

#Nbr of variants in the SNP file
grep -v "##" Exome_Norm_HC_calls.snps.filtered.vcf | wc -l

#Nbr of variants in the INDEL file
grep -v "##" Exome_Norm_HC_calls.indels.filtered.vcf | wc -l

#Nbr of variants in the merged file
grep -v "##" Exome_Norm_HC_calls.filtered.vcf | wc -l
```

#### Extract PASS variants only

It would also be convenient to have a vcf with only passing variants. For this, we can go back the `GATK SelectVariants` tool. This will be run much as above except with the `--exlude-filtered` option instead of `-select-type`.

```bash
gatk --java-options '-Xmx12g' SelectVariants \
    -R ~/workspace/inputs/references/genome/ref_genome.fa \
    -V ~/workspace/germline/Exome_Norm_HC_calls.filtered.vcf \
    -O ~/workspace/germline/Exome_Norm_HC_calls.filtered.PASS.vcf \
    --exclude-filtered
```

### Perform annotation of filtered variants
Now that we have high-confidence, filtered variants, we want to start understanding which of these variants might be clinically or biologically relevant. Ensembl's Variant Effect Predictor (VEP) annotation software is a powerful tool for annotating variants with a great deal of biological features. This includes such information as protein consequence (non-coding or coding), population frequencies, links to external databases, various scores designed to estimate the importance of individual variants on protein function, and much more.

Note, the first time you run VEP it will create a fasta index for the reference genome in VEP cache. Therefore, it will take longer to run that first time but should speed up substantially for subsequent runs on files with similar numbers of variants.


```bash
#VEP annotate hard-filtered exome results

####### Already done code below
#cd ~/workspace/vep_cache/
#wget https://github.com/Ensembl/VEP_plugins/archive/refs/heads/release/104.zip
#unzip 104.zip 
#ls -halt ~/workspace/vep_cache/
#Already done
#mv VEP_plugins-release-104 Plugins
#ls -halt ~/workspace/vep_cache/Plugins

#get annotation fasta for homo_sapiens 
#cd ~/workspace/vep_cache/homo_sapiens
#Already done
#wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

####### Start here
cd ~/workspace/germline
#Remember, we are using the REF ~/workspace/inputs/references/genome/ref_genome.fa which only contains chr6 and chr17.
nohup vep --cache --dir_cache ~/workspace/vep_cache \
    --dir_plugins ~/workspace/vep_cache/Plugins \
    --fasta ~/workspace/inputs/references/genome/ref_genome.fa \
    --fork 8 --assembly=GRCh38 --offline --vcf \
    -i ~/workspace/germline/Exome_Norm_HC_calls.filtered.PASS.vcf \
    -o ~/workspace/germline/Exome_Norm_HC_calls.filtered.PASS.vep.vcf \
    --check_existing --total_length --allele_number --no_escape --everything \
    --use_given_ref --force_overwrite --coding_only &

#Open the VEP vcf, row 39 and 40, contains information about the VEP annotation.
less -SN Exome_Norm_HC_calls.filtered.PASS.vep.vcf
#Or just grep the VEP info.
grep "##" Exome_Norm_HC_calls.filtered.PASS.vep.vcf | grep "VEP" | less -SN
#Copy the html VEP summary to your local computer
scp -ri ~/PEM_KEY_ID.pem ubuntu@AWS_ADDRESS_HERE:~/workspace/germline/Exome_Norm_HC_calls.filtered.PASS.vep.vcf_summary.html .
```

Open the HTML file.

***Answer the following:***

- How many variants were processed by VEP?
- What is the most common kind of consequence?
- Why are only variants the coding regions detected?
- Open the Exome_Norm_HC_out.bam and Exome_Norm_sorted_mrkdup_bqsr.bam in IGV 
    - Find a missense variant and examine it in IGV, does the variant look real?
    - E.g. take chr17 744946
        - Find the variant in the VEP VCF as well as IGV.
- VEP assigns each variant a consequence types, as described [here](https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html).
    - What is the difference between splice region variants and splice donor/acceptor variants?
    - The --pick option is used, as described [here](https://m.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick), do you think it makes sense to use? 
----
## Somatic variants

We will run multiple variant callers and merging the callset.

### Varscan

First off is [VARSCAN](http://varscan.sourceforge.net/), that employs a robust heuristic/statistic approach to call variants that meet desired thresholds for read depth, base quality, variant allele frequency, and statistical significance.

As seen below varscan uses the mpileup command from samtools. What the mpileup command does can be explored [here](https://www.htslib.org/doc/samtools-mpileup.html)

```bash
mkdir -p ~/workspace/somatic/varscan
cd ~/workspace/somatic/varscan

#Have a look at the input data to the varscan caller
samtools mpileup -l ~/workspace/inputs/references/exome/exome_regions.bed \
    --no-BAQ -f ~/workspace/inputs/references/genome/ref_genome.fa \
    ~/workspace/align/Exome_Norm_sorted_mrkdup_bqsr.bam \
    ~/workspace/align/Exome_Tumor_sorted_mrkdup_bqsr.bam | less -SN

# Run varscan
java -Xmx12g -jar ~/workspace/bin/VarScan.v2.4.2.jar somatic \
    <(samtools mpileup -l ~/workspace/inputs/references/exome/exome_regions.bed \
        --no-BAQ -f ~/workspace/inputs/references/genome/ref_genome.fa \
        ~/workspace/align/Exome_Norm_sorted_mrkdup_bqsr.bam \
        ~/workspace/align/Exome_Tumor_sorted_mrkdup_bqsr.bam) \
    ~/workspace/somatic/varscan/exome --min-var-freq 0.05 --mpileup 1 --output-vcf

ls -halt

java -Xmx12g -jar ~/workspace/bin/VarScan.v2.4.2.jar processSomatic \
    exome.snp.vcf exome.snp --min-tumor-freq 0.05 --max-normal-freq 0.01

java -Xmx12g -jar ~/workspace/bin/VarScan.v2.4.2.jar processSomatic \
    exome.indel.vcf exome.indel --min-tumor-freq 0.05 --max-normal-freq 0.01

find ~/workspace/somatic/varscan -name '*.vcf' -exec bgzip -f {} \;
find ~/workspace/somatic/varscan -name '*.vcf.gz' -exec tabix -f {} \;

gatk VariantFiltration -R ~/workspace/inputs/references/genome/ref_genome.fa \
    -V exome.snp.Somatic.vcf.gz --mask exome.snp.Somatic.hc.vcf.gz \
    --mask-name "processSomatic" --filter-not-in-mask -O exome.snp.Somatic.hc.filter.vcf.gz

gatk VariantFiltration -R ~/workspace/inputs/references/genome/ref_genome.fa \
    -V exome.indel.Somatic.vcf.gz --mask exome.indel.Somatic.hc.vcf.gz \
    --mask-name "processSomatic" --filter-not-in-mask -O exome.indel.Somatic.hc.filter.vcf.gz

bcftools concat -a -o exome.vcf.gz -O z exome.snp.Somatic.hc.filter.vcf.gz exome.indel.Somatic.hc.filter.vcf.gz
tabix -f ~/workspace/somatic/varscan/exome.vcf.gz
```

### Strelka

Now we are going to run the second variant caller, [STRELKA](https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md). Strelka calls germline and somatic small variants from mapped sequencing reads and is optimized for rapid clinical analysis of germline variation in small cohorts and somatic variation in tumor/normal sample pairs. Both germline and somatic callers include a final empirical variant rescoring step using a random forest model to reflect numerous features indicative of call reliability which may not be represented in the core variant calling probability model.

```bash
mkdir -p ~/workspace/somatic/strelka/exome
cd ~

source activate strelka-env
~/workspace/bin/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam=workspace/align/Exome_Norm_sorted_mrkdup_bqsr.bam \
    --tumorBam=workspace/align/Exome_Tumor_sorted_mrkdup_bqsr.bam \
    --referenceFasta=workspace/inputs/references/genome/ref_genome.fa \
    --exome --runDir=workspace/somatic/strelka/exome

#Please specify according to the number of cpus available or how many you would like to allocate to this job. In this case, four were given.
# Runtime: ~ 3min
python2 ~/workspace/somatic/strelka/exome/runWorkflow.py -m local -j 8

conda deactivate

cd ~/workspace/somatic/strelka/exome/results/variants
zcat somatic.snvs.vcf.gz | \
    awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - \
    > somatic.snvs.gt.vcf
zcat somatic.indels.vcf.gz | \
    awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - \
    > somatic.indels.gt.vcf
find ~/workspace/somatic/strelka/exome/results/variants/ -name "*.vcf" -exec bgzip -f {} \;
find ~/workspace/somatic/strelka/exome/results/variants/ -name "*.vcf.gz" -exec tabix -f {} \;

bcftools concat -a -o exome.vcf.gz -O z somatic.snvs.gt.vcf.gz somatic.indels.gt.vcf.gz

tabix exome.vcf.gz
```

### Mutect2

Last up is [MuTect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php). MuTect2 is a somatic SNP and indel caller that combines the DREAM challenge-winning somatic genotyping engine of the original MuTect (Cibulskis et al., 2013) with the assembly-based machinery of HaplotypeCaller.

### Exome data commands
```bash
#Obtaining germline resource from GATK
cd ~/workspace/inputs/references
#Already done
#gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz .
#gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi .

mkdir -p ~/workspace/somatic/mutect
cd ~/workspace/somatic/mutect

#Creating a panel of normals
# Runtime: ~ 17min
nohup gatk --java-options "-Xmx12G" Mutect2 \
    -R ~/workspace/inputs/references/genome/ref_genome.fa \
    -I ~/workspace/align/Exome_Norm_sorted_mrkdup_bqsr.bam \
    -tumor-sample HCC1395BL_DNA -O Exome_Norm_PON.vcf.gz &

##### WAIT FOR THIS NOHUP STEP TO FINISH BEFORE CONTINUING #####
##### WAIT FOR THIS NOHUP STEP TO FINISH BEFORE CONTINUING #####
##### WAIT FOR THIS NOHUP STEP TO FINISH BEFORE CONTINUING #####

#Running Mutect2 Using latest version of GATK
# Runtime: ~20m
nohup gatk --java-options "-Xmx12G" Mutect2 \
    -R ~/workspace/inputs/references/genome/ref_genome.fa \
    -I ~/workspace/align/Exome_Tumor_sorted_mrkdup_bqsr.bam \
    -tumor HCC1395_DNA -I ~/workspace/align/Exome_Norm_sorted_mrkdup_bqsr.bam \
    -normal HCC1395BL_DNA --germline-resource ~/workspace/inputs/references/af-only-gnomad.hg38.vcf.gz \
    --af-of-alleles-not-in-resource 0.00003125 --panel-of-normals \
    ~/workspace/somatic/mutect/Exome_Norm_PON.vcf.gz -O ~/workspace/somatic/mutect/exome.vcf.gz \
    -L chr6 -L chr17 &

##### WAIT FOR NOHUP FOR THIS STEP TO FINISH BEFORE CONTINUING #####
##### WAIT FOR NOHUP FOR THIS STEP TO FINISH BEFORE CONTINUING #####
##### WAIT FOR NOHUP FOR THIS STEP TO FINISH BEFORE CONTINUING #####


# Filtering mutect variants
gatk --java-options "-Xmx12G" FilterMutectCalls -V ~/workspace/somatic/mutect/exome.vcf.gz \
    -O ~/workspace/somatic/mutect/exome_filtered.vcf.gz \
    -R ~/workspace/inputs/references/genome/ref_genome.fa

#Running mutect2 using gatk version 3.6
#java -Xmx12g -jar /usr/local/bin/GenomeAnalysisTK.jar -T MuTect2 --disable_auto_index_creation_and_locking_when_reading_rods -R ~/workspace/data/raw_data/references/ref_genome.fa -I:tumor ~/workspace/data/DNA_alignments/chr6+chr17/final/Exome_Tumor_sorted_mrkdup_bqsr.bam -I:Normal ~/workspace/data/DNA_alignments/chr6+chr17/final/Exome_Norm_sorted_mrkdup_bqsr.bam --dbsnp ~/workspace/data/raw_data/references/Homo_sapiens_assembly38.dbsnp138.vcf.gz --cosmic ~/workspace/data/raw_data/references/Cosmic_v79.dictsorted.vcf.gz -o ~/workspace/data/results/somatic/mutect/exome.vcf.gz -L ~/workspace/data/results/inputs/SeqCap_EZ_Exome_v3_hg38_primary_targets.v2.interval_list

echo ~/workspace/somatic/mutect/exome_filtered.vcf.gz > ~/workspace/somatic/mutect/exome_vcf.fof
bcftools concat --allow-overlaps --remove-duplicates \
    --file-list ~/workspace/somatic/mutect/exome_vcf.fof --output-type z \
    --output ~/workspace/somatic/mutect/mutect_exome_filtered.vcf.gz

tabix ~/workspace/somatic/mutect/mutect_exome_filtered.vcf.gz
```

### Now we are going to merge the variants detected from all three variant callers

The reason for merging is that the variant callers are working differently internally to identify somatic alterations, with different strengths and weaknesses. With outputs from all three algorithms, we can now merge the variants to generate a comprehensive list of detected variants:

```bash
# Unzip the vcf.gz files before combining Variants
cd ~/workspace/somatic
ls -halt ~/workspace/somatic/varscan/exome.vcf.gz
ls -halt ~/workspace/somatic/strelka/exome/results/variants/exome.vcf.gz
ls -halt ~/workspace/somatic/mutect/mutect_exome_filtered.vcf.gz

gunzip ~/workspace/somatic/varscan/exome.vcf.gz
gunzip ~/workspace/somatic/strelka/exome/results/variants/exome.vcf.gz
gunzip ~/workspace/somatic/mutect/mutect_exome_filtered.vcf.gz

#Need to change header sample names in vcf file produced by mutect2 in order to combine variants with those from other algorithms
sed -i 's/HCC1395BL_DNA/NORMAL/' ~/workspace/somatic/mutect/mutect_exome_filtered.vcf
sed -i 's/HCC1395_DNA/TUMOR/' ~/workspace/somatic/mutect/mutect_exome_filtered.vcf

# (UNIQUIFY command) java -Xmx4g -jar /usr/local/bin/GenomeAnalysisTK.jar -T CombineVariants -R ~/workspace/data/raw_data/references/ref_genome.fa -genotypeMergeOptions UNIQUIFY --variant:varscan ~/workspace/data/results/somatic/varscan/exome.vcf --variant:strelka ~/workspace/data/results/somatic/strelka/exome/results/variants/exome.vcf --variant:mutect ~/workspace/data/results/somatic/mutect/new_gatk_files/exome.vcf -o ~/workspace/data/results/somatic/exome.unique.vcf.gz
#java -Xmx24g -jar ~/workspace/bin/GenomeAnalysisTK.jar -T  CombineVariants -R ~/workspace/inputs/references/genome/ref_genome.fa -genotypeMergeOptions PRIORITIZE --rod_priority_list mutect,varscan,strelka --variant:varscan ~/workspace/somatic/varscan/exome.vcf --variant:strelka ~/workspace/somatic/strelka/exome/results/variants/exome.vcf --variant:mutect ~/workspace/somatic/mutect/exome.vcf -o ~/workspace/somatic/exome.merged.vcf.gz

java -Xmx12g -jar ~/workspace/bin/picard.jar MergeVcfs \
    -I ~/workspace/somatic/varscan/exome.vcf \
    -I ~/workspace/somatic/strelka/exome/results/variants/exome.vcf \
    -I ~/workspace/somatic/mutect/mutect_exome_filtered.vcf \
    -O ~/workspace/somatic/exome.merged.vcf

bgzip -c ~/workspace/somatic/exome.merged.vcf > ~/workspace/somatic/exome.merged.vcf.gz
tabix -p vcf ~/workspace/somatic/exome.merged.vcf.gz
```

### Left Align and Trim
The reason for left align the variants and trim then is explained [here](
https://genome.sph.umich.edu/wiki/Variant_Normalization#Left_alignment).
```bash
cd ~/workspace/somatic/

gatk --java-options "-Xmx12G" LeftAlignAndTrimVariants \
    -V ~/workspace/somatic/exome.merged.vcf.gz \
    -O exome.merged.leftalignandtrim.vcf \
    -R ~/workspace/inputs/references/genome/ref_genome.fa
```
Note that when running on chromosome 6 and 17 merged variants file, this gave 0 variants aligned.

### We will split multi-allelic variants into multiple records

```bash
cd ~/workspace/somatic/

vt decompose -s ~/workspace/somatic/exome.merged.leftalignandtrim.vcf \
    -o ~/workspace/somatic/exome.merged.leftalignandtrim.decomposed.vcf
```

### Basic Filtering on Somatic Variants

First, let's do a basic filtering for `PASS` only variants on our merged and normalized vcf file:
```bash
cd ~/workspace/somatic
gatk --java-options "-Xmx12G" SelectVariants -R ~/workspace/inputs/references/genome/ref_genome.fa \
    --exclude-filtered -V ~/workspace/somatic/exome.merged.leftalignandtrim.decomposed.vcf \
    -O ~/workspace/somatic/exome.merged.norm.pass_only.vcf
```

### Annotation with VEP. 
Again, we will use VEP to annotate the somatic variants as we did for the germline variants. 

```bash
cd ~/workspace/somatic
# Runtime: ~4min

#ssh -o ServerAliveInterval=300 -i course-setup-student-key.pem ubuntu@ec2-52-23-206-90.compute-1.amazonaws.com
#source .bashrc
#cd workspace/somatic/
nohup vep --cache --dir_cache ~/workspace/vep_cache \
    --dir_plugins ~/workspace/vep_cache/Plugins \
    --fasta ~/workspace/inputs/references/genome/ref_genome.fa --fork 8 \
    --assembly=GRCh38 --offline --vcf -i ~/workspace/somatic/exome.merged.norm.pass_only.vcf \
    -o ~/workspace/somatic/exome.merged.norm.annotated.vcf \
    --check_existing --total_length --allele_number  --no_escape --everything \
    --use_given_ref --force_overwrite &


#Copy the html VEP summary to your local computer
scp -ri ~/PEM_KEY_ID.pem ubuntu@AWS_ADDRESS_HERE:~/workspace/somatic/exome.merged.norm.annotated.vcf_summary.html .
```

Open the HTML file.

***Answer the following:***

- How many variants were processed by VEP?
- Any diffrence in the distribution of variants vs. the germline variants? 
- The variant categories e.g. intron_variant and regulatory_region_variant should be approached carefully, why?

## Inspecting variants in IGV
- Download the procdessed DNA bam files to your local machine if you have not done so yet, see the session "Introduction to IGV"

- Open the Exome_Tumor_sorted_mrkdup_bqsr.bam and Exome_Norm_sorted_mrkdup_bqsr.bam in IGV.

- Keep the VCF file open in a terminal window on AWS
- Try to find a relevant variant in TP53 in the VCF file on AWS
```bash
#Can be done in many ways, this is just an example.
cd ~/workspace/somatic
grep "TP53" exome.merged.norm.annotated.vcf | less -SN
#### Or 
less -SN exome.merged.norm.annotated.vcf
#In the less window type /TP53
```
- Why is the TP53 variant in three rows in the VCF file? 
    - lead: how many somatic variant callers were run?
- Go to the position (chr17 7675088) on IGV in your local machine
![](https://i.imgur.com/TQ8KKZ8.png)

    - Is it present in the normal DNA?
        - Do you still think it is valid? How come the variant can be present in the germline DNA?

- Do the same thing for a stop_gained and a frameshift variant
```bash
cd ~/workspace/somatic
grep "stop_gained" exome.merged.norm.annotated.vcf | less -SN
```
- The top variant is seen only once, this is just identified correctly by one caller because it changes two bases in a row
    - open chr6 1930220 in IGV

![](https://i.imgur.com/DbVcGjS.png)

- let us convert the vcf-file to a MAF file format
    - note all information is not kept but it simplifies looking at the variants
```bash
cd ~/workspace/somatic
perl ~/workspace/bin/mskcc-vcf2maf/vcf2maf.pl --input-vcf ~/workspace/somatic/exome.merged.norm.annotated.vcf --output-maf ~/workspace/somatic/exome.merged.norm.annotated.maf  --tumor-id TUMOR --normal-id NORMAL --inhibit-vep --ref-fasta ~/workspace/inputs/references/genome/ref_genome.fa

less -SN ~/workspace/somatic/exome.merged.norm.annotated.maf

#Count the number of variant types
cut -f9 ~/workspace/somatic/exome.merged.norm.annotated.maf | sort | uniq -c
#As you see the variant nomenclature is not the same in VCF and MAF

grep Nonsense_Mutation ~/workspace/somatic/exome.merged.norm.annotated.maf | less -SN 
```
- Have a look at the nosense variant in CCDC40 in IGV
    - chr17 80050160
    - What is the VAF?
    - What is the VAF of the TP53 variant?
    - Reflections?
        - Lead: clonality

- Can you find a variant in BRCA1? 
    - What is the impact? Is it relevant?

- There was a BRCA 1 variant in the germline that we did not dicsuss earlier
```bash
#Let us grep "BRCA1" and send the output to another grep command and take the "HIGH" impact variants
grep BRCA1 ~/workspace/germline/Exome_Norm_HC_calls.filtered.PASS.vep.vcf | grep "HIGH" | less -SN 
```
- let us inspect this variant in IGV, why is one allele missing in the tumor?
    - This is the variant location: chr17   43057078 
![](https://i.imgur.com/PXKeG7C.png)



















