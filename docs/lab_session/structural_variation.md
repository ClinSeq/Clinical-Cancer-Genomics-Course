# Calling structural variation

In this practical session we will look into structural rearrangements. The cfDNA sequencing was performed on commercial biomaterial from a patient with mCRPC. 

The in-solution hybridisation based capture assay, developed for the clinical trial (ProBio)[https://www.probiotrial.org] was applied.

## The ProBio Assay
---
![](https://i.imgur.com/Onu5vYJ.png)

**Assay explained:**
![](https://i.imgur.com/HLDUR32.png)

---


Go to the folder that contain the ProBio ctDNA data from the mCRPC commercial reference sample

## Create mini-bam files and run Delly
```bash
cd ~/workspace/svs

# First create smaller bam files for faster 
# processing since we already know where the structural variants are (on chromsome 10 and 21).

# To create the smaller bam file we will be using bedtools, 
# one of the most useful tools in genomics. 
# The mapped sequencing is fed into intersectBed and all 
# reads overlapping with chr10 or chr21 are kept. The rest are discarded.
# This is an example where the genome build needs to be known. 
# Are we using chr10 or just 10? How can you find out?
# Where is this info to be found at the terminal?
# Clue - its a very quick thing to check in the current folder.


# Create a bed file that can be used to slice out the desired chromsomes

#Create an empty file and then add tab.delimited lines
echo -ne "" > chr.bed
echo -e "10\t1\t135534747" >> chr.bed
echo -e "21\t1\t48129895" >> chr.bed

#check the file
cat chr.bed

#OBS - the commands below will be using nohup which
# means that the commands will be sent and running
# in the backgrond. You should immediately get the 
# prompt back after running nohup

# This takes a couple of minutes
BAM=./bams/RB-P-ResBio12-CFDNA-sample12-KH-C3-nodups.bam
nohup intersectBed -a $BAM -b chr.bed > T_mini.bam  2> nohup_Tmini.log &

BAM2=./bams/PB-P-HD2-N-1811-KH-C3-nodups.bam
nohup intersectBed -a $BAM2 -b chr.bed > N_mini.bam 2> nohup_Nmini.log &

#If you want to check what is running in the background type
# "htop" and you will see an overview of the processes running
# and the CPU load of your AWS instance.
# Quit "htop" by pressing "q".

#After the chr selection has been performed, check which chromosomes that re in the mini-bam file
samtools view N_mini.bam | cut -f3 | sort | uniq -c
samtools view T_mini.bam | cut -f3 | sort | uniq -c


#Index the bam files
samtools index T_mini.bam
samtools index N_mini.bam

# If you want to run dell on the entire bam file - takes one hour
#nohup delly call -o ./RB-P-ResBio12-CFDNA-sample12-KH-C3_PB-P-HD2-N-1811-KH-C3_somatic_delly.bcf -g ./human_g1k_v37_decoy.fasta ./bams/RB-P-ResBio12-CFDNA-sample12-KH-C3-nodups.bam ./bams/PB-P-HD2-N-1811-KH-C3-nodups.bam > delly_nohup.log &

#Now run Delly on the samll bam files, write to a .bcf file
nohup delly call -o ./mini.bcf -g ./human_g1k_v37_decoy.fasta ./T_mini.bam ./N_mini.bam > delly_nohup.log &

#On the entire bam file
#bcftools convert -O v -o ./RB-P-ResBio12-CFDNA-sample12-KH-C3_PB-P-HD2-N-1811-KH-C3_somatic_delly.vcf ./RB-P-ResBio12-CFDNA-sample12-KH-C3_PB-P-HD2-N-1811-KH-C3_somatic_delly.bcf

bcftools convert -O v -o ./mini.vcf ./mini.bcf
```

## Filtering using bcftools

Now we will filter with some basic filters to remove false positive variants
```bash
cd ~/workspace/svs/
# RB-P-ResBio12-CFDNA-sample12-KH-C3_delly.vcf 
# to check total number of variants 

#If processing the entire bam files.
#VCF_FILE=RB-P-ResBio12-CFDNA-sample12-KH-C3_PB-P-HD2-N-1811-KH-C3_somatic_delly.vcf
VCF_FILE=mini.vcf

grep -v "^#" $VCF_FILE  -c
# 49

# 1. Filter only "PASS" variants 
bcftools view -f "PASS"  $VCF_FILE | grep -v "^#" -c 
#37
bcftools view -f "PASS"  $VCF_FILE > VCF_NEW.vcf

# 2. Filter by INFO column PE - Paired-End read SR - Split Read support 
bcftools filter -i 'INFO/PE>15 | INFO/SR>15' VCF_NEW.vcf | grep -v "^#" -c 
#5
bcftools filter -i 'INFO/PE>15 | INFO/SR>15' VCF_NEW.vcf > VCF_NEW_FILTERED.vcf
#have a look at the output
less -SN VCF_NEW_FILTERED.vcf
```

## Download and view the structural variants


Now we will filter with some basic filters to remove false positive variants
```bash
# Download the output to your computer
scp -ri ~/PEM_KEY_ID.pem ubuntu@AWS_ADDRESS_HERE:~/workspace/svs/*.bam .
scp -ri ~/PEM_KEY_ID.pem ubuntu@AWS_ADDRESS_HERE:~/workspace/svs/*.bai .
scp -ri ~/PEM_KEY_ID.pem ubuntu@AWS_ADDRESS_HERE:~/workspace/svs/VCF_NEW_FILTERED.vcf .

#Now we will download a couple of files with data supporting the variants from an internally developed tool for identifying structural rearrangements
wget https://course-5534.s3.amazonaws.com/svs/evidence_svs.sort.bam
wget https://course-5534.s3.amazonaws.com/svs/evidence_svs.sort.bam.bai
wget https://course-5534.s3.amazonaws.com/svs/sample-cfdna-svcaller-DEL.gtf
wget https://course-5534.s3.amazonaws.com/svs/sample-cfdna-svcaller-TRA.gtf
```


## Investigate the variants in IGV

- Remember - these files are mapped to hg19
- Open the following files in IGV:
    - T_mini.bam
    - N_mini.bam
    - VCF_NEW_FILTERED.vcf
    - evidence_svs.sort.bam
        - this bam file contains reads that are the output from a internally developed variant caller (SVcaller)
    - sample-cfdna-svcaller-DEL.gtf
    - sample-cfdna-svcaller-TRA.gtf
        - The gtf-files mark the variants identified by SVcaller


It is possible to select multiple files at the same time by pressing CMD (mac, use other key on windows, likely CTRL)
![](https://i.imgur.com/EiG2V7J.jpg)


You should now see something like this:
![](https://i.imgur.com/PqkzmDx.png)

This sample was analysed using multiple structural variant callers and two variants are known

- One variant in PTEN
- One variant on chr21 resulting in the TMPRSS2-ERG gene fusion
- Write PTEN and press "GO" or just "ENTER"

![](https://i.imgur.com/HGY0g7R.png)

- Was the variant detected by Delly?
- Zoom in on the variant in PTEN 

![](https://i.imgur.com/hDNzQGG.png)

- For the T_mini.bam
    - Color alignments "no color"
    - Group alignments "none"
    - Sort alignments "base"

![](https://i.imgur.com/uVy0UZX.png)

- Right click on a read in the evicence_svs.sort.bam
    - Select "View mate region in split screen".
    - This is a very useful way of seing both ends of a structural variant.

![](https://i.imgur.com/cBZCmD3.png)


![](https://i.imgur.com/Kos5iz7.png)


- Try to answer the question in the screenshot above.


- Now, have a look at the TMPRSS2-ERG fusion
    - Go here: chr21:39,766,196-42,940,331
    - Is this variant found by Delly? How do you see that?
    - Zoom in on the left end, how many reads are supportin the variant according to Delly?
        - Click on the VCF file and look for PE (paired-end) and SR (split-read) info.
![](https://i.imgur.com/lvw5Tt1.png)

