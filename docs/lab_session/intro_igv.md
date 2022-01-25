
# Introduction to IGV

For various resons there is often a need to manually inspect sequencing data. For this purpose, we will use [IGV](https://software.broadinstitute.org/software/igv/home), the Integrative Genomics Viewer. IGV is in a nutshell an easy-to-use and interactive tool for the visual exploration of genomic data (not just sequening data). In IGV it is possible to integrate multiple data types in a straight forward manner. 

Download the IGV version that includes java.

```bash
#First download IGV to your local computer. 
#Using a mac, execute the igv.sh bash script (obs - change to the correct path)
/PATH_GOES_HERE/IGV_DIRECTORY/igv.sh

#Go to your favourite directory and download a bam file for the IGV intro:
scp -ri ~/PEM_KEY_ID.pem ubuntu@AWS_ADDRESS_HERE:~/workspace/igv_intro_bams/* .
```

## HCC1143 cell line

The cell line sequencincg data that will be used to explore the features of IGV comes from a 52 year old caucasian woman with breast cancer. Additional information can be found at the [ATCC](https://www.atcc.org/products/crl-2321) (American Type Culture Collection) website. The bam-file contains only reads on chromosome 21:19,000,000-20,000,000 in order to reduce file sizes.

## Getting familiar with IGV

- IGV comes pre-loaded with the genome build HG19. During this introduction we will use **hg19**.
- Remote data files (tracks) can be loaded from the IGV server. 
    - Select "File" & "Load from Server...":
        - Ensembl Genes 
        - GC Percentage
        - dbSNP (latest version)
        - Repeat masker

![](https://i.imgur.com/hReiSrn.png)

- Load the bam file
    - Select "File" & "Load from File...":
        - Locate "HCC1143.normal.21.19M-20M.bam" on the hard drive and press "Open".

### Familiarize the different sections of IGV

- The top genome ruler.
- The top panel data tracks.
- The mid panel sequencing data tracks.
- The bottom panel data- and annotation tracks.

- Note the gene model.
    - line with arrows: the direction/strand of the gene
    - the thin boxes: untranslated regions
    - the thick boxes: coding exons
 
![](https://i.imgur.com/9XX7ncG.png)

### Investigate read alignments

- Go to chr21:19,479,237-19,479,814.
- Right click on the aligned reads in the mid panel sequencing data tracks.
- Test to:
    - Right click on a read and sort alignments by start location.
    - group alignments by pair orientation.
![](https://i.imgur.com/umnqsyZ.png)

- Experiment with the settings.
- Click on an individual read and assess the information.
![](https://i.imgur.com/VWfQbuk.png)
- Enable the ruler, it simplifies investigating variants. Push the button to the upper right, see figure below:
![](https://i.imgur.com/9BXgPzm.png)

### Investigate a SNP

- Go to chr21:19479237-19479814
- Note a heterozygous variant at chr21:19479321.
- Zoom in to observe the individual base changes.
- Sort aglinment according to base.
- Color alignments according to strand.
- Try to answer the following:
    - Reflect on mapping and base qualities, is this a high quality SNP?
    - How does the shaded non-reference bases help?
    - How does read strand information help?
![](https://i.imgur.com/4981oRj.png)

### Investigate a homopolymer region with variability

- Go to chr21:19518412-19518497.
    - Group alignments by read strand.
    - Sort alignments by base.
- Try to answer the following:
    - Is this a region with variability that can be trusted, motivate?
- Go to chr21:19518470
    - Sort alignments by base.
- Try to answer the following:
    - Is this a region with variability that can be trusted, motivate?
![](https://i.imgur.com/IS7JXrO.png)

### GC precentage and coverge

- Go to chr21:19611925-19631555.
    - In the menue bar, go to View > Preferences > Alignments tab > set the visibility range threshold to 20 kb.
    - Use Collapsed view for the read-pairs.
    - Group alignments by "none".
    - Sort alignments by start position.
    - Color alignments by -> insert size and pair orientation.
- Is GC precentage and coverage correlated?

### SNPs on different alleles

- Go to chr21:19666833-19667007.
- Sort by base at position chr21:19666901
- Try to answer the following:
    - Are these valid SNPs?
    - Explan what you are observing.
![](https://i.imgur.com/UwAjaMe.png)
### Problematic region

- Load a the repeat masker remote data track.
    - Select "File" & "Load from Server...":
        - Repeat masker
- Go to chr21:19800320-19818162.
- Try to answer the following:
    - Why are some reads white?
    - Explan what you are observing.
![](https://i.imgur.com/GQrfdLF.png)

### Deleted region

- Go to chr21:19324469-19331468.
- Select expanded view.
- View reads as pairs.
- Color alignments by insert size and pair orientation.
- Sort reads by insert size.
- Try to answer the following:
    - Click on a red read pair and a grey read pair, can you identify information supporting the deletion in the red pairs? 
    - Is the deletion hetero- or homozygous?
![](https://i.imgur.com/IohithX.png)



### Mis-alignment

- Go to chr21:19102154-19103108.
- Group alignment by "Chromosome of mate".
- Zoom out in a stepwise manner and compare the region of interest with surrounding regions.
- Try to answer the following:
    - Click on one of the brown pairs, where does the mate map?
    - What is the reason for an increase in poorly mapped reads?

## Download the tumor and normal bam files processed today

Go to your favourite directory
```bash
# Tumor bam and bai files
scp -ri ~/PEM_KEY_ID.pem ubuntu@AWS_ADDRESS_HERE:~/workspace/align/Exome_Tumor_sorted_mrkdup_bqsr.ba* ./

# Normal bam and bai files
scp -ri ~/PEM_KEY_ID.pem ubuntu@AWS_ADDRESS_HERE:~/workspace/align/Exome_Norm_sorted_mrkdup_bqsr.ba* .




# The in solution hybridization based capture definition files
scp -ri ~/PEM_KEY_ID.pem ubuntu@AWS_ADDRESS_HERE::~/workspace/inputs/references/exome/probe_regions.bed ./
scp -ri ~/PEM_KEY_ID.pem ubuntu@AWS_ADDRESS_HERE:~/workspace/inputs/references/exome/exome_regions.bed ./
```

### Open the files in IGV

- Start a new IGV session
    - Select File -> New Session
    - Select hg38
    - Load the bam files
        - Exome_Norm_sorted_mrkdup_bqsr.bam
        - Exome_Tumor_sorted_mrkdup_bqsr.bam
    - Load the cature definion files
        - exome_regions.bed
        - probe_regions.bed
- Change viewing preferences
    - In the menue bar, go to View > Preferences > Alignments tab > set the visibility range threshold to 25 kb.
    - Use Expanded view for the read-pairs.
    - Group alignments by “none”.
    - Sort alignments by start position.
    - Color alignments by -> no color.
    - Type TP53 in the search box and press enter.
- Try to answer the following:
    - Do you observe any variants?
        - If yes, what kinds?
        - Are they coding/non-coding?
        - Any somatic variants?
        - Try to assess the consequence of a selected variant.
    - Why do you think there is a difference between the probe- and exome regions?
        - Assess the coverge pattern relatively the probe regions. Try to explain the pattern.

![](https://i.imgur.com/DXha8Ct.png)














































