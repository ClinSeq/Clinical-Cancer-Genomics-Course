---
title: Introduction to IGV

---

## Integrative Genomics Viewer (IGV)

As we perform bioinformatics analyses, there is often a need to manually inspect sequencing data and our analysis results. For this purpose, we can use the [Integrative Genomics Viewer (IGV)](http://igv.org/). IGV is an easy-to-use and interactive tool for the visual exploration of genomic data. IGV also allows us to integrate multiple data types in a straight forward manner, as we will see later.

### Getting started

If you have not done so already, we need to start by downloading the IGV Desktop Application onto our local computers. Go to <https://igv.org/doc/desktop/#DownloadPage/> and select the correct IGV version for your computer (MacOS, Windows, or Linux) and complete the installation. **Make sure to select a version that includes Java.**

We will use sequencing data from the HCC1143 cell line to explore the features of IGV. This cell line comes from a 52-year-old Caucasian woman with breast cancer. Additional information on the cell line can be found at the [American Type Culture Collection (ATCC)](https://www.atcc.org/products/crl-2321). The sequencing data is in the form of a [binary alignment map (BAM) file](https://en.wikipedia.org/wiki/Binary_Alignment_Map) and contains sequencing reads only from chromosome 21:19,000,000-20,000,000 in order to reduce file size.

Start by downloading a BAM file of HCC1143 DNA-sequencing data aligned to the reference genome hg19 to your local computer.

```bash
# Go to a directory on your computer where you would like to download the file using ´cd´
# Then, download the BAM file and its index file for this IGV introduction
wget https://rnabio.org/assets/module_2/HCC1143.normal.21.19M-20M.bam
wget https://rnabio.org/assets/module_2/HCC1143.normal.21.19M-20M.bam.bai
```

### Getting familiar with IGV

IGV comes pre-loaded with the human reference genome assembly hg19. During this introduction we will therefore use **hg19**.

Remote data files (tracks) can be loaded from the IGV server for visualization. Let's load some data tracks as shown in the image below. Select "File" and then "Load from Server...". Select the checkboxes as shown below and press OK to load the tracks for Ensembl Genes, GC Percentage, dbSNP (latest version), and Repeat Masker.

![](https://i.imgur.com/hReiSrn.png)

Then, let's load the BAM file from above into our IGV session by selecting "File" and then "Load from File...". Locate the "HCC1143.normal.21.19M-20M.bam" file on your computer, click on it, and press "Open".

### IGV sections

You'll notice that the IGV view consists of several sections. Locate and explore the following:

-   The top genome ruler.
-   The top panel data tracks.
-   The mid panel sequencing data tracks.
-   The bottom panel data- and annotation tracks.
-   The gene model:
    -   line with arrows: the direction/strand of the gene
    -   the thin boxes: untranslated regions
    -   the thick boxes: coding exons

![](https://i.imgur.com/9XX7ncG.png)

### Investigating read alignments

-   Go to chr21:19,479,237-19,479,814.
-   Click on an individual read and see what information is provided.
-   Right click on the aligned reads in the mid panel sequencing data tracks.
-   Try out the following and see what changes in the IGV view:
    -   Right click on a read and sort alignments by start location.
    -   Group alignments by pair orientation.

![](https://i.imgur.com/umnqsyZ.png)

-   Enable the ruler- this simplifies investigating variants in IGV. Push the button in the upper right, see figure below.

![](https://i.imgur.com/9BXgPzm.png)

### Investigating a single nucleotide polymorphism (SNP)

Single nucleotide polymorphisms (SNP) are germline substitutions of a single nucleotide at a specific position in the genome. Our cell line has a SNP at **chr21:19479321**. Go to this location in IGV and do the following:

-   If needed, zoom in to observe the individual base changes at chr21:19479321.
-   Sort the alignments according to base.
-   Color alignments according to strand.
-   Answer the following:
    -   Is the SNP heterozygous or homozygous?
    -   A good quality variant should be supported by reads on both the plus and minus strands. Is this the case for the SNP?

![](https://i.imgur.com/4981oRj.png)

### Investigating a homopolymer region with variability

Homopolymers are stretches of mono nucleotides (DNA bases) greater than two bases long which occur together. For example, "ATCCCCCCCCCA" has a homopolymer of length 9 (base "C"). One homopolymer region in the human genome is located around **chr21:19518470**. Go to this region in IGV and do the following:

-   Group alignments by read strand.
-   Sort alignments by base.
-   Answer the following: Is this a region with variability that can be trusted and why? In other words, is the variability we see in this region real or an artefact?

![](https://i.imgur.com/IS7JXrO.png)

### Investigating a deleted region

In genetics, a deletion is a mutation (a genetic aberration) in which a part of a chromosome or a sequence of DNA is left out during DNA replication. Our cell line has a deleted region at **chr21:19324469-19331468**. Go to this region in IGV and then do the following:

-   Select expanded view.
-   View reads as pairs.
-   Color alignments by insert size and pair orientation.
-   Sort reads by insert size.
-   Answer the following:
    -   Look at the red read pairs at the top of the screen- can you identify information supporting the deletion in the red read pairs? What about some of the grey read pairs?
    -   Is the deletion hetero- or homozygous?

![](https://i.imgur.com/IohithX.png)

### Investigating a problematic region

Finally, let's take a look at one problematic region in the genome and what the data looks like in the region for our cell line. This region is located at **chr21:19800320-19818162**. Go there in IGV and answer the following questions:
-   Why are some of the reads white in this region?
-   Why is this a problematic region?

![](https://i.imgur.com/GQrfdLF.png)
