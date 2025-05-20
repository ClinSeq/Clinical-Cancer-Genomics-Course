---
title: Clinical Cancer Genomics Project Exercise - Analysis of Gene Panel Sequencing Data
---

## Overview

This exercise aims to provide you with practical experience in analyzing clinical cancer genomics data. You will investigate six cases of gene panel sequencing, focusing on copy number variations (CNVs) and small somatic variants (snvs and indels). The goal is to visualize these data, identify tumor-relevant variants and phenotypes, and address quality control (QC) issues within the context of each sample.

## Objectives

- Explore results of gene panel sequencing in clinical cancer genomics.
- Develop skills in data visualization and interpretation.
- Identify tumor-relevant phenotypes and genomic alterations and relate these findings to potential clinical implications.

## Tools and Resources

- **Gene Panel Sequencing Data**: Six anonymized cell-free DNA (cfDNA) cases are provided [here](https://kise-my.sharepoint.com/:u:/g/personal/johan_lindberg_ki_se/EZQlPqn2JkJEhkO0QTLt7koBLIGLG9rXa1jhDPR-LxnVeQ?e=W5MwXY), each including copy number and small somatic variant data:
    * Annotated somatic variant calls in the [MAF](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) format,
    * Copy number analysis data as two files, **bins** and **segments**.
- **Processing environment**: Choose a suitable coding environment for your work, such as [R/Quarto](https://quarto.org/) or [Python/Jupyter](https://jupyter.org/).
- **Relevant databases**: You may find these resources useful: [cancerhotspots.org](http://cancerhotspots.org), [OnkoKB](https://www.oncokb.org), [COSMIC](https://cancer.sanger.ac.uk/cosmic), [Tumorportal](http://tumorportal.org/) and [cBioBortal](https://www.cbioportal.org/).
- **Computational Resources**: This work will not be computationally intensive. You can work on your own computer.
- **Sensitive Data**: Although these data are anonymized, you must not redistribute the data files and you must not retain a copy of the data files after completing the course.


## Instructions

For each sample you will need to write a summary of your observations, with the following points addressed:

- Sample name
- A visualization of the result
- Any apparent QC issues with result (not including basic sequencing QC, as you will not have the raw sequencing data available)
- Estimated tumor cell or tumor DNA content
- Obvious or suspected mutational phenotypes
- Apparent "driver" variants and for each such variant whether it appears to be
    - clonal
    - featuring a second hit
    - activating or deactivating the gene(s)
    - treatment-predictive in prostate cancer
    
Collect the result into a single pdf file. Make sure your name(s) are included, and upload the result into your directory on the server. 

## Hints

* Take a look at the databases listed above. You can either check individual genes or variants manually, or you can download for example the hotspot specification as an excel (don't take the MAF) file.

* You can code your analysis similarly to how you worked in the copy number exercise. All files can be parsed as tables. You can choose whether to code everything and generate a pdf report, or just put results in e.g. Word and export the final version as pdf.

* If you have problems figuring out how to do something, feel free to ask the teachers so that you do not get stuck. But first try ChatGTP or something similar!

