# GENOMIC ALTERATIONS INTERPRETATION -- EXERCISES

By David Tamborero.

## (1) Exercise with genomic variants knowledgebases: ClinVar

There are several international projects devoted to gather knowledge associated to gene variants according to the results that are continuously generated by clinical and preclinical studies. These efforts are consolidated in different databases (aka knowledgebases). Each knowledgebase may have a different scope and a distinct process to gather, annotate and access the data. Understanding these details is important for an accurate use of the resources.

Some of these knowledgebases are more focused on **germline variants and their relevance to cause (or predispose) to diseases** (including cancer, i.e. variants that increase hereditary risk of cancer). One of these databases is [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), a project under the umbrella of the (American) National Center for Biotechnology.

The ClinVar website allow to query information associated with a particular variant* of your interest by using different formats, including variant identifiers that are provided by third-party associations (e.g. dbSNP identifiers). The more common method to query for a variant is by introducing the exact genomic coordinates. ClinVar classifies variants by their pathogenicity, in five levels (pathogenic, likely pathogenic – benign, likely benign – variant of unknown significance) following the recommendations of the American College of Medical Genetics and Genomics and the Association for [Molecular Pathology](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/). These guidelines include criteria that are specific to the germline context, such as the assessment of the variant inheritance model in families of the individual presenting the phenotype of interest.

* the search box also allows to look for all the genomic variants in a given gene or for a given disease.

**Exercise**: search for the gene **CFTR**, which encodes a protein that functions as a channel across the cell membrane that produce mucus, sweat, saliva, tears, and digestive enzymes. Germline pathogenic variants in this gene are demonstrated to lead to cytstic fibrosis. As a result of querying this gene in the ClinVar database, you will see a new page with a table listing all the variants observed this gene and having phenotypes annotated in the resource.

Select only those variants that have been classified as **‘Pathogenic’** by using the corresponding filter in the ‘Clinical significance’ section located at the left side of the page (see image below).

![](https://i.imgur.com/LDehTTJ.jpg)

Importantly, ClinVar acts as a **repository** of assertions, but there is no in-house team devoted to manually review these assertions (as opposite to other databases in which ‘editors’ review and ‘approve’ the content before releasing it). Instead, ClinVar provides a score (from zero to four stars) that can be used as a proxy on the **‘quality consensus’ of the entry**. For instance, one star is automatically given by the system when several curators have submitted an assertion for the same variant with contradictory interpretations. Two stars are given if the interpretation of the variant coincides across the curators that have submitted information for that variant. And three stars means that the assertion has been submitted by a certified committee of experts, and thus is considered a consensus expert recommended interpretation. Maximum score (four stars) is given for variant effects that are currently used in the standard clinical practice (which can be only submitted by the corresponding authorities in the ClinVar repo).


In the CFTR results page provided by ClinVar, select now the filter **‘Practice guidelines’** in the ‘Review status’ section located at the left side of the page. This should limit the number of displayed variants to ~20-25. Note in the ‘Condition(s)’ column of he results table that most of these variants are indeed described as pathogenic for Cystic Fibrosis. If you click in any of them, you will open a new page with the details of this variant and associated information (see image pasted below). Note also how this ClinVar page includes links to other resources with more information of the variant. The information provided by ClinVar includes (you need to scroll down in the page) who has submitted the variant assertion and the published references supporting the variant effect.

![](https://i.imgur.com/aFpKleR.jpg)

ClinVar is a repository widely employed to look for information of variants associated with Mendelian diseases. It also contains information of **germline variants associated to cancer**, e.g. cancer predisposing variants. You can repeat the search for **BRCA1** or **BRCA2** (genes whose pathogenic germline variants predispose to ovarian and breast cancer), or **MLH1, MSH2** or **MSH6** (genes whose germline pathogenic variants predispose to colorectal cancer), as examples of the latter.

(NOTE about the page interface: be aware that when you search for a new gene in ClinVar, the results are likely going to be filtered by the last applied filter conditions (so it s not a new search from scratch)-- which can be a bit unnoticed as the active left bar options are not that obvious to notice).

Please note that this example is focused on querying variants annotated in the ClinVar database for a given gene, but it is also possible to query for information (if any) of a specific variant of your interest. For that, you will need to use standarized variant nomenclature, as HGVS (e.g. NM_000314.4:c.395G>T); for more details, visit this [link]( https://www.ncbi.nlm.nih.gov/clinvar/docs/help/#quick_start)

## (2) Exercise with genomic variants knowledgebases: OncoKB

**OncoKB** is a knowledgebase curated by Memorial Sloan Kettering members and focused on cancer mutations (https://www.oncokb.org/). This resource provides classification of the oncogenic effects of variants in cancer genes and their value (if any) as biomarkers of prognosis, diagnosis and drug response. This database is curated manually and thus all the entries have been approved by their Scientific Content Management Team. Therefore, OncoKB does not provide a measure of the quality of the assertions, as all the entries passed their quality criteria. This can not confused with the level of evidence supporting the described effect; for instance, one biomarker can be supported by clinical studies or by pre-clinical data, which represents two distinct levels of clinical evidence; but in both cases, the quality of the supporting data is sound (according to their criteria). In detail, the classification that OncoKB follows to classify the actionability of the biomarkers can be found here: https://www.oncokb.org/levels

The searches can be done at the level of a particular variant, or you can query for all the variants with associated information available for a given gene (see image below).

![](https://i.imgur.com/JjfB4WE.png)


You can look for your gene of interest, e.g. FGFR3; you will see the assertions divided in annotated alterations, therapeutic, prognostic and FDA approved (in this case, there is no “diagnostic” entry associated to this gene).


![](https://i.imgur.com/B3tGo97.png)


Annotated alterations summarizes the variants in this gene with known effect, classified as oncogenic/likely oncogenic or Likely neutral. Note that for oncogenic variant, the annotation includes the particular effect (gain-of-function in this case, as FGFR3 is an oncogene). Therapeutic entries are the subset of oncogenic variants that are also known to be biomarkers of drug response, following the aforementioned classification rank. Note that a biomarker must include the cancer type(s) and the drug(s) in which that effect has been described. Note that the FDA-recognized content is a recent (October 2021), incorporation in oncoKB, in which the entries labelled as so are recognized as FDA approved information. Click in an entry to see more details, including the supporting literature employed to curate the information.

![](https://i.imgur.com/DSOPAe4.png)

## (3) Use of a (lightweight) clinical decision support system
Molecular profiles of tumors can inform treatment decisions; most of the therapeutic interventions are currently guided by genomic profiling (specially mutations) and -to a lesser extent- expression signatures. Performance of NGS assays is common in the oncology setting, both in investigational as well as routine care. The variants observed by the assay must be then interpreted to first identify those that are functionally relevant for the tumor individual, and then which of them are also clinically relevant according to current evidence. To do this process manually is time consuming and prone to inaccuracies or even errors. Clinical decision support systems (CDSS) implement efficient processes to annotate and classify the variants using a variety of computational tools and databases. As the complexity of molecular-guided therapies continues to grow, the use of CDSS is increasingly important.

There are several commercial CDSS available at present, but larger academic centers also develop their own in-house solutions, which may better adapt to their specific needs and investigational initiatives. One of this academic solutions, the **Molecular Tumor Board Portal** (MTBP), released an open version available for the community (http://mtbp.org) This public portal is a lightweight* version of the ‘production’ instances used in real clinical project. Among others, the public portal provides a limited analysis of the data (accepts only mutations, with no information about mutation signatures, and do not incorporate information about origin/clonality of the mutations neither generate clinical flags such as allocation to clinical trials which requires additional clinical/pathology information). This public resource is limited to only provide a generic functional and predictive interpretation of a list of variants uploaded by the user for research purposes only.

*note the disclaimer in the website: as the public portal is considered a research tool and must be used in that context!*

Open the MTBP public version website http://mtbp.org; remind that this is a free tool, and you do not even need to login to be able to make an analysis (however, if you want to keep the results of your analysis in a repository, you need to register with a valid email – so the system can associate the


repository to that account). With or without being logged in, you can click the ‘See examples’ or the ‘Analyse your variants’ button (see next image):

![](https://i.imgur.com/han8Oen.png)


‘See examples’, as the name indicates, will show the reports that are obtained by analysing some variants observed in three different tumor types used as examples here. If you click ‘Analyse your varaints’, you will be able to perform an analysis for your own variants of interest. This will open a new window; you will need to introduce a name to identify your analysis in the system (step 1), select the tumor type of the sequenced sample (step 2), and then upload the variants (step 3).

Note that there is a box in step 3 to include an OncoKB token – this is an optional step, which allows the system to query more comprehensively that database if the user is registered to that resource (and thus obtained a token; see the question mark for more details). Keep this field empty if you do not have that token. In addition, the variants can be passed by either uploading a VCF file or by typing the variants in a free text box. The VCF file must follow standards (including the use of the ‘PASS’ label in the filter column to indicate the variants that are to be included in the analysis) that can be annotated in hg19/GRCh37 or hg38/GRCh38. Variants typed in the free text box must follow HGVS standarized syntax. Please use the question marks to open popups with more information.

When ready, click the ‘Run’ button:

![](https://i.imgur.com/m72ePgr.png)


If everything goes well (including the validation of the variants’ format), the information will be sent to the system and processed during a while; the results of the analysis will be available in the system when finished. Note that the resulting report contains:
- some descriptives of the gene and the variant (such as the protein coordinates and the consequence type)
- the evidence (if any) supporting the functional relevance classification
- the evidence (if any) supporting the role of the variant as cancer biomarker
    
Note that the report is interactive, and you can access to furtehr details, including the original sources of evidence, by clicking around in the document. Remember also that the evidence supporting the variant’s functional relevance is classified in three families:

- evidence of type A: the effect is supported by manually curated data (i.e. results of clinical or experimental studies)
- evidence of type B: the effect can be assumed by using well sedimented biological knowledge (i.e. to assume that a premature stop codon before the last exon-exon junction will trigger nonsense-mediated decay)
- evidence of type C: the effect is based on a bioinformatics prediction, which is the lowest level of evidence.

Also, the evidence supporting the biomarker is tiered by several factors, including the match of the patient and biomarker cancer type and the clinical status of the biomarker effect (e.g. it is accepted in clinical guidelines versus is supported by preliminary results of clinical trials) following ESMO recommendations.

![](https://i.imgur.com/9GLqfJ0.png)

Use the VCF files distributed for the course (or use variants of your own interest) to test the system and discuss the results. Again, remember that the open MTBP is a lightweight version of the system, and thus the analyses are limited to the generic interpretation provided here.

VCFs can be downloaded here:

[BREAST_sample_muts.hg19.vcf](https://course-cg-5534.s3.amazonaws.com/vcf/BREAST_sample_muts.hg19.vcf)

[LUAD_sample_muts.hg19.vcf](https://course-cg-5534.s3.amazonaws.com/vcf/LUAD_sample_muts.hg19.vcf)








