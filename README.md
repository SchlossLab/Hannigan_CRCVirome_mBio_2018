# Hannigan-2016-ColonCancerVirome
Investigating the gut virus communities associated with colon cancer.

##Study Outline
###Introduction
Sample source is the colon cancer 30-30-30 study that was originally published by Zackular JP.

> Zackular, J. P., Rogers, M. A. M., Ruffin, M. T. & Schloss, P. D. The human gut microbiome as a screening tool for colorectal cancer. Cancer Prev Res (Phila) 7, 1112–1121 (2014).

###Sample Summary
The samples were archived at -20C and are able to be used for virome and whole metagenome sequencing. The bacterial communities have already been characterized by 16S in the cited paper.

* Healthy: 30
* Colonic Adenoma: 30
* Colonic Adenocarcinoma: 30

###Processing and Sequencing Approach
I used my optimized virome purification protocol with some modifications for human samples. This approach will only look at genomic DNA from virus like particles (VLPs) and will not detect genomes integrated into host cell genomes. In other words, it is evaluating the *actively infecting* virome. The **Powersoil** plate extraction protocol was used for the whole metagenome processing, which is all of the DNA contents within the stool.

[Virome Extraction Protocol](https://github.com/Microbiology/HanniganNotebook/blob/master/notebook/protocols/protocols/ColonCancerFecalViromePurification.md)

Sequencing library preparations were done using the Illumina **NexteraXT** prep kit. The samples will be sequenced on the HiSeq 4000 with one lane devoted to the virome library and the other devoted to the whole metagenome library.

###Aims
#### 1. Identify Specific Viruses Associated with Colon Cancer
*Rationale*: Viruses have causative roles in a variety of cancers, including human papillomavirus and human polyomavirus. This is because viruses can silence and/or mutate host genes that alter host functionality. Colon cancer may also have a virus association.

*Hypothesis*: A subset of gut viruses are associated with colon cancer.

*Approach*: I will identify viruses using assembled contigs from the total virus dataset. Viruses will be identified as eukaryotic viruses and phages using both my ORF taxonomic annotation approach and my newest netowrk prediction model. "Virus-Cancer Association" will be defined as a statistically significant increase in relative abundance in the cancer/ademona samples comapared to the healthy controls. Abundance will be calculated by quantifying reads that map to contigs.

#### 2. Identify virus diversity signatures that are associated with colon cancer.

*Rationale*: Bacterial community signatures, including diversity, are associated with colon cancer developement and progression. Previous studies have indicated that altered bacterial diversity indicates an altered diversity in other microbial communities, including fungi and viruses. Often decreased bacterial diversity is associated with increased phage/virus diversity. It is worth noting that my preliminary data actually suggests the opposite to be true, and disease is associated with decreased viral diversity, which suggests a unique virus diversity signature (compared to the few other studies).

*Hypotheses*: Virome diversity is greater in colon cancer specimines compared to healthy controls.

*Approach*: I will evaluate virus diversity using three techniques. I will first be using a relative abundance table of contig abundance in each sample to calculate the alpha diversity associated with each sample, and the beta diversity dissimilarities between samples. Because this is a reference independent approach, I will supplement it with a reference dependent approach where contig "OTUs" will be defined by their taxonomic identifications. Finally I will use operational protein families (OPFs) instead of contigs to calculate the functional diversity, which is becoming th estandard in the field anyways.

#### 3. Classify colon cancer samples using a virus-based machine learning algorithm

*Rationale*: Bactieral community signatures have been used to classify stool samples as healthy, cancerous, or pre-cancerous with high accuracy in previous studies from our lab. By comparing the performance of virus models to bacterial models, we can begin to understand which microbial communities are most highly associated with cancer. That high association will suggest the microbial community is more involved in cancer biologically, although it is an admitedly tenuous link that will inform further studies. A better algorithm that incorporates virus information will also set a precedence for using virus signatures in predictive models, and an improved model would have therapeutic benefits.

*Hypothesis*: A virome-based classification model is more accurate than the bacterial alternative, and the two together will outperform an individual model.

*Approach*: I will use two apporaches to constructing classification models from the dataset. I will use the existing random forest protocol as well as generating and evaluating a set of classification models using the [Superlearner](https://cran.r-project.org/web/packages/SuperLearner/vignettes/SuperLearnerPresent.pdf) CRAN package (evaluated many approaches to dataset classification). I will compare the performance of the virus and bacterial models, as wel as compare those two to a model containing both datasets.

###Other Points
I also want to run through this data using my new network modeling approach, and compare community dynamics between the disease states. This might be worth also including in the predictive models.
