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

Sequencing library preparations were done using the Illumina **NexteraXT** prep kit. The samples will be sequenced on the HiSeq 4000 with one lane devoted to the virome library and the other devoted to the whole metagenome library.

###Analytical Questions
1. **What specific gut viruses are associated with colon cancer?** We know that viruses have causative roles in a variety of cancers, and this may be the case in colon cancer. Here we will evaluate what viruses are in higher abundance in colon cancer cases and may therefore play a role in cancer developement.

2. **What virus community signatures are associated with colon cancer, and how do these compare to the bacterial communities?** We know that there are bacterial community signatures associated with colon cancer developement and presence. What is the virome diversity associated with colon cancer? It will also be important to assess how these trends compare to the bacterial community signatures. In many cases, lower bacterial diversity has been associated with higher viral diversity in disease states.

3. **Can we construct an accurate predictive model using virome data, and can this supplement models using bacterial data?** We can construct models for predicting and classifying colon cancer based on bacterial community information. We will try this with the virus community data, compare it's accuracy to that of bacterial community models, and evaluate the extent to which they supplement each other.
