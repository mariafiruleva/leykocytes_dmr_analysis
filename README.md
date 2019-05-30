# Effect of smoking on human leukocyte epigenome
__Students:__
Polina Pavlova, Maria Firuleva

__Supervisors:__
Oleg Sergeev, Yulia Medvedeva, Julia Kornienko

# Project Description

The Russian Children’s Study is a prospective cohort of 516 boys who were enrolled at 8–9 years of age and provided semen samples at 18–19 years of age. RRBS of sperm was conducted to identify the methylation level of CpG dinucleotides. At the moment of enrollment into the study, the TCDD dioxine (which is one of the most harmful endocrine disrupting chemicals) concentration in the blood of each boy was measured for further longitudinal study of its influence on the reproductive health. Moreover, each boy visited the clinic biennially - for blood sampling; annually - for urine sampling, follow up of growth and puberty and interviewing => 20 000+ sample aliquots and 1000+ analyzing parameters were collected in total for further analysis.

##### What is already known?
- 52 differentially methylated regions (DMRs) in sperm were identified that distinguished lowest and highest peripubertal serum TCDD concentrations (RCS, Pilsner et al., 2018)
- Peripubertal exposure to toxicants (like dioxins) and smoking affects the methylation of sperm DNA at the age of 18 (last year project, Y.Kornienko)

##### What is needed to be known?
- How does smoking influence the DNA methylation in peripheral blood leukocytes at the age of 18? 
- How are associated epigenome changes of sperm and leukocytes in the same study participants?

##### Aims of the project:
- Analysis of smoking influence on the DNA methylation level of peripheral blood leukocytes at the age of 18
- Comparison of DNA methylation level in sperm (results of last year project, Julia Kornienko) and leukocytes in the same study participants

##### What were the data to analyse:
- 36 samples: 11 smoking and 25 don’t
- Data regarding lifestyle habits of each of 36 chosen participants
- Methylation levels of 2277623 CpGs present in at least one of 36 samples

##### Scripts Description
1. data_preparation.ipynb - script for data preprocessing 
2. visualization.ipynb - script for drawing plots presented here
3. Aclust_GEE.R - script for A-clustering and GEE analysis
4. DMRcate_RRBS_leuk.R - script for DMRcate analysis

##### Brief summary of the results:
1. Similar methylation distribution revealed in leukocytes and sperms - for dataset with CpGs presented at least sample
2. Methylation distribution is different in sperms and leukocytes in dataset used for A-clustering - CpGs presented in all samples
3. Comparison of CpGs distribution across samples shows similar distribution in sperms and leukocytes
4. Comparison of methylation range also shows similar distibution
5. Using A-clustering method found 217 clusters. After GEE identified 77 significant clusters (p-value < 0.05)
6. A-clustering+GEE comparison shows similar distribution in sperms and leukocytes. In sperms data were identified 814 clusters (A-clust) and 136 significant clusters after GEE. 
7. DMRcate package (R) were used for search of DMR by smoking. Factor - smoking last 6 months in binary classification (yes/no). 34 significant DMRs revealed, of them 19 overlap with at least one promoter (reference - hg38), 145 significant CpGs. 
8. Genes associated with significant CpGs encode pseudogenes, lincRNAs, antisense RNAs, zinc fingers, miRNA, regulatory elements, proteins of cell adhesion, amino-acid transporter.

##### Detailed description of results:
*Samples selection:*
- 51 out of 516 samples were chosen by the following criterias: 
- Prepubertal TCDD exposure at 8-9 years old
- Semen quality at 18-19 years old 
- Frozen semen samples at 18-19 years old 
- Buffy coat samples at 18-19 years old closest to date of semen samplingю Buffy coat samples – source of leukocyte libraries.

To start with, 10 IDs with either the highest or the lowest TCDD concentration were selected. Then were selected one ID with the highest semen quality and one with the lowest. Later, 39 more IDs were chosen by random selection 13 IDs belonging to each tercile in terms of TCDD concentration. Therefore the data set of 51 samples was created. 
45 leukocyte libraries were sequenced. 4 samples were excluded because of low quality. Then additionally 5 samples were excluded because of coverage less than 10. Finally, 36 samples with satisfactory sequencing quality were selected for this project .

**_Data preprocessing:_** 

From the file containing information regarding all CpGs for 41 samples (with CpGs coverage >= 10x) we extracted methylation levels for the selected 36 samples. So, we got dataset with methylation level of 2277623 CpGs for 36 samples (dataset with information of CpGs presented in at least one sample). Then we prepared dataset for A-clustering algorithm (it works with methylation levels of CpGs presented in all samples). So, for A-clustering+GEE approach we had dataset with methylation level of 19466 CpGs. 

**_A-clustering+GEE:_**

Next, we made A-clustering and GEE regression using R packages “Aclust” and “gee”. This algorithm allows to detect sets of neighboring CpGs sites that are correlated with each other. With this approach we identified 217 A-clusters, from them 77 were significant (p-value < 0.05). But this approach has a significant drawback. Because of this algorithm works with information presented in all samples (without any NAs), we needed to restrict data from 2277623 to 19466 CpGs. So, a lot of important information could have been lost. 

Thus, another approach were suggested.  

_**DMRcate:**_

DMRcate - R package for search of DMRs associated with exposure to a factor. In our case, we used smoking in binary classification (smoke or not) to find DMRs associated with smoking influence. 145 significant CpGs and 34 significant DMRs (p-value < 0.05) were found in our data. From them 19 DMRs overlap with at least one promoter (reference - hg38). We found 23 genes associated with significant DMRs. These genes are associated with antisense RNAs, lincRNAs, miRNA, pseudogenes, zinc fingers, transcriptional factor, spliceosome, cell adhesion and migration, kinase, metalloprotease, electron transport chain and amino-acid transporter. 

_**Comparison with sperms**_ (last year project - Julia Kornienko):

- Comparison of methylation distribution for datasets with methylation level for CpGs presented in at least one sample. We can see similar methylation distribution - the most of CpGs have very low methylation level, and some CpGs have high methylation level

![GitHub Logo](/plots/violinplot_leuk_all_CpGs.png)
![GitHub Logo](/plots/violinplot_sperm_all_CpGs.png)

- Comparison of methylation distribution for dataset with methylation level for CpGs presented in all samples. Here we can see different distribution. In sperms the most of CpGs have very low level of methylation and some CpGs with high methylation, but in leukocytes vise versa. We suppose, the reason is that global demethylation occurs in sex cells

![GitHub Logo](/plots/violinplot_leuk_dataset_for_Aclust.png)
![GitHub Logo](/plots/violinplot_sperms_dataset_for_Aclust.png)

- Comparison of CpGs distribution across sample. In both type of cells we can see geometric distribution. Most often the CpGs is found in one sample, 1.5 times less in two samples and etc

![GitHub Logo](/plots/number_samples.png)
![GitHub yukornienko/mediation_analysis](Figure4.png)

- A-clustering+GEE approach comparison. In general, we see similar distribution of number of CpGs in cluster. In leukocytes we found 217 clusters, in sperms - 874. Number of significant clusters (p-value < 0.05) in leukocytes 77, in sperms - 136
