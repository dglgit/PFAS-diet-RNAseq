# PFAS diet RNASeq

To learn R I am re-analyzing the RNASeq data from https://www.mdpi.com/2072-6643/13/11/3902
In the paper, the authors feed mice (Mus musculus) diets with high levels of PFOS (perfluorooctanoic acid).
I got the raw reads from GSE185183. 

I assembled the raw reads using the Galaxy platform. I used HISAT2 to align raw reads. I Gffcompare on stringtie-assembled transcriptomes and the Mus musculus annotated transcriptome (mm39) to generate the reference transcriptome. Then I used featurecounts to get counts of each transcript. The folders `pfos-mm/` and `control-mm/` contain the outputs of featurecounts for the PFOS-diet group and the cotnrol group, respectively. 

Due to github's size limitations, I cannot upload the reference .gtf file I used. Go to https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/ to download the mouse reference genome (refGene.gtf).

I also tried re-running the program using RNA Star feature counts run on Galaxy and a human annotated genome. 

`tryLoading.r` concatenates all featurecounts files and performs DESeq2. 
`process-star.r` processes and concatenates the featurecounts from RNA star and performs DESeq2.
`analyzeDS.r` runs fgsea on the output of DESeq2 to get pathways that are up- and down-regulated. 

