# Introduction
Acute myeloid leukaemia (AML) is an aggressive cancer of the hematopoietic system, characterised by poorly differentiated myeloid cells. A preliminary, unpublished study extracted a phenolic compound from extra virgin olive oil, which was used to treat three samples of HL60 (an AML cell-line) for 1, 6 and 12 hours against a negative control. The total RNA of each sample was extracted and subcontracted to EMBL for sequencing, which returned a single-end read FASTQ file per sample. A bulk RNA-seq pipeline was constructed to perform QC checks on the data files at various steps of the analysis, align the reads to a human reference genome, normalise the data, conduct Differential Gene Expression (DGE) analysis, and at a higher level, Gene Set Enrichment Analysis (GSEA).


<img width="342" alt="RNA-seq pipeline" src="https://user-images.githubusercontent.com/74971601/184557197-db66695b-48de-472a-b88b-d1506029575f.png">
