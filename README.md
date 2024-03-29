## Code repository for "RNA-seq Analysis of an Acute Myeloid Leukaemia Cell Line Treated with Phenolic Compounds"
Acute myeloid leukaemia (AML) is an aggressive cancer of the hematopoietic system, characterised by poorly differentiated myeloid cells. A preliminary, unpublished study extracted a phenolic compound from extra virgin olive oil, which was used to treat three samples of HL60 (an AML cell-line) cells for 1, 6 and 12 hours against a negative control. The total RNA of each sample was extracted and subcontracted to EMBL for sequencing, which returned a single-end read FASTQ file per sample. A bulk RNA-seq pipeline was developed to perform QC checks on the data files at various steps of the analysis, align the reads to a human reference genome, conduct differential gene expression analysis, and Gene Set Enrichment Analysis (GSEA). The flowchart below is a summary of the pipeline.

<img width="360" alt="Dissertation_pipeline" src="https://user-images.githubusercontent.com/74971601/194959190-9f29c72c-b205-42af-b736-9400c5bb08bd.png">

For reproducibility purposes the versions of the tools used have been provided:

| Tool          | Version |
|---------------|---------|
| FastQC        | 0.11.9  |
| FastQ Screen  | 0.14.0  |
| MultiQC       | 1.11    |
| Cutadapt      | 3.1     |
| Trim Galore!  | 0.6.7   |
| STAR          | 2.7.9a  |
| RSEM          | 1.3.3   |
| Bioconductor  | 3.13    |
| edgeR         | 3.14    |
| GAGE          | 3.15    |
| Pathview      | 1.30.1  |
| org.Hs.eg.db  | 3.1.0   |
| AnnotationDbi | 1.58.0  |
| devtools      | 2.4.3   |


