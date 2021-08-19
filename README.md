# Introduction
Acute myeloid leukaemia (AML) is an aggressive cancer of the hematopoietic system, characterised by poorly differentiated myeloid cells. A preliminary, unpublished study extracted a phenolic compound from extra virgin olive oil, which was used to treat three samples of an AML cell-line (HL60) for 1, 6 and 12 hours against a negative control. The total RNA of each sample was extracted and subcontracted to EMBL for sequencing, which returned a single-end read FASTQ file per sample. Commonly used tools for bulk RNA-Sequencing (RNA-Seq) will be used to perform QC checks on the data files at various steps of the analysis, align the reads to a human reference genome, normalise the data, and conduct differential expression analysis. After acquiring a list of genes which are up or down regulated after treatment, Gene Ontology (GO) terms may be assigned to them. Differential Gene Expression (DGE) analysis, mostly consisting of custom scripts will then be used to visualise and determine the effects of the treatment on gene expression.

# Sample information
- Single-end reads
- Sanger/Illumina 1.9 encoding (ASCII+33 quality scores as Phred scores)
- 51bp reads
- 30x depth
- HiSeq 2000 Sequencing System (Illumina) used by EMBL in Heidelberg, Germany
- Total RNA was extracted from the treated cells using the RNeasy Mini kit (Qiagen) procedure

| File name                                                             | Treatment         | Size (GB) | Sample name |
|-----------------------------------------------------------------------|-------------------|-----------|-----------|
| C9LMRACXX_5_16s004919-1-1_Zammit-Mangion_lane216s004919_sequence.txt  | 1 hour phenol  | 6.3      | A      |
| C9LMRACXX_10_16s004924-1-1_Zammit-Mangion_lane316s004924_sequence.txt | 6 hour phenol  | 5.2      | B      |
| C9LMRACXX_15_16s004929-1-1_Zammit-Mangion_lane416s004929_sequence.txt | 12 hour phenol | 4.4      | C      | 
| C9LMRACXX_16_16s004930-1-1_Zammit-Mangion_lane416s004930_sequence.txt | Negative control  | 5.0      | D      |
