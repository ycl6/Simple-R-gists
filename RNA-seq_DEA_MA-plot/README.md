# Visualise differential expression analysis results with MA plot

I-Hsuan Lin

University of Manchester

February 08, 2021

## Introduction
This notebook shows readers how to use ggplot2 package to create various plots to show the magnitude of the change in expression between two groups (M) in comparison to the magnitude of average expression over all samples (A).

## About this dataset

The RNA-seq dataset [E-MTAB-8411](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8411) contains 5 wild-type (GK1, GK3, GK5, GK7 and GK9) and 4 macrophage-specific Bmal1 knockout samples (GK2, GK4, GK6 and GK10). The fastq files were retreived from [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8411/samples/) and processed with the follow steps:

Adapter and QC trimming with BBDuk
Mapping to Mouse genome (assemble MM10 and [GENCODE M24](https://www.gencodegenes.org/mouse/release_M24.html) annotation) with STAR
Differential expression analysis with DESeq2. The log2 fold change shrinkage was not applied.
