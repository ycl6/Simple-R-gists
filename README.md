# Simple R gists

### Showcase short R examples

-----

**License:** GPL-3.0

## National survey

#### - [Deaths registered monthly in England and Wales](ONS_Monthly_Deaths_Data_England_and_Wales)

This notebook shows how to use `readxl` package to retreive the *Deaths registered monthly in 
England and Wales* Dataset from [Office for National Statistics](https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/monthlyfiguresondeathsregisteredbyareaofusualresidence) 
and create various plots to show the number of deaths with `ggplot2`.

#### - [COVID-19 Vaccination Statistics in England](NHS_England_COVID-19_Vaccination)

This notebook shows how to use `readxl` package to retreive the *Monthly COVID-19 Vaccinations* Dataset 
from [NHS England](https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-vaccinations/) 
and create various plots to show key statistics with `ggplot2`.

## Gene expression

#### - [Visualise differential expression analysis results with MA plot](RNA-seq_DEA_MA-plot)

This notebook shows readers how to use ggplot2 package to create various plots to show the 
magnitude of the change in expression between two groups (M) in comparison to the magnitude 
of average expression over all samples (A).

#### - [Normalisation of bulk RNA-seq counts](RNA-seq_Count_Normalisation)

This is a reproducible demo that shows readers how to use perform count normalisation of 
RNA-seq data using several common methods:

- DESeq2’s median-of-ratios
- edgeR’s TMM (trimmed mean of M-values)
- RPM/FPM/CPM (reads/fragments/counts per million mapped reads/fragments)
- TPM (transcripts per million mapped reads)
- RPKM/FPKM (reads/fragments per kilobase of transcript length per million mapped reads/fragments)

#### - [Plot ligand-receptor expression on UMAP and t-SNE projection](plotReducedDimLR)

This is a reproducible demo that shows readers how to use the function `plotReducedDimLR` which 
make use of the [ggnewscale](https://cran.r-project.org/package=ggnewscale) R package to ovelay 
the expressions of ligand and receptor genes on a single plot of UMAP or t-SNE projection.

## Others

#### - [Visualise interaction between 2 variables with parallel sets diagrams](Plot_Parallel_Sets_Diagrams)

This notebook shows readers how to use `ggforce` package to create parallel sets diagram to 
show the interaction between 2 categorical variables.

#### - [Testing of the Spatial Image Analysis of Tissues (SPIAT) R package](Test_SPIAT)

This notebook shows how to use SPIAT to perform spatial data processing, quality control, 
visualization and analysis.
