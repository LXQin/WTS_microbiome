# Leveraging Host Whole Transcriptome Sequencing Data for Microbiome Studies: A Comparative Analysis with 16S rRNA Sequencing and Pre-processing Evaluation

This repository contains the data and scripts (.sh, .rmd) to reproduce the analysis in the original paper. To reproduce the analysis, following the steps below.

## Data access

The 16S rRNA sequencing dataset is publicly available and can be accessed through R package HMP2Data. The Whole transcriptomics data can be  downloaded as controlled-access sequences from the Genotypes and Phenotypes database (accession number: phs001626). 

## Data preparation

The script to get whole transcriptome sequencing (WTS) derived abundance table and the corresponding summary statistics are in /FASTQ_Processing/. The WTS-derived abundance table and the corresponding meta data should be in /raw data/ for later processing. The WTS-derived abundance table used in the manuscript is available upon request. 

## Run Rmarkdown code

- Run P1_Data_preprocessing.Rmd to pre-process the abundance tables and get Raw, Matched, Filtered, and Normalized datasets

- Run P2_Data_structure_visualization.Rmd to get descriptive statistics comparison of abundance

- Run P3_Diversity.Rmd to compute alpha and beta diversity and get visualization results for comparison.

- Run P4_DA.Rmd to get differential abundance analysis (DAA) results. Note you need to specify the contrast. And make sure you have the /DAResults/ folder to save the results for visualization.

- Run P5_DA_visualization.Rmd to visualize the DAA results for comparison.

## Folder structure for reproducing the results

raw data/

P1_Data_preprocessing.Rmd

P2_Data_structure_visualization.Rmd

P3_Diversity.Rmd

P4_DA.Rmd

P5_DA_visualization.Rmd

DAResults/

 └──CDUCvsNonIBD
 
 └──IBDvsNonIBD
 
