# Genetic Susceptibility to Schizophrenia Through Neuroinflammatory Pathways is Associated with Retinal Thickness: Findings from the UK Biobank

## Author
Finn Rabe <finn dot rabe at bli dot uzh dot ch>
Philipp Homan<philipp dot homan at bli dot uzh dot ch>

## Overview
This repository contains the code, data, and supplementary materials for the study titled "Genetic Susceptibility to Schizophrenia Through Neuroinflammatory Pathways is Associated with Retinal Thickness: Findings from the UK Biobank". The study investigates the relationship between genetic predisposition to schizophrenia, measured via polygenic risk scores (PRS), and retinal thickness in healthy individuals. The findings suggest that genetic susceptibility to schizophrenia, particularly through neuroinflammatory pathways, is associated with retinal thinning, providing insights into potential early biomarkers for schizophrenia.

## Repository Structure
1. src/
retinflam_load.R: Script to load and preprocess data.
retinflam_do.R: Main analysis script performing statistical computations, including robust regression and mediation analyses.
retinflam_do.py: Python script for generating filtered datasets and figures.
2. output/figures
figures/: Generated figures and visualizations.
Heatmaps of retinal subfield associations.
Mediation analysis tables in PDF format.
3. pub/
References (references.bib): Bibliography file for citations in the manuscript.
4. Manuscript Files
retinflam_ms.Rmd: R Markdown file for generating the manuscript with embedded code and results.

## Data Sources
Data used in this study were obtained from the UK Biobank (application ID: 102266).
Genetic data were processed following QC protocols to ensure high-quality SNPs and imputed genotypes.
Retinal imaging data were acquired using optical coherence tomography (OCT).

### Prerequisites
- R (version 4.0.5 or later)
- Python (version 3.9 or later)

## Installing
Clone the repository: **git clone https://github.com/homanlab/prsoct.git**

## Producing all analyses and manuscript
**Rscript -e "rmarkdown::render('retinflam_ms.Rmd')"**
