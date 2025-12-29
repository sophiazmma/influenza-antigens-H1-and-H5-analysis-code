# Influenza Antigens H1 and H5 Analysis

This repository contains the R code used for the analysis of single-cell RNA sequencing (scRNA-seq) and B cell receptor (BCR) repertoire data from mice infected with H1 and H5 influenza viruses. The analysis covers data preprocessing, clustering, differential expression analysis, and BCR lineage reconstruction.

## Overview

The study aims to characterize the B cell response to different influenza hemagglutinin antigens (H1 and H5). The code provided here allows for the reproduction of the key findings and figures presented in the manuscript.

## Repository Structure

- `analysis.R`: The main analysis script containing the entire workflow.
- `install_dependencies.R`: Script to install all required R packages.
- `data/`: Directory structure for input data (see **Data Availability**).
- `results/`: Directory for output figures and tables.

## Dependencies

The analysis relies on a variety of R packages and external tools.

### System Requirements
- IgBLAST (NCBI) 1.22.0 or compatible
- Build tools (e.g., `build-essential`, `libcurl4-openssl-dev`)

### R Packages
To install the required R packages, run the `install_dependencies.R` script:

```r
source("install_dependencies.R")
```

The key packages include:
*   **CRAN**: `Seurat`, `ggplot2`, `dplyr`, `tidyverse`, `patchwork`, `ggpubr`, `pheatmap`, `RColorBrewer`, `stringr`, `readr`, `tibble`, `future`, `clustree`, `eulerr`, `ggforce`, `harmony`, `vegan`, `circlize`, `treemap`.
*   **Bioconductor**: `BiocManager`, `msigdbr`, `clusterProfiler`, `GSVA`, `AUCell`, `monocle`, `SingleR`, `celldex`, `Biostrings`, `alakazam`, `shazam`, `tigger`, `dowser`, `ggtree`, `treeio`.
*   **GitHub**: `DoubletFinder`, `copykat`, `SCENIC`, `monocle3`, `scRepertoire`.

## Workflow Description

The `analysis.R` script is organized into the following major sections:

1.  **Environment Setup**: Loading libraries and setting up parallel processing.
2.  **BCR Annotation**: Instructions for running local IgBLAST on raw FASTA files to generate AIRR-compliant formatting.
3.  **Clonal Thresholding**: Using `shazam` to calculate the nearest neighbor distance threshold for defining B cell clones based on junction region similarity.
4.  **scRNA-seq Processing**:
    *   Loading 10x Genomics data.
    *   Quality control and filtering.
    *   Dimensionality reduction (PCA, UMAP) and clustering using `Seurat`.
5.  **Differential Expression Analysis**: Identifying marker genes for different B cell subpopulations (e.g., Naive, GC, Plasma cells) and comparing H1 vs. H5 responses.
6.  **BCR Repertoire Analysis**:
    *   Lineage tree reconstruction using `dowser` and `IgPhyML`.
    *   Visualization of phylogenetic trees with isotype and SHM annotations.
    *   Analysis of CDR3 length distributions and V-gene usage.

## Usage

1.  **Clone the repository**:
    ```bash
    git clone https://github.com/sophiazmma/influenza-antigens-H1-and-H5-analysis-code.git
    ```
2.  **Prepare Data**: Place your raw data files (Cell Ranger outputs, FASTA files) in the `data/` directory. *Note: You may need to adjust file paths in `analysis.R` to match your local directory structure if it differs.*
3.  **Run Analysis**: Open `analysis.R` in RStudio or run from the command line:
    ```bash
    Rscript analysis.R
    ```

## Data Availability

[Insert specific information about where the raw sequencing data is deposited, e.g., NCBI GEO Accession Number]

## License

[Insert License Information]
