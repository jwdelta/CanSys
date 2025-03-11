# CanSys
CanSys is tool that quantifies the contributions of germline and somatic alterations to biological pathway-level disturbances in individual cancer samples. It annotates VCF files with CADD scores, calculates gene-level impact scores using CADD and DepMap scores, and maps affected genes onto biological pathways using specified databases (GO and KEGG).

## Requirements
Before running this script, ensure the following software packages and dependencies are installed in your environment.

Software Requirements:
 - Python (version 3): Required for running Python scripts involved in the pipeline.
 - R (version 4): Required for executing R scripts for statistical analysis.
 - tabix: Essential for indexing and rapid retrieval of compressed TAB-delimited genomic data files.

Python Packages:
- pandas: Provides robust data structures for efficient data manipulation and analysis.

R Packages:
- ANNOVAR: Utilized for annotating which genes are associated with variations. Please also download the gene-based annotation file "refGene" (build: hg38).
- rGMMtest: Utilized for fitting Gaussian Mixture Model (GMM) to 1D data and assigning data points to individual Gaussian components. (This R package is currently under submission.)
- optparse: Facilitates the parsing of command-line options in R scripts.
- fgsea: Fast Gene Set Enrichment Analysis tool.

## Setup
CanSys requires no installation. However, users need to execute `bash download.sh` to download the necessary databases before running the CanSys tool.

