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
CanSys requires no installation. However, users must execute `bash download.sh` to download the required databases prior to running the CanSys tool. Additionally, users need to specify the ANNOVAR directory path by modifying the `annovar_Dir` variable in the `run_CanSys.sh` script.


## Usage
./run_pipeline.sh [OPTIONS]

Options:

 -Vcf <path>: (Required) Path to the VCF file containing somatic or germline variants for annotation and analysis.
 
 -SampleName <name>: (Required) Name of the sample being analyzed. Used for naming output files.
 
 -Expression <path>: (Optional) Path to the gene expression file. Read count data is recommended. Plasea refer to `test/input/expression.txt` for an example.
 
 -Database <name>: (Required) Specifies the database for pathway-level score calculation. Valid options are GO, KEGG, or ALL.

 -Cancer <name>: (Required) Specifies the cancer type. Valid options are Ampulla_of_Vater, Biliary_Tract, Bladder_or_Urinary_Tract, Bone, Bowel, Breast, Cervix, CNS_or_Brain, Esophagus_or_Stomach, Eye, Head_and_Neck, Kidney, Liver, Lung, Lymphoid, Myeloid, Ovary_or_Fallopian_Tube, Pancreas, Peripheral_Nervous_System, Pleura, Prostate, Skin, Soft_Tissue, Testis, Thyroid, Uterus and Vulva_or_Vagina.
 
 -Cutoff_CADD <numeric>: (Optional) Cutoff value for CADD scores. Default: 0.
 
 -nPermSimple <numeric>: (Optional) Number of permutations used in the permutation test to estimate P-values. Default: 1000.

Examples:
#Run the CanSys tool using the somatic VCF file. Not filter out unexpressed genes.
./run_pipeline.sh -Vcf ./test/input_files/somatic.vcf -SampleName example -Database ALL -Cancer Breast
#Run the CanSys tool using the somatic VCF file. Filter out unexpressed genes. Set the number of permutations to 6000.
./run_pipeline.sh -Vcf ./test/input_files/somatic.vcf -SampleName example -Expression ./test/input_files/expression.txt -Database ALL -Cancer Breast -nPermSimple 8000
#Run the CanSys tool using the germline VCF file. Filter out unexpressed genes. Filter out variants with CADD scores less then 15.
./run_pipeline.sh -Vcf ./test/input_files/germline.vcf -SampleName example -Expression ./test/input_files/expression.txt -Database ALL -Cancer Breast -Cutoff_CADD 15

 

