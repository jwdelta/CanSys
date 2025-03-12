# CanSys
CanSys is tool that quantifies the contributions of germline and somatic alterations to biological pathway-level disturbances in individual cancer samples. It annotates Variant Call Format (VCF) files with Combined Annotation Dependent Depletion (CADD) scores, calculates gene-level impact scores using CADD and Cancer Dependency Map (DepMap) scores, and maps affected genes onto biological pathways using specified databases (GO and KEGG). We also offer a web-based application [cansysplot](https://cansysplot.com/) that allows users to analyze and interactively visualize pathway disturbances in networks. Additionally, the application facilitates the visualization of disturbance networks using data derived from The Cancer Genome Atlas (TCGA).

## Requirements
Before running this script, ensure the following software packages and dependencies are installed in your environment.

**Software Requirements:**
 - Python (version 3): Required for running Python scripts involved in the pipeline.
 - R (version 4): Required for executing R scripts for statistical analysis.
 - tabix: Essential for indexing and rapid retrieval of compressed TAB-delimited genomic data files.

**Python Packages:**
- pandas: Provides robust data structures for efficient data manipulation and analysis.

**R Packages:**
- ANNOVAR: Utilized for annotating which genes are associated with variations. Please also download the gene-based annotation file "refGene" (build: hg38).
- rGMMtest: Utilized for fitting Gaussian Mixture Model (GMM) to 1D data and assigning data points to individual Gaussian components. (This R package is accompanied by a manuscript currently under submission)
- optparse: Facilitates the parsing of command-line options in R scripts.
- fgsea: Fast Gene Set Enrichment Analysis tool.

## Setup
CanSys requires no installation. However, users must execute `bash download.sh` to download the required databases prior to running the CanSys tool. This process may take ～2 hours to complete, and requires ～83GB of available storage. Additionally, users need to specify the ANNOVAR directory path by modifying the `annovar_Dir` variable in the `run_CanSys.sh` script.

## Usage
```shell
./run_CanSys.sh [OPTIONS]
```

**Options**:

 -Vcf <path>: (Required) Path to the VCF file containing somatic or germline variants for annotation and analysis.
 
 -SampleName <name>: (Required) Name of the sample being analyzed. Used for naming output files.
 
 -Expression <path>: (Optional) Path to the gene expression file. Read count data is recommended. Plasea refer to `test/input/expression.txt` for an example.

 -Output_Dir <path>: (Optional) Specifies the output directory for saving results. If not specified, files will be saved in the current working directory by default.
 
 -Database <name>: (Required) Specifies the database for pathway-level score calculation. Valid options are GO, KEGG, or ALL.

 -Cancer <name>: (Required) Specifies the cancer type. Valid options are Ampulla_of_Vater, Biliary_Tract, Bladder_or_Urinary_Tract, Bone, Bowel, Breast, Cervix, CNS_or_Brain, Esophagus_or_Stomach, Eye, Head_and_Neck, Kidney, Liver, Lung, Lymphoid, Myeloid, Ovary_or_Fallopian_Tube, Pancreas, Peripheral_Nervous_System, Pleura, Prostate, Skin, Soft_Tissue, Testis, Thyroid, Uterus and Vulva_or_Vagina.
 
 -Cutoff_CADD <numeric>: (Optional) Cutoff value for CADD scores. Default: 0.
 
 -nPermSimple <numeric>: (Optional) Number of permutations used in the permutation test to estimate P-values. Default: 1000.

**Examples**:
```shell
#Run the CanSys tool using the somatic VCF file. Not filter out unexpressed genes.
./run_CanSys.sh -Vcf /path/to/test/input_files/somatic.vcf -SampleName example -Output_Dir /path/to/test/output_files -Database ALL -Cancer Breast

#Run the CanSys tool using the somatic VCF file. Filter out unexpressed genes. Set the number of permutations to 6000.
./run_CanSys.sh -Vcf /path/to/test/input_files/somatic.vcf -SampleName example -Expression /path/to/test/input_files/expression.txt -Output_Dir /path/to/test/output_files -Database ALL -Cancer Breast -nPermSimple 6000

#Run the CanSys tool using the germline VCF file. Filter out unexpressed genes. Filter out variants with CADD scores less then 15.
./run_CanSys.sh -Vcf /path/to/test/input_files/germline.vcf -SampleName example -Expression /path/to/test/input_files/expression.txt -Output_Dir /path/to/test/output_files -Database ALL -Cancer Breast -Cutoff_CADD 15
```

## Outputs
Two output files will be generated for both GO and KEGG analyses: one containing the complete results and the other including only the significant results with Benjamini-Hochberg (BH) adjusted P-values < 0.25. The output files contain detailed information across several columns, each providing specific insights into the biological pathways analyzed. Here is a brief overview of the columns included in the output files:
 - goid/entry: This column contains the ID of the biological pathway.
 - Name: Lists the names of the biological pathways.
 - N: Represents the number of genes within each biological pathway.
 - DE: Indicates the number of affected genes present in each biological pathway.
 - ES: Refers to the enrichment score.
 - NES: Refers to the normalized enrichment score.
 - GeneScoreMean: Refers to the average gene-level impact scores in each biological pathway. 
 - PDS: Refers to the pathway disturbance score (CanSys score).
 - pval: Contains the P-values.
 - padj: Contains the P-values adjusted using the Benjamini-Hochberg (BH) method.
 - Affected Genes: Lists the genes that have a gene-level impact score greater than 0 in each biological pathway.

## Citing this work
If you use the CanSys tool or its web-based application in your research, please cite: Common and rare germline variants together with somatic mutations alter the integrity of cancer hallmark regulatory networks. (Currently under submission)

## Acknowledgments
We would like to express our sincere gratitude to the developers of CADD, DepMap, and all other algorithms and dependencies integrated into our tool.

## Contributing
We welcome any bug reports, enhancement requests, and other contributions. To submit a bug report or enhancement request, please use the [GitHub issues tracker](https://github.com/JiaweiDai-create/CanSys/issues).

## License
This project is licensed under the terms of the MIT license.
