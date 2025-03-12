#!/bin/bash

# Script Name: run_CanSys.sh
# Description: Automates the quantification of biological pathway anomalies.
# Author: Jiawei Dai
# Data Created: 2024-2-5
# Last Modified: 2025-1-8

# ==================================================================
# Section: Input/Output Configuration
# Description: Specifies the paths for input files, output files, and execution parameters essential for the pipeline process.
# Note: Please ensure the variant positions are aligned to the hg38 reference genome.
# ==================================================================

# Initialize variables for input/output files and parameters
Vcf=""
SampleName=""
Expression=""
Output_Dir=""
Database=""
Cancer=""
Cutoff_CADD="0"
nPermSimple="1000"
while [[ $# -gt 0 ]]; do
    case "$1" in
	   -Vcf) Vcf="$2"; shift; shift ;;
	   -SampleName) SampleName="$2"; shift; shift ;;
	   -Expression) Expression="$2"; shift; shift ;;
	   -Output_Dir) Output_Dir="$2"; shift; shift ;;
	   -Database) Database="$2"; shift; shift ;;
	   -Cancer) Cancer="$2"; shift; shift ;;
	   -Cutoff_CADD) Cutoff_CADD="$2"; shift; shift ;;
	   -nPermSimple) nPermSimple="$2"; shift; shift ;;
	*) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Validate required input parameters
if [[ -z "$Vcf" ]]; then
    echo "Error: -Vcf is a required option."
    exit 1
fi

if [[ -z "$SampleName" ]]; then
    echo "Error: -SampleName is a required option."
    exit 1
fi

if [[ -z "$Database" ]]; then
    echo "Error: -Database is a required option."
    exit 1
fi

if [[ -z "$Cancer" ]]; then
    echo "Error: -Cancer is a required option."
    exit 1
fi

# Check if Database is one of the allowed values
if [[ ! "$Database" =~ ^(GO|KEGG|ALL)$ ]]; then
    echo "Error: -Database must be GO, KEGG, or ALL."
    exit 1
fi


# Specify the directory containing the pipeline scripts and resources
Pipeline_Dir=$(pwd)

# Specify the download directory for ANNOVAR
annovar_Dir=/path/to/ANNOVAR/annovar

# Resources
CADD_prescored=$Pipeline_Dir/resources/CADD_GRCh38_v1.7_no_anno
DepMap=$Pipeline_Dir/resources/DepMap_scores_23Q2/${Cancer}.txt

# Specifies the output directory for saving results
if [[ -z "$Output_Dir" ]]; then
    Output_Dir=$Pipeline_Dir
fi

# Create a temporary directory for execution
Temp_ExecDir="${Pipeline_Dir}/test/$(date +%s%N${RANDOM} | sha256sum | head -c 16).tmp"
echo "Temporary execution directory created at: ${Temp_ExecDir}"
mkdir -p "${Temp_ExecDir}" && cd "${Temp_ExecDir}" || exit


# ==================================================================
# Section: Variant Annotation
# Description: Annotates VCF files with CADD scores.
# ==================================================================

annotate_variants() {
    local vcf="$1"
    local sample="$2"
    echo "Annotating ${vcf} for sample ${sample}..."
    # Run ANNOVAR
    $annovar_Dir/table_annovar.pl $vcf $annovar_Dir/humandb/ -buildver hg38 -out $sample -remove -protocol refGene -operation g -nastring . -polish -vcfinput
    # Process annotation here
    grep -v '^#' "${sample}.hg38_multianno.vcf" | cut -f 1,2 | while read -r chr pos; do
    	if [[ "$chr" =~ ^chr ]]; then
    		chr=${chr:3}
    	fi
    	tabix "$CADD_prescored/whole_genome_SNVs.tsv.gz" "$chr:$pos-$pos" >> "${sample}.CADD_annotation.txt"
    	tabix "$CADD_prescored/gnomad.genomes.r4.0.indel.tsv.gz" "$chr:$pos-$pos" >> "${sample}.CADD_annotation.txt"
    done
    if [ -e "${sample}.CADD_annotation.txt" ];then
        python "$Pipeline_Dir/src/annotate_CADD_scores_in_VCF.py" "${sample}.CADD_annotation.txt" "${sample}.hg38_multianno.vcf" "${sample}.annotated.vcf"
    else
        cp "$vcf" "${sample}.annotated.vcf"
    fi
}

# Annotate VCFs
annotate_variants "$Vcf" "$SampleName"


# ==================================================================
# Section: Gene-Level Impact Score Calculation
# Description: Calculates gene-level impact scores using CADD and DepMap scores.
# ==================================================================

# Binary representation of gene expression status if Expression file is provided
if [[ -n "$Expression" ]]; then
	echo 'Initiating Binary Representation of Gene Expression Status.'
    Rscript --slave "$Pipeline_Dir/src/GeneExprBinaryRep.R" "$Expression" "${SampleName}.GeneExprBinaryRep.txt"
fi

# Calculation of Gene-Level Impact Score
echo "Calculating gene-level impact scores for ${SampleName}..."
cmd="python \"$Pipeline_Dir/src/cal_gene_level_scores.py\" --Vcf \"${SampleName}.annotated.vcf\" --SampleName \"$SampleName\" --DepMap \"$DepMap\" --cutoff_CADD \"$Cutoff_CADD\""

# Add optional --GeneExprBinaryRep if provided
if [[ -n "$Expression" && -f "${SampleName}.GeneExprBinaryRep.txt" ]]; then
    cmd+=" --GeneExprBinaryRep \"${SampleName}.GeneExprBinaryRep.txt\""
fi

# Execute the command
eval $cmd


# ==================================================================
# Section: Pathway-Level Score Calculation
# Description: Maps affected genes onto biological pathways
# ==================================================================

echo "Calculating pathway-level scores for ${SampleName}..."
if [[ "$Database" == "GO" ]]; then
    Rscript --slave "$Pipeline_Dir/src/cal_pathway_level_scores.R" --gene_level_impact_scores "${SampleName}.gene_level_impact_scores.txt" --database_name GO --database_resource "$Pipeline_Dir/resources/go_bp.rds" --universe_genes "$Pipeline_Dir/resources/universe.Homo_sapiens.rds" --nPermSimple "${nPermSimple}" --full_output_file_path "${SampleName}.pathway_level_scores.GO.txt" --signif_output_file_path "${SampleName}.signif.pathway_level_scores.GO.txt"
fi

if [[ "$Database" == "KEGG" ]]; then
    Rscript --slave "$Pipeline_Dir/src/cal_pathway_level_scores.R" --gene_level_impact_scores "${SampleName}.gene_level_impact_scores.txt" --database_name KEGG --database_resource "$Pipeline_Dir/resources/kegg_exclude_OS_and_HD_pathways.rds" --universe_genes "$Pipeline_Dir/resources/universe.Homo_sapiens.rds" --nPermSimple "${nPermSimple}" --full_output_file_path "${SampleName}.pathway_level_scores.KEGG.txt" --signif_output_file_path "${SampleName}.signif.pathway_level_scores.KEGG.txt"
fi

if [[ "$Database" == "ALL" ]]; then
    Rscript --slave "$Pipeline_Dir/src/cal_pathway_level_scores.R" --gene_level_impact_scores "${SampleName}.gene_level_impact_scores.txt" --database_name GO --database_resource "$Pipeline_Dir/resources/go_bp.rds" --universe_genes "$Pipeline_Dir/resources/universe.Homo_sapiens.rds" --nPermSimple "${nPermSimple}" --full_output_file_path "${SampleName}.pathway_level_scores.GO.txt" --signif_output_file_path "${SampleName}.signif.pathway_level_scores.GO.txt"
    Rscript --slave "$Pipeline_Dir/src/cal_pathway_level_scores.R" --gene_level_impact_scores "${SampleName}.gene_level_impact_scores.txt" --database_name KEGG --database_resource "$Pipeline_Dir/resources/kegg_exclude_OS_and_HD_pathways.rds" --universe_genes "$Pipeline_Dir/resources/universe.Homo_sapiens.rds" --nPermSimple "${nPermSimple}" --full_output_file_path "${SampleName}.pathway_level_scores.KEGG.txt" --signif_output_file_path "${SampleName}.signif.pathway_level_scores.KEGG.txt"
fi


# Check if GO pathway-level score file exists and move it to the outputs directory
mkdir -p Output_Dir
if [[ -e "${SampleName}.pathway_level_scores.GO.txt" ]]; then
    mv "${SampleName}.pathway_level_scores.GO.txt" "$Output_Dir"
    echo "GO pathway-level scores moved to outputs directory."
fi

if [[ -e "${SampleName}.signif.pathway_level_scores.GO.txt" ]]; then
    mv "${SampleName}.signif.pathway_level_scores.GO.txt" "$Output_Dir"
    echo "Significant GO pathway-level scores moved to outputs directory."
fi

# Check if KEGG pathway-level score file exists and move it to the outputs directory
if [[ -e "${SampleName}.pathway_level_scores.KEGG.txt" ]]; then
    mv "${SampleName}.pathway_level_scores.KEGG.txt" "$Output_Dir"
    echo "KEGG pathway-level scores moved to outputs directory."
fi

if [[ -e "${SampleName}.signif.pathway_level_scores.KEGG.txt" ]]; then
    mv "${SampleName}.signif.pathway_level_scores.KEGG.txt" "$Output_Dir"
    echo "Significant KEGG pathway-level scores moved to outputs directory."
fi

# Clean up: Remove the temporary execution directory to tidy up after script execution
/bin/rm -rf "${Temp_ExecDir}"
echo "Temporary directory ${Temp_ExecDir} has been removed."


echo "Pipeline execution completed."






