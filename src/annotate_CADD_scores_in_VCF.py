#! /usr/bin/env python

# Script Description:
# This script is used to annoate CADD scores in VCF Files based on Tabix results.

# Import necessary libraries
import os
import re
import sys

# Command-line arguments for input and output files
CADD_annotation = sys.argv[1]  # Path to the CADD annotation file
unannotated_vcf = sys.argv[2]  # Path to the VCF file to be annotated
annotated_vcf = sys.argv[3]  # Output path for the annotated VCF file

# Initialize dictionaries to store CADD scores
dict_RawScore = {}  # Stores raw scores
dict_PHRED = {}  # Stores PHRED scores

# Check if the CADD annotation file exists
if os.path.isfile(CADD_annotation):
    # Open and read the CADD annotation file
    with open(CADD_annotation, 'r') as f_CADD_anno:
        lines = f_CADD_anno.readlines()
        for line in lines:
            tags = line.strip().split('\t')
            Chrom_Pos_Ref_Alt = '_'.join(tags[0:4])  # Concatenate key identifiers
            dict_RawScore[Chrom_Pos_Ref_Alt] = tags[-2]  # Raw score
            dict_PHRED[Chrom_Pos_Ref_Alt] = tags[-1]  # PHRED score
else:
    print('File not exist: '+CADD_annotation)
    sys.exit(1)  # Exit if file not found

# Process the unannotated VCF file and add annotations
with open(unannotated_vcf, 'r') as f_in_vcf, open(annotated_vcf, 'w') as f_out_vcf:
    lines = f_in_vcf.readlines()
    # Identify the last INFO line to insert new INFO definitions correctly
    last_info_line = ''
    for line in lines:
        if line.startswith('##INFO='):
            last_info_line = line
    for line in lines:
        if line.startswith('#'):
            f_out_vcf.write(line)
            if line == last_info_line:
                # Add new INFO field definitions for CADD scores
                f_out_vcf.write('##INFO=<ID=CADD_RawScore,Number=1,Type=Float,Description="CADD score: RawScore">\n')
                f_out_vcf.write('##INFO=<ID=CADD_PHRED,Number=1,Type=Float,Description="CADD score: PHRED">\n')
        else:
            tags = line.strip().split('\t')
            alt_idx = int(re.split(r'[\/|]', tags[-1].split(':')[0])[1]) - 1
            Alt = tags[4].split(',')[alt_idx]
            # Adjust Chrom_Pos_Ref_Alt based on the presence of 'chr' prefix
            Chrom_Pos_Ref_Alt = '_'.join([tags[0].lstrip('chr'), tags[1], tags[3], Alt])
            RawScore = dict_RawScore.get(Chrom_Pos_Ref_Alt, '.')  # Retrieve RawScore, default to '.'
            PHRED = dict_PHRED.get(Chrom_Pos_Ref_Alt, '.')  # Retrieve PHRED, default to '.'
            # Append the CADD scores to the INFO column
            tags[7] += f';CADD_RawScore={RawScore};CADD_PHRED={PHRED}'
            f_out_vcf.write('\t'.join(tags) + '\n')

