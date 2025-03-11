#! /usr/bin/env python

import argparse
import pandas as pd

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Load the required input files.')
    parser.add_argument('-i', '--Vcf', required=True, help='Path to the VCF file containing somatic or germline variants.')
    parser.add_argument('-s', '--SampleName', required=True, help='Sample name.')
    parser.add_argument('-d', '--DepMap', required=True, help='Path to the DepMap file.')
    parser.add_argument('-g', '--GeneExprBinaryRep', help='Path to the gene expression status file (optional).')
    parser.add_argument('-c', '--cutoff_CADD', required=True, help='Cutoff value for CADD scores')
    return parser.parse_args()

def get_variant_level_AF_and_CADD(input_vcf, cutoff_CADD):
    """Extract variant-level allele frequency and CADD scores from a VCF file.
    
    Args:
        input_vcf (str): Path to the VCF file.

    Returns:
        dict: A dictionary with variant identifiers as keys and allele frequency, CADD scores as values.
    """

    dict_AF_and_CADD = {}
    with open(input_vcf, 'r') as f:
        for line in f:
            if line.startswith('#CHROM'):
                headers = line.strip().split('\t')
                INFO_index = headers.index("INFO")
                FORMAT_index = headers.index("FORMAT")
            elif not line.startswith('#'):
                parts = line.strip().split('\t')
                INFO = {k: v for k, v in (item.split('=') for item in parts[INFO_index].split(';') if '=' in item)}
                FORMAT_parts = parts[FORMAT_index+1].split(':')
                FORMAT_dict = {fmt: val for fmt, val in zip(parts[FORMAT_index].split(':'), FORMAT_parts)}
                AD_str = FORMAT_dict.get('AD', '0')
                DP_str = FORMAT_dict.get('DP', '0')
                AD_list = [int(x) for x in AD_str.split(',')]
                if len(AD_list) > 1:
                    DP = sum(AD_list)
                else:
                    DP = int(DP_str)
                AF = AD_list[-1] / DP if DP > 0 else 0
                gene = INFO.get('Gene.refGene', '').replace('\\x3b', ';').split(';')[0]
                if INFO['CADD_PHRED'] != '.':
                    CADD_PHRED = float(INFO.get('CADD_PHRED', '0')) / 100
                    if CADD_PHRED > (float(cutoff_CADD)/100):
                        key = '_'.join([gene, parts[0], parts[1], parts[3], parts[4]])
                        dict_AF_and_CADD.setdefault(key, []).extend([round(AF,4), CADD_PHRED])
    return dict_AF_and_CADD

def calculate_gene_level_impact_score(input_vcf, DepMap_file, cutoff_CADD, GeneExprBinaryRep_file=None):
    """Calculate gene-level impact score, considering optional files.

    Args:
        input_vcf (str): Path to the input VCF file.
        DepMap_file (str): Path to the DepMap file.
        GeneExprBinaryRep_file (str, optional): Path to the gene expression status file.

    Returns:
        dict: A dictionary with genes as keys and calculated impact scores as values.
    """

    # Initial processing for the input VCF
    dict_AF_and_CADD = get_variant_level_AF_and_CADD(input_vcf, cutoff_CADD)
    dict_variant_level_CADD = {v: scores[1] for v, scores in dict_AF_and_CADD.items()}
    df_variant_level_CADD = pd.DataFrame(list(dict_variant_level_CADD.items()), columns=['Variant', 'CADD_score'])
    df_variant_level_CADD['Gene'] = df_variant_level_CADD['Variant'].apply(lambda x: x.split('_')[0])
    
    # Load DepMap data and calculate impact score
    df_DepMap = pd.read_csv(DepMap_file, sep='\t')
    dict_DepMap = dict(zip(df_DepMap['Gene_Symbol'], df_DepMap['all_GeneDependency_avg']))
    
    # Filter genes based on expression data if provided
    if GeneExprBinaryRep_file:
        df_GeneExprBinaryRep = pd.read_csv(GeneExprBinaryRep_file, sep='\t')
        expressed_genes = set(df_GeneExprBinaryRep.loc[df_GeneExprBinaryRep['Expressed'] == 1, 'Gene'])
        selected_genes = expressed_genes & set(df_variant_level_CADD['Gene'])
    else:
        selected_genes = set(df_variant_level_CADD['Gene'])
    
    df_filtered = df_variant_level_CADD[df_variant_level_CADD['Gene'].isin(selected_genes)]

    # Calculate gene-level impact score
    avg_cadd_score = df_filtered.groupby('Gene')['CADD_score'].mean()
    mapped_scores = avg_cadd_score.index.map(lambda gene: dict_DepMap.get(gene, 0))
    dict_gene_level_impact_score = avg_cadd_score.mul(mapped_scores).to_dict()

    # Filter out genes with an impact score of 0
    dict_gene_level_impact_score = {gene: score for gene, score in dict_gene_level_impact_score.items() if score != 0}
    
    return dict_gene_level_impact_score


def main():
    args = parse_args()
    dict_gene_level_impact_score = calculate_gene_level_impact_score(input_vcf=args.Vcf, DepMap_file=args.DepMap, GeneExprBinaryRep_file=args.GeneExprBinaryRep, cutoff_CADD=args.cutoff_CADD)
    
    # Output results
    output_file = f"{args.SampleName}.gene_level_impact_scores.txt"
    pd.DataFrame(list(dict_gene_level_impact_score.items()), columns=['Gene', 'ImpactScore']).to_csv(output_file, sep='\t', index=False)
    print(f"Gene-level impact scores have been saved to {output_file}")

if __name__ == "__main__":
    main()









        



