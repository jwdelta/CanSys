#!/usr/bin/env bash
set -euo pipefail

# Script Description:
# Filter an ANNOVAR-annotated VCF by gnomAD v4.1 population specific allele frequency.

# Usage:
# bash filter_vcf_by_popAF.sh <vcf_file> <pop> <cutoff>

if [[ $# -lt 3 ]]; then
    echo "Usage: bash $0 <vcf_file> <pop> <cutoff>"
    exit 1
fi

vcf_in="$1"
pop_input="$2"
cutoff="$3"

# Extract sample name for output naming
sample=$(basename "$vcf_in" .hg38_multianno.vcf)

# Map population code to ANNOVAR gnomAD v4.1 exome AF key
pop_uc="${pop_input^^}"
case "$pop_uc" in
    ALL ) AF_KEY="gnomad41_exome_AF" ;;
    AFR ) AF_KEY="gnomad41_exome_AF_afr" ;;
    AMR ) AF_KEY="gnomad41_exome_AF_amr" ;;
    ASJ ) AF_KEY="gnomad41_exome_AF_asj" ;;
    EAS ) AF_KEY="gnomad41_exome_AF_eas" ;;
    FIN ) AF_KEY="gnomad41_exome_AF_fin" ;;
    MID ) AF_KEY="gnomad41_exome_AF_mid" ;;
    NFE ) AF_KEY="gnomad41_exome_AF_nfe" ;;
    SAS ) AF_KEY="gnomad41_exome_AF_sas" ;;
    OTH ) AF_KEY="gnomad41_exome_AF_remaining" ;;
    * ) echo "Error: unsupported pop '$pop_input' (valid: ALL, AFR, AMR, ASJ, EAS, FIN, MID, NFE, SAS, OTH)"; exit 1 ;;
esac

vcf_out="${sample}.hg38_multianno.filtered.vcf"

# Print a small header to stderr for readability
#echo -e "CHR\tPOS\tREF\tALT\t${AF_KEY}\tSTATUS" >&2

awk -v AF_KEY="$AF_KEY" -v CUTOFF="$cutoff" '
BEGIN{ FS=OFS="\t"; }
/^##/ { print; next; }
/^#/  { print; next; }
{
    chr=$1; pos=$2; ref=$4; alt=$5;
    info=$8; af_max=0.0; target=AF_KEY "="; found=0;

    n=split(info, arr, ";");
    for(i=1;i<=n;i++){
	if(index(arr[i], target)==1){
            split(arr[i], kv, "=");
	    v=kv[2];
	    if(v=="" || v=="."){
		af_max=0.0;
	    } else {
	        m=split(v, vals, ",");
		for(j=1;j<=m;j++){
		    val=(vals[j]=="" || vals[j]==".") ? 0.0 : vals[j]+0;
		    if(val>af_max) af_max=val;
	        }
	    }
	    found=1;
	    break;
	}
    }
    if(!found) af_max=0.0;
    
    status = (af_max > CUTOFF+0) ? "DROP" : "KEEP";

    # show ALL ALTs as given to ease debugging
    #printf("%s\t%s\t%s\t%s\t%.6g\t%s\n", chr, pos, ref, alt, af_max, status) > "/dev/stderr";
    
    # write to filtered VCF only when kept
    if(status=="KEEP") print;
}
' "$vcf_in" > "$vcf_out"
