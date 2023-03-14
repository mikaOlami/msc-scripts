#!/bin/bash
# Mika Olami 17/03/2022
# Updated 05/04/2022 - Added significant pairs analysis
# Updated 23/05/2022 - Changed line 51 to detect unannotated genes starting with "Chr:"
# Updated 24/05/2022 - Removed the expected and observed vs expected lists output at line 81,84

#   Run the analysis pipeline for a chimeras file
#   USAGE:
#        ~/scripts/ChimerasAnalysisPipeline.sh [chimeras file] [optional: annotations bedfile]
#
# NEXT SCRIPT -> ~/scripts/ChimerasFilterIntras.py
#
###############################################################################################

display_usage() { 
  echo -e "\n~~~ Run the analysis pipeline for a chimeras file ~~~";
	echo -e "~/scripts/ChimerasAnalysisPipeline.sh [chimeras file] [optional: annotations bedfile]\n"; 
	} 
 
# if not enough arguments supplied, display usage 
	if [  $# -lt 1 ] 
	then 
		display_usage;
		exit 1;
	fi 

chimeras_file=$1;  # chimeras file
if [ -z "$2" ]; then annots_bed="/home/ls/mikao/LD/Annot_Ld1S/LD_annots.copies.bed"; else annots_bed=$(realpath $2); fi;  # optional annotations bedfile

# output file names
links_out=$(basename $chimeras_file .out)_links_genes.out;
relative_coords=$(basename $chimeras_file .out)_relative_coords.out;
matrix=$(basename $chimeras_file .out)_matrix.csv;
gene_count=$(basename $chimeras_file .out)_gene_count_no_intra.out;

# Number of columns in chimeras file
ColNum=$(awk '{print NF; exit}' $chimeras_file);

# Convert chimeras to links with annotations
if [ $ColNum -gt 6 ]
  then
    ~/scripts/GetChimerasLinksAndCoords.sh $chimeras_file $annots_bed;
  else
    ~/scripts/GetChimerasCoords.sh $chimeras_file $annots_bed;
fi

# Save chimeras links - only geneA and geneB annotations (unannotated genes are named "unannotated")
awk '{if ($1~/^Chr:.*/ && $4~/^Chr:.*/) print "unannotated\tunannotated"; else if ($1~/^Chr:.*/) print "unannotated\t" $4; else if ($4~/^Chr:.*/) print $1 "\tunannotated"; else print $1 "\t" $4}' $relative_coords > $links_out;

# Save the chimeras matrix as a csv file
python3 ~/scripts/ChimerasMatrix.py $links_out $matrix;

# Save the non-self chimeras count for every gene in the library
awk -F , '{ for(i=1; i<=NF;i++) if (NR!=i) j+=$i; print $1 "\t" j; j=0 }' $matrix | tail -n +2 > $gene_count;

# Total non-self chimeras
total=$(awk '{s+=$2}END{print s}' $gene_count);


chimeras_count=$(cat $relative_coords | wc -l);
intra_genes_count=$(awk '$1==$2 {print}' $links_out | wc -l);
intra_percent=$(echo "scale=2;($intra_genes_count/$chimeras_count)*100" | bc);
inter_genes_count=$(awk '$1!=$2 {print}' $links_out | wc -l);
inter_percent=$(echo "scale=2;($inter_genes_count/$chimeras_count)*100" | bc);

echo -e "########### REPORT FOR: $chimeras_file ###########";
echo -e "Chimeras Count:\t\t$chimeras_count";
echo -e "Intra genes Count:\t$intra_genes_count ($intra_percent %)";
echo -e "Inter genes Count:\t$inter_genes_count ($inter_percent %)";