#!/bin/bash

# run genomecov on the SORTED bam file

# PARAMETERS
input_genome_bam=$1;                              # input - sorted bam file
input_genome=$(basename $input_genome_bam .bam);  # bam file name without a suffix
genomecov_res=$input_genome"_genomecov.out";      # output file name

display_usage() { 
  echo -e "\n~~~ Run genomecov on a sorted bam file ~~~";
	echo -e "Usage: ~/scripts/get_genomecov.sh [sorted bam file] \n"; 
	} 
 
# if not enough arguments supplied, display usage 
	if [ $# -lt 1 ] 
	then 
		display_usage;
		exit 1;
	fi 

echo "Running genomecov...";
bedtools genomecov -ibam "$input_genome_bam" -d > "$genomecov_res";
echo "Saved genomecov results into $genomecov_res";