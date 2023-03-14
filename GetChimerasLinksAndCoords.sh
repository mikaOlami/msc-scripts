#!/bin/bash

# get chimeras links file and convert coords to annotations
#  USAGE: ~/scripts/GetChimerasLinksAndCoords.sh [chimeras file]

display_usage() { 
  echo -e "\n~~~ Get chimeras links file and convert coordinates to annotations ~~~";
	echo -e "\nUsage: ~/scripts/GetChimerasLinksAndCoords.sh [chimeras file] [optional: annotations bedfile]\n";
	} 
 
# if not enough arguments supplied, display usage 
	if [  $# -lt 1 ] 
	then 
		display_usage;
		exit 1
	fi 

chimeras_file=$1;  # chimeras file
if [ -z "$2" ]; then annots_bed="~/LD/Annot_Ld1S/LD_annots.copies.bed"; else annots_bed=$(realpath $2); fi;  # optional annotations bedfile


links_file=$(basename ${chimeras_file%.*})_links.out;  # links output file
relative_coords_output=$(basename ${chimeras_file%.*})_relative_coords.out;  # relative coords output file

~/scripts/ChimerasToLinks.sh $chimeras_file;
python3 ~/scripts/ConvertCoords.py $annots_bed $links_file > $relative_coords_output;

