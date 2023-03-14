#!/bin/bash

# updated on 06/07/2022

# calculate unknown transcripts of L.donovani using the output of genomecov


# PARAMETERS

genomecov=$1;      # genomecov output file
transcripts=$2;    # transcripts output name
min_coverage=$3;   # minimun coverage for detection of transcripts
max_gap=$4;        # maximun gap in coverage for detection of transcripts

genome=$(basename "$transcripts" .bed);
bed_file="/home/ls/mikao/LD/Annot_Ld1S/LD_annots.copies.bed";

display_usage() { 
  echo -e "\n~~~ Calculate L.donovani transcripts from both directions ~~~";
	echo -e "Usage: ~/scripts/Pipeline_Ld.sh [genomecov output] [transcripts name] [min coverage] [max gap] \n"; 
	} 
 
# if not enough arguments supplied, display usage 
	if [ $# -lt 4 ] 
	then 
		display_usage;
		exit 1;
	fi 


if [ ! -d "./both_directions_script_Transcripts" ]; then mkdir -p ./both_directions_script_Transcripts; fi;
dir="./both_directions_script_Transcripts";
genomecov=$(realpath "$genomecov");
cd "$dir";
tac "$genomecov" > $genome"_rev_genomecov.out";
genomecov_reverse=$genome"_rev_genomecov.out";
output_file="$genome.csv";

echo "creating $output_file";
echo "Min Coverage,Max Gap,Total Transcripts,Number of Non-overlapping transcripts (with known annotations),Number of overlapping transcripts (with known annotations),Number of Regions after Filtering" >> "$output_file";
echo "Running script and calculating transcripts on $genome. Min coverage: $min_coverage ,Gap: $max_gap ...";
python3 ~/scripts/calcTranscriptsRelativeGaps.py "$genomecov" $genome"_front.bed" $min_coverage $max_gap;
echo "Finished creating front transcripts";
python3 ~/scripts/calcTranscriptsRelativeGaps.py "$genomecov_reverse" $genome"_back.bed" $min_coverage $max_gap;
echo "Finished creating back transcripts";
awk '{print $1 "\t" $3 "\t" $2}' $genome"_back.bed" > $genome"_back_rev.bed"
multiIntersectBed -i $genome"_front.bed" $genome"_back_rev.bed" -names front back | grep front,back | cut -f1,2,3 > $genome"_transcripts.bed";
transcripts=$genome"_transcripts.bed";
total_transcripts=$(cat $transcripts | wc -l);
non_overlapping=$(intersectBed -a $transcripts -b $bed_file -v | wc -l);
overlapping=$(intersectBed -a $transcripts -b $bed_file -u | wc -l);

echo "Filtering overlapping regions with known annotations (cds, snos, tRNA, vtRNA, rRNA, 5S, 7SL, SLRNA, UsnRNA)...";
genome_no_overlaps="$genome"_"no_overlaps.bed";
subtractBed -a $transcripts -b $bed_file > "$genome_no_overlaps";
no_overlaps_regions=$(cat $genome_no_overlaps | wc -l);
echo "$min_coverage","$max_gap","$total_transcripts","$non_overlapping","$overlapping","$no_overlaps_regions" >> "$output_file";
cat $genome"_transcripts.bed" | sort -k1,1n -k2,2n -o $genome"_transcripts.bed"{,};
cat "$genome_no_overlaps" | sort -k1,1n -k2,2n -o "$genome_no_overlaps"{,};
rm $genome"_front.bed" $genome"_back.bed" $genome"_back_rev.bed" "$genomecov_reverse";
echo "Script Finished.";
