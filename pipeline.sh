#!/bin/bash
# Calculate T.brucei transcripts from both directions
# Mika Olami

# PARAMETERS
genomecov=$1;                            # genomecov file
transcripts=$2;                          # name of transcripts    
min_coverage=$3;                         # minimun coverage
max_gap=$4;                              # maximum gap allowed
genome=$(basename "$transcripts" .bed);  # transcripts output file name

display_usage() { 
  echo -e "\n~~~ Calculate T.brucei transcripts from both directions ~~~";
	echo -e "Usage: ~/scripts/pipeline.sh [genomecov output] [transcripts name] [min coverage] [max gap] \n"; 
	} 
 
# if not enough arguments supplied, display usage 
	if [  $# -lt 4 ] 
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
echo "Min Coverage,Max Gap,Total Transcripts,Number of Non-overlapping transcripts (with known annotations),Number of overlapping transcripts (with known annotations),TBsRNA Found,Number of Regions after Filtering,TBsRNA Found after Filtering" >> "$output_file";
echo "Running script and calculating transcripts on $genome Cov: $min_coverage ,Gap: $max_gap ...";
python3 ../calcTranscriptsRelativeGaps.py "$genomecov" $genome"_front.bed" $min_coverage $max_gap;
echo "Finished creating front transcripts";
python3 ../calcTranscriptsRelativeGaps.py "$genomecov_reverse" $genome"_back.bed" $min_coverage $max_gap;
echo "Finished creating back transcripts";
awk '{print $1 "\t" $3 "\t" $2}' $genome"_back.bed" > $genome"_back_rev.bed"
multiIntersectBed -i $genome"_front.bed" $genome"_back_rev.bed" -names front back | grep front,back | cut -f1,2,3 > $genome"_transcripts.bed";
transcripts=$genome"_transcripts.bed";
total_transcripts=$(cat $transcripts | wc -l);
non_overlapping=$(intersectBed -a $transcripts -b ~/Tb_bed_files/*.bed -v | wc -l);
overlapping=$(intersectBed -a $transcripts -b ~/Tb_bed_files/*.bed -u | wc -l);
total_tbsrna=$(intersectBed -wa -wb -a $transcripts -b ~/Tb_sRNA_todate_updated.bed | cut -f 7 | sort | uniq | wc -l);

echo "Filtering overlapping regions with known annotations (cds, snos, tRNA, miscRNA, rRNA, UsnRNA)...";
genome_no_overlaps="$genome"_"no_overlaps.bed";
subtractBed -a $transcripts -b ~/Tb_bed_files/*.bed > "$genome_no_overlaps";
no_overlaps_regions=$(cat $genome_no_overlaps | wc -l);
tbsrna_after_filtering=$(intersectBed -wa -wb -a $genome_no_overlaps -b ~/Tb_sRNA_todate_updated.bed | cut -f7 | sort | uniq | wc -l);
echo "$min_coverage","$max_gap","$total_transcripts","$non_overlapping","$overlapping","$total_tbsrna","$no_overlaps_regions","$tbsrna_after_filtering" >> "$output_file";
rm $genome"_front.bed" $genome"_back.bed" $genome"_back_rev.bed" "$genomecov_reverse";
echo "Script Finised.";
