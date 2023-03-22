#!/bin/bash

#  Run_Chimera_pipeline.sh
#  FIND CHIMERAS USING:
#  trim_galore
#  cutadapt (trim SL sites)
#  flash
#  fastp (trim polyA)
#  prinseq-lite (trim polyA/T)
#  bwa
#  splash (find chimeras)

#  USAGE: nohup ~/scripts/Run_Chimera_pipeline.sh [R1 fastq file] >> log &

# updated 13/03/2022 - aligned unmerged R1&R2 reads with bwa aln and pairing them with bwa sampe
#                      and changed cutadapt to paired-end mode
# updated 17/05/2022 - added polyA/T trimming to reads 
#
# 12/12/2022 - filter secondary alignments
#
# 22/12/2022 - merge reads while allowing outies --><-- with flash option -O
#
# NEXT SCRIPT -> ~/scripts/ChimerasAnalysisPipeline.sh
#
#################################################################################################

display_usage() { 
  echo -e "\n~~~ Run chimera pipeline to find chimeras ~~~" ;
	echo -e "\nUsage: nohup ~/scripts/Run_Chimera_pipeline.sh [R1 fastq file] >> log & \n" ;
	} 
 
# if not enough arguments supplied, display usage 
	if [  $# -lt 1 ] 
	then 
		display_usage;
		exit 1;
	fi 

#  PARAMETERS

Lib_R1=$1;                             # R1 fastq file
base=$(basename $Lib_R1 R1.fastq.gz);  # Library basename
Lib_R2=$(echo $Lib_R1 | sed 's/R1/R2/g');               # R2 fastq file
flash_output=$(basename $Lib_R1 | sed -E 's/R1.*gz$/flash.out/g');  # flash output file names
val1=$(basename $Lib_R1 | sed -E 's/\.fastq.gz$/_val_1.fq.gz/g');   # trim galore R1 file name
val2=$(basename $Lib_R1 | sed -E 's/R1.*gz$/R2_val_2.fq.gz/g');     # trim galore R2 file name
SL_db="~/LD_Chimera_Libs/Sequencing_19Oct2021/SL_DB.out";         # SL DB file
SL_rev_db="~/LD_Chimera_Libs/Sequencing_19Oct2021/SL_rev_DB.out"  # reverse complement SL DB file
declare -i align_length=5;  # cutadapt parameter: Require MINLENGTH overlap between read and adapter for an adapter to be found

#######################

# trim adaptors with trim galore

echo "Trim galore ...";
echo "";

trim_galore --retain_unpaired --paired $Lib_R1 $Lib_R2 --stringency 5 --output_dir . ;

# trim with flash

echo "flash ...";
echo "";

flash -M 60 -O -o $flash_output $val1 $val2;

# trim SL sites with cutadapt

echo "cutadapt - trim SL sites (R1, R2 not combined)...";
echo "";
cutadapt -g file:$SL_db -a file:$SL_rev_db -O $align_length -m 15 -o "$base"trimmed_notCombined_1.fastq -p "$base"trimmed_notCombined_2.fastq "$flash_output.notCombined_1.fastq" "$flash_output.notCombined_2.fastq" >>"$base"cutadapt_log.txt 2>> log.txt;

echo "cutadapt - trim SL sites (R1+R2 combined)...";
echo "";
cutadapt -g file:$SL_db -a file:$SL_rev_db -O $align_length -m 15 -o "$base"trimmed_merged.fastq "$flash_output.extendedFrags.fastq" >>"$base"cutadapt_log.txt 2>> log.txt;

# trim polyA with fastp

echo -e "fastp - trim polyA ...\n";

fastp -x -i "$base"trimmed_merged.fastq -o "$base"trimmed_merged_fastp.fastq;

fastp -x -i "$base"trimmed_notCombined_1.fastq -I "$base"trimmed_notCombined_2.fastq -o "$base"trimmed_notCombined_1_fastp.fastq -O "$base"trimmed_notCombined_2_fastp.fastq;

# trim polyA and polyT tails with prinseq

echo -e "prinseq - trim polyA and polyT tails ...\n";

perl prinseq-lite.pl -fastq "$base"trimmed_merged_fastp.fastq -trim_tail_left 10 -trim_tail_right 10 -out_good "$base"trimmed_merged_prinseq -out_bad null -lc_method entropy -lc_threshold 60 -min_len 40 >> "$base"prinseq.log &

perl prinseq-lite.pl -trim_tail_left 10 -trim_tail_right 10 -out_good "$base"trimmed_notCombined_prinseq -out_bad null -lc_method entropy -lc_threshold 60 -min_len 40 -fastq "$base"trimmed_notCombined_1_fastp.fastq -fastq2 "$base"trimmed_notCombined_2_fastp.fastq >> "$base"prinseq.log &

wait;

# filter homopolymers with script

echo -e "Filtering homopolymers from reads ...\n";

python3 ~/scripts/FastqFilterHomopolymers.py "$base"trimmed_notCombined_prinseq_1.fastq > "$base"trimmed_notCombined_1_homopolymers.fastq~;

python3 ~/scripts/FastqFilterHomopolymers.py "$base"trimmed_notCombined_prinseq_2.fastq > "$base"trimmed_notCombined_2_homopolymers.fastq~;

python3 ~/scripts/FastqFilterHomopolymers.py "$base"trimmed_merged_prinseq.fastq > "$base"trimmed_merged_homopolymers.fastq;

# repair the paired-end reads with bbmap repair.sh

echo -e "repairing paired-end reads ...\n";

bbmap/repair.sh in1="$base"trimmed_notCombined_1_homopolymers.fastq~ in2="$base"trimmed_notCombined_2_homopolymers.fastq~ out1="$base"trimmed_notCombined_1_homopolymers.fastq out2="$base"trimmed_notCombined_2_homopolymers.fastq 2> "$base"repairPE.log;

# align reads to genome with bwa

echo "bwa - R1 not merged ...";
echo "";

nohup bwa aln -t 4 -f "$base"aln_sa1.sai ~/LD/LD_bwa_index "$base"trimmed_notCombined_1_homopolymers.fastq &

echo "bwa - R2 not merged ...";
echo "";

nohup bwa aln -t 4 -f "$base"aln_sa2.sai ~/LD/LD_bwa_index "$base"trimmed_notCombined_2_homopolymers.fastq &

echo "bwa - R1+R2 merged ...";
echo "";

nohup bwa-mem2 mem -t 4 ~/LD/LD_bwa_mem_index "$base"trimmed_merged_homopolymers.fastq > "$base"merged_bwa_mem.sam 2> "$base"merged_bwa.log &

wait;

echo "bwa sampe - pairing unmerged R1/R2 reads ...";
echo "";

bwa sampe -f "$base"bwa_mem_vs_LD_Genome.sam ~/LD/LD_bwa_index "$base"aln_sa1.sai "$base"aln_sa2.sai "$base"trimmed_notCombined_1_homopolymers.fastq "$base"trimmed_notCombined_2_homopolymers.fastq 2> "$base"sampe_trimmed_polyA.log;

# convert to bam with samtools ##################

echo "Converting to bam ...";
echo "";

samtools view -Sb -F 4 "$base"merged_bwa_mem.sam -o "$base"merged_bwa_mem.bam;

samtools view -Sb -F 4 "$base"bwa_mem_vs_LD_Genome.sam -o "$base"bwa_mem_vs_LD_Genome.bam;

# filter secondary alignments with samtools ###############

samtools view -Sb -F 256 "$base"merged_bwa_mem.bam | samtools sort -o "$base"merged_bwa_mem_filtered_sorted.bam;

samtools view -Sb -F 256 "$base"bwa_mem_vs_LD_Genome.bam | samtools sort -o "$base"bwa_mem_vs_LD_Genome_filtered_sorted.bam;

# index bam files with samtools 

echo "Indexing bam files ...";
echo "";

samtools index "$base"merged_bwa_mem_filtered_sorted.bam;

samtools index "$base"bwa_mem_vs_LD_Genome_filtered_sorted.bam;

# find chimeras on merged reads with splash

echo "splash - find chimeras on merged file ...";
echo "";

python splash-master/src/find_chimeras.py -i "$base"merged_bwa_mem_filtered_sorted.bam --min-chim-dist 40 > "$base"merged_chimeras.out;

# sort bam files by read names with samtools

echo "Sorting united bam file by read name...";
echo "";

samtools sort -n "$base"bwa_mem_vs_LD_Genome_filtered_sorted.bam -o "$base"notCombinedUnited_filtered_sorted.bam;

# find chimeras on unmerged reads with python script

echo "ExtractChimeras.py - find non-overlapping R1 & R2 chimeras ...";
echo "";

python3 ~/scripts/ExtractChimeras.py "$base"notCombinedUnited_filtered_sorted.bam > "$base"notCombinedUnited_chimeras.out;

# delete temporary files

echo "deleting sam files ...";
echo "";

rm "$base"merged_bwa_mem.sam "$base"bwa_mem_vs_LD_Genome.sam;

echo -e "deleting unsorted bam files ...\n";

rm "$base"merged_bwa_mem.bam "$base"bwa_mem_vs_LD_Genome.bam;

echo "deleting sai files ...";
echo "";

rm "$base"aln_sa1.sai "$base"aln_sa2.sai;

echo -e "deleting singleton reads ...\n";

rm "$base"trimmed_notCombined_prinseq_1_singletons.fastq "$base"trimmed_notCombined_prinseq_2_singletons.fastq;

echo -e "deleting temporary fastq files ...\n";

rm "$base"trimmed_notCombined_1_homopolymers.fastq~ "$base"trimmed_notCombined_2_homopolymers.fastq~

echo -e "Done! next script to run: ~/scripts/ChimerasAnalysisPipeline.sh\n";
