######################################
# Calculate Transcipts Using Genomecov
# Author: Mika Olami
# Date: 25/04/2021
######################################

#!/usr/bin/env python
import sys
import math


current_start=0  # transcript start
current_end=0    # transcipt end
current_chr=''   # current chromosome
cov_peak=0       # current maximum coverage (peak)
bed_dict = {}    # dictionary that maps chromosomes to a list of [start, end] transcripts
gap_count = 0    # gap counter

# how to use
if len(sys.argv) < 4:
    exit("Usage: python3 get_genomecov_transcripts.py [genomecov results] [output file] "
         "[min coverage # per tanscript] [max gap length]")

min_coverage = int(sys.argv[3])  # minimal coverage allowed for transcripts detection
max_gap = int(sys.argv[4])       # max gap length allowed for a transcript

# read the genomecov results
with open(sys.argv[1],"r") as genomecov_coverge:
    # read each line
    for base in genomecov_coverge:
        base=base.strip()
        base=base.split("\t")
        if float(base[2]) >= min_coverage:  # check the coverage
            cov_peak=max(float(base[2]), cov_peak)
            min_coverage=max(math.ceil(cov_peak*0.05), min_coverage)
            if len(current_chr) == 0 or current_chr != base[0]:  # check the chromosome
                if len(current_chr)>0 and [current_start,current_end] not in bed_dict[current_chr]:
                    bed_dict[current_chr].append([current_start, current_end])
                current_chr = base[0]  # current chromosome
                cov_peak=float(base[2])
                current_start = int(base[1])  # transcript start
            current_end = int(base[1])        # transcript end
            gap_count=0  # reset gap counter
            continue
        elif (len(current_chr) > 0 and gap_count < max_gap):  # allowed gaps
            gap_count += 1
            continue
        else:  # passed the max gap length
            gap_count += 1
            if gap_count > max_gap:
                if len(current_chr) > 0:
                    if current_chr not in bed_dict.keys():
                        bed_dict[current_chr] = []  # add a list to the dictionary
                    if current_start != current_end:
                        bed_dict[current_chr].append([current_start, current_end]) # add the transcript [start, end] to the dictionary
                    current_chr = ''  # reset current chromosome
                    gap_count = 0  # reset gap counter
                    cov_peak = 0
                    min_coverage = int(sys.argv[3])

# input file closed

# write the output to a BED file (chr   start   end)
with open(sys.argv[2],"w") as output_file:
    for key in sorted(bed_dict.keys()):
        for line in sorted(bed_dict[key]):
            output_file.write("%s\t%d\t%d\n" %(key,line[0],line[1]))

# output file closed