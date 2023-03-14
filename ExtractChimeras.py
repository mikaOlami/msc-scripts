#!/usr/bin/env python
import sys
import pysam

#   Extract chimeras from united file of non-overlapping R1 & R2
#   Input must be sorted by readname !
#   USAGE:
#       python3 ExtractChimeras.py [united R1&R2 bam file]

#   how to use
if len(sys.argv) < 2:
    exit("\nExtract chimeras from united file of non-overlapping R1 & R2\nUsage:\n\t python3 ExtractChimeras.py [united R1&R2 bam file]\n")
    
bampath = sys.argv[1]

# load BAM file
bamfile = pysam.AlignmentFile(bampath, 'rb')
bam_iter = bamfile.fetch(until_eof=True)

# iterate through reads
read1 = None
read2 = None
for read in bam_iter:
    read2 = read
    # check for read pairs
    if read1 is not None and read1.query_name == read2.query_name:
        chrom1 = bamfile.getrname(read1.tid)
        chrom2 = bamfile.getrname(read2.tid)
        start1 = read1.reference_start
        start2 = read2.reference_start
        end1 = read1.reference_end
        end2 = read2.reference_end
        dist = min(abs(end2-start1),abs(end1-start2))
        # check if the pair is in different chromosomes or more than 1000 nt apart - report chimeras
        if chrom1 != chrom2 or dist>1000:
            print(chrom1, start1, end1, chrom2, start2, end2, sep='\t')     
    read1 = read

bamfile.close()
