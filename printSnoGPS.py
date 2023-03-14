#!/usr/bin/env python
# Mika Olami 02/05/2021
import sys

#   Print best results from snoGPS
#   USAGE:
#       python3 printSnoGPS.py [snoGPS intermediate file]

#   how to use
if len(sys.argv) < 2:
    exit("\nPrint best results from snoGPS\nUSAGE:\t\n\t"
         "python3 printSnoGPS.py [snoGPS intermediate file]\n")

#  parameters
gene = None   # result's chromosome
coord = None  # result's coordinates
score = None  # result's snoGPS score
bulge = None  # result's number of bulged nucleotides in the internal stem
target = None # result's expected target site
box = None    # result's ACA/AGA box
rRNA_match = None      # result's length of matching bps to target rRNA
rRNA_mismatch = None   # result's length of mismatching bps to target rRNA
ext_stem_match = None  # result's length of matching bps to expected external stem
int_stem_match = None  # result's length of matching bps to expected internal stem
last_line_in_current_result = False  # flag indicating last line of current result


# reset the current result
def resetParams():
    global gene, coord, score, bulge, target, box, rRNA_match, rRNA_mismatch, ext_stem_match, int_stem_match, strand, \
        sequence, last_line_in_current_result
    gene = None
    coord = None
    score = None
    bulge = None
    target = None
    box = None
    rRNA_match = None
    rRNA_mismatch = None
    ext_stem_match = None
    int_stem_match = None
    strand = None
    sequence = None
    last_line_in_current_result = False


# print all results passing the filter - an AGA box and external stem length of at least 5 bps
def printResults():
    if not [x for x in (gene, coord, target, score, bulge, box, ext_stem_match, int_stem_match, rRNA_match,
                        rRNA_mismatch, strand, sequence) if x is None]:
        if box.startswith('AGA') and ext_stem_match >= 5:
            print(gene + ":" + str(coord[0]) + "-" + str(coord[1]), target, score, bulge, box,
                  ext_stem_match, int_stem_match, rRNA_match, rRNA_mismatch, strand, sequence, sep='\t')


# main - parse all results of snoGPS, print results if they pass the filter
with open(sys.argv[1], "r") as snoGPS:
    for res in snoGPS:
        res = res.strip()
        # first line of results - contains result information
        if res.startswith('>'):
            resetParams()
            gene = res.split()[0]
            gene = gene.strip('>').split('.')
            gene = gene[0] if len(gene) == 2 else '.'.join(gene[0::1])
            score = res.split()[1]
            coords = res.split()[2]
            p_start = int(coords.strip('()').split('-')[0])
            p_end = int(coords.strip('()').split('-')[1])
            g_start = int(gene.split(':')[1].split('-')[0])
            gene = gene.split(':')[0]
            coord = [g_start + p_start, g_start + p_end]
            target = res.split()[4]
            rRNA_match = int(res.split()[7].split('/')[0])
            rRNA_mismatch = int(res.split()[7].split('/')[1])
            int_stem_match = int(res.split()[7].split('/')[2])
            ext_stem_match = int(res.split()[7].split('/')[3])
            strand = res.split()[10].strip("()")

        # ACA/AGA box line
        elif 'HACA:' in res:
            box = res.split()[2]

        # sequence line
        elif '#' not in res and 'X' not in res and len(res) > 5:
            sequence = res

        # predicted H/ACA snoRNA matches to rRNA
        elif res.startswith('X'):
            left_side_rRNA_idx = res.rindex("L")
            first_internal_stem_idx = res.index("I")
            right_side_rRNA_idx = res.index("R")
            last_internal_stem_idx = res.rindex("I")
            bulge = (first_internal_stem_idx - left_side_rRNA_idx - 1) + (
                    right_side_rRNA_idx - last_internal_stem_idx - 1)
            last_line_in_current_result = True

        # last line of current result
        if last_line_in_current_result:
            printResults()
            last_line_in_current_result = False
