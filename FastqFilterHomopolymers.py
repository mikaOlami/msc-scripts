#!/usr/bin/env python
# Mika Olami 18/05/2022
import sys
import re

#   Filter reads that contains homopolymers of at least X nts from a fastq file
#   USAGE:
#       python3 FastqFilterHomopolymers.py [fastq file] [homopolymer minimum length (default: 15)]

#   how to use
if len(sys.argv) < 2:
    exit("\nFilter reads that contains homopolymers of at least X length from a fastq file\n\tUSAGE:\t\n\t"
         "python3 FastqFilterHomopolymers.py [fastq file] [homopolymer minimum length (default: 15)]\n")

MinLength = 14 if len(sys.argv) == 2 else int(sys.argv[2]) - 1  # homopolymer minimum length allowed

if MinLength <= 0:
    exit("\nError. Homopolymer minimum length has to be greater than 1. Exiting...\n")

# parameters
ReadName = None
Sequence = None
ThirdLine = None
ReadScore = None
ParsedThirdLine = False


# reset the current parameters
def resetParams():
    global ReadName, ReadScore, Sequence, ThirdLine, ParsedThirdLine
    ReadName = None
    Sequence = None
    ThirdLine = None
    ReadScore = None
    ParsedThirdLine = False


# print all reads passing the filter of homopolymers
def printFiltered():
    if not [x for x in (ReadName, Sequence, ThirdLine, ReadScore) if x is None]:
        print(ReadName, Sequence, ThirdLine, ReadScore, sep='\n')

# parse fastq file
with open(sys.argv[1], "r") as fastq:
    for line in fastq:
        line = line.strip()
        # first line - read name
        if line.startswith('@'):
            ReadName = line
        # second line - sequence
        elif '+' not in line and '@' not in line and not ParsedThirdLine:
            # check for homopolymers of at least X (default 15)
            if re.search(r"([ATCG])\1{%s,}" % MinLength, line) is None:
                Sequence = line
        # third line
        elif '+' in line:
            ThirdLine = line
            ParsedThirdLine = True
        # fourth line - score
        elif '+' not in line and ParsedThirdLine:
            ReadScore = line
            printFiltered()
            ParsedThirdLine = False
            resetParams()
