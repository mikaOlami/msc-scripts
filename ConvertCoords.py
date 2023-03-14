#!/usr/bin/env python
import sys

#   Convert chimera links coordinates to annotations
#   USAGE:
#   python3 ConvertCoords.py [annotations bed file] [chimeras links file]

#   annotations bed file (filewcoords): a bed file contatining all genes annotations
#   chimeras links file (filetotranslate): a file containing the chimeras links genes (chrom1,start1,end1,chrom2,start2,end2)

# updated 23/05/2022 - Added "Chr:" at the start of genes without annotation at lines 80,91,100,111

#   how to use
if len(sys.argv) < 3:
    exit("Usage:\n\tpython3 ConvertCoords.py [annotations bed file] [chimeras links file]")

RNA_array = {}  # Chrom[position] -> Annotation
conversion = {}  # Chrom[position] -> annotation's relative position
annotations = {}  # Chrom[position] -> annotation's [start, end] positions


# update dictionaries
def updateDicts(s, e, chrom, nm):
    if s < e:  # strand plus
        for p in range(s, e + 1):
            if p not in RNA_array[chrom]:
                RNA_array[chrom][p] = nm
                conversion[chrom][p] = (p - s) + 1
            if annotations[chrom].get(p) is None:
                annotations[chrom][p] = [s, e]
            else:
                existing_dist = p - annotations[chrom][p][1]
                current_dist = s - p
                if current_dist > existing_dist:
                    RNA_array[chrom][p] = nm
                    conversion[chrom][p] = (p - s) + 1
                    annotations[chrom][p] = [s, e]
    else:  # strand minus
        for p in range(e, s + 1):
            if p not in RNA_array[chrom]:
                RNA_array[chrom][p] = nm
                conversion[chrom][p] = (s - p) + 1
            if annotations[chrom].get(p) is None:
                annotations[chrom][p] = [s, e]


# create RNA_array and conversion dictionaries
with open(sys.argv[1], "r") as filewcoords:
    for line in filewcoords:
        l = line.strip()
        l = l.split("\t")
        RNA = l[0]  # chromosome
        start = int(l[1])
        endp = int(l[2])
        strand = l[5]  # gene strand
        name = l[3]    # gene annotation
        if RNA not in annotations:
            annotations[RNA] = {}
        if RNA not in RNA_array:
            RNA_array[RNA] = {}
        if RNA not in conversion:
            conversion[RNA] = {}
        updateDicts(start, endp, RNA, name)

# convert chimeras coordinates to annotations and print results
with open(sys.argv[2], "r") as filetotranslate:
    for line in filetotranslate:
        l = line.strip()
        l = l.split("\t")
        chrom1 = l[0]
        start1 = int(l[1])
        end1 = int(l[2])
        chrom2 = l[3]
        start2 = int(l[4])
        end2 = int(l[5])
        if chrom1 not in RNA_array:  # unannotated gene1
            print("Chr:" + chrom1, start1, end1, end='\t', sep='\t')
        elif chrom1 in RNA_array and RNA_array[chrom1].get(start1) is not None:
            if conversion[chrom1].get(end1) is not None:
                print(RNA_array[chrom1][start1], conversion[chrom1][start1], conversion[chrom1][end1], end='\t',
                      sep='\t')
            else:
                annot_end = annotations[chrom1][start1][1]
                downstream = end1 - annot_end
                print(RNA_array[chrom1][start1], conversion[chrom1][start1], end='\t', sep='\t')
                print(conversion[chrom1][annot_end], '+', downstream, end='\t', sep='')
        else:
            if RNA_array[chrom1].get(start1) is None and RNA_array[chrom1].get(end1) is None:
                print("Chr:" + chrom1, start1, end1, end='\t', sep='\t')
            elif RNA_array[chrom1].get(start1) is None and RNA_array[chrom1].get(end1) is not None:
                annot_start = annotations[chrom1][end1][0]
                upstream = annot_start - start1
                print(RNA_array[chrom1][end1], end='\t', sep='\t')
                print("-", upstream, end='\t', sep='')
                print(conversion[chrom1][end1], end='\t', sep='\t')

        if chrom2 not in RNA_array: # unannotated gene2
            print("Chr:" + chrom2, start2, end2, sep='\t')
        elif chrom2 in RNA_array and RNA_array[chrom2].get(start2) is not None:
            if conversion[chrom2].get(end2) is not None:
                print(RNA_array[chrom2][start2], conversion[chrom2][start2], conversion[chrom2][end2], sep='\t')
            else:
                annot_end = annotations[chrom2][start2][1]
                upstream = end2 - annot_end
                print(RNA_array[chrom2][start2], conversion[chrom2][start2], end='\t', sep='\t')
                print(conversion[chrom2][annot_end], '+', upstream, sep='')
        else:
            if RNA_array[chrom2].get(start2) is None and RNA_array[chrom2].get(end2) is None:
                print("Chr:" + chrom2, start2, end2, sep='\t')
            elif RNA_array[chrom2].get(start2) is None and RNA_array[chrom2].get(end2) is not None:
                annot_start = annotations[chrom2][end2][0]
                upstream = annot_start - start2
                print(RNA_array[chrom2][end2], end='\t', sep='\t')
                print("-", upstream, end='\t', sep='')
                print(conversion[chrom2][end2], sep='\t')
