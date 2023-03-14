#!/usr/bin/env python
# Mika Olami 31/01/2022
import sys
import pandas as pd

#   Create a Chimeras matrix from chimeras links file
#   Chimeras Matrix: geneA, geneB, count of geneA+geneB chimeras
#   Output is saved as a csv file
#   USAGE:
#       nohup python3 ChimerasMatrix.py [chimeras links file] Optional: output file name (default: output.csv) &

#   how to use
if len(sys.argv) < 2:
    exit("\tUSAGE:\n\tpython3 ChimerasMatrix.py [chimeras links file (geneA geneB)] Optional: output file name (default: output.csv)\n")
    
outfile = sys.argv[2] if len(sys.argv) == 3 else 'output.csv'  # define output file name  
chimera_dict = {} # chimeras dictionary: chimera_dict[geneA] = {geneB: count, geneC: count, ...}

# parse the chimeras links file and update each geneA and geneB
# in dictionary with their appearance count
with open(sys.argv[1],"r") as chimeras_links:
    for link in chimeras_links:
        l=link.strip()
        l=l.split("\t")
        geneA=l[0]
        geneB=l[1]
        if geneA not in chimera_dict:  # initialize chimeras {geneA, geneB} with one appearance
            chimera_dict[geneA] = {}
            chimera_dict[geneA][geneB] = 1
        else:                            # add 1 for every chimeras {geneA, geneB} appearance
            if geneB not in chimera_dict[geneA]:
              chimera_dict[geneA][geneB] = 1
            else:
              chimera_dict[geneA][geneB] = int(chimera_dict[geneA][geneB])+1
        if geneB != geneA and geneB not in chimera_dict:  # for symmetry, initialize chimeras {geneB, geneA} with one appearance (geneA!=geneB)
            chimera_dict[geneB] = {}
            chimera_dict[geneB][geneA] = 1
        elif geneB != geneA:
            if geneA not in chimera_dict[geneB]:
              chimera_dict[geneB][geneA] = 1
            else:                                         # for symmetry, add 1 for every chimeras {geneB, geneA} appearance
              chimera_dict[geneB][geneA] = int(chimera_dict[geneB][geneA])+1

            
# convert dictonary to matrix (data frame)
matrix=pd.DataFrame.from_dict(chimera_dict, orient='index')
# sort by row and column headers
matrix=matrix.sort_index(axis=0)
matrix=matrix.sort_index(axis=1)
matrix.to_csv(outfile, index=True)
print("Output matrix saved at",outfile)