#!/usr/bin/env python
# Mika Olami 27/03/2022
# Updated 25/05/2022 - Added unannotated&gene chimeras
import sys

#   Tag genes found in the chimeras links file from the annotations type file
#   USAGE:
#       python3 tagChimeras.py [chimeras links file] [annotations to type file]

#   how to use
if len(sys.argv) < 3:
    exit("\nTag genes found in the chimeras links file from the annotations type file\nUsage:\n\t python3 "
         "tagChimeras.py [chimeras links file] [annotations to type file]\n")

GenesDict = {}        # chimeras genes distribution dictionary: GeneA[GeneB] => # of geneA-geneB chimeras
AnnotationsDict = {}  # Annotations dictionary: Annotation => gene type


# initialize the annotations dictionary
def initAnnotDict():
    with open(sys.argv[2], "r") as annots_file:
        for annot in annots_file:
            ann = annot.strip()
            ann = ann.split('\t')
            AnnotationsDict[ann[0]] = ann[1]


# calculate the chimera genes distribution and update in the dictionary
def calcChimeraGenes():
    # read chimeras links (annot1, annot2 in a chimera)
    with open(sys.argv[1], "r") as chimeras_links:
        for chimera in chimeras_links:
            chim = chimera.strip()
            chim = chim.split('\t')
            # get annotation gene type
            geneA = AnnotationsDict[chim[0]]
            geneB = AnnotationsDict[chim[1]]
            if geneA not in GenesDict:
                GenesDict[geneA] = {}
            if geneB not in GenesDict:
                GenesDict[geneB] = {}
            # initialize count to 1, and add 1 for every chimera
            if geneB not in GenesDict[geneA]:
                GenesDict[geneA][geneB] = 1
            else:
                GenesDict[geneA][geneB] += 1
            if geneA not in GenesDict[geneB] and geneA != geneB:
                GenesDict[geneB][geneA] = 1
            elif geneA != geneB:
                GenesDict[geneB][geneA] += 1


# remove duplicates from dictionary
def removeDuplicateKeys():
    for geneA in GenesDict:
        for geneB in GenesDict[geneA]:
            if geneA != geneB:
                GenesDict[geneB].pop(geneA, None)


# print the chimeras genes distribution
def printChimeraGenesDistribution():
    for geneA in GenesDict:
        for geneB in GenesDict[geneA]:
            print(geneA, "-", geneB, ":\t", GenesDict[geneA][geneB], sep="")


# main
initAnnotDict()
calcChimeraGenes()
removeDuplicateKeys()
printChimeraGenesDistribution()
