#!/usr/bin/env python
# Mika Olami 27/04/2022
# Updated 03/05/2022 - Added UTRs to clusters with mRNAs
import sys

# Filter out intra chimeras and print chimeras links
# Input: chimeras links file (geneA, geneB) and chimeras copies IDs (cluster of same gene)
# USAGE: python3 ChimerasFilterIntras.py [chimeras links] [genes clusters]

#   how to use
if len(sys.argv) < 3:
    exit(
        "\tUSAGE:\n\tpython3 ChimerasFilterIntras.py [chimeras links] [genes clusters]\n")

#  parameters
clustersDict = {}  # gene ID -> all gene IDs in same cluster (copy)


# create the genes clusters dictionary
def createClustersDict():
    with open(sys.argv[2], "r") as clusters:
        for cluster in clusters:
            cluster = cluster.replace("(-)", "")
            cluster = cluster.replace("(+)", "")
            c = cluster.strip()
            c = c.split()
            # for every gene cluster, set the other genes in the cluster in dictionary
            # initialize dictionary for every gene in the cluster
            for idx in range(len(c)):
                clustersDict[c[idx]] = c[0:idx] + c[idx + 1:]
                clustersDict[c[idx]].append(c[idx])
                # add 5/3pUTRs to mRNA genes in the cluster
                for gene in clustersDict[c[idx]]:
                    if gene.startswith("Ld1S") and "UTR" not in gene:
                        clustersDict[c[idx]].append(gene + "_5pUTR")
                        clustersDict[c[idx]].append(gene + "_3pUTR")
                    elif gene.startswith("Ld1S") and "5pUTR" not in gene and "3pTUR" in gene:
                        clustersDict[c[idx]].append(gene + "_3pUTR")
                    elif gene.startswith("Ld1S") and "3pUTR" not in gene and "5pTUR" in gene:
                        clustersDict[c[idx]].append(gene + "_5pUTR")


# update the clusters of UTR genes in the clusters dictionary
def updateClustersDictUTRs():
    for key in clustersDict.keys():
        UTRs = [gene for gene in clustersDict[key] if "UTR" in gene]
        for u in UTRs:
            if u in clustersDict:
                for c in clustersDict[key]:
                    if c not in clustersDict[u]:
                        clustersDict[u].append(c)


# filter the intra/copies chimeras pairs from the chimeras links
def filterIntraChimeras():
    with open(sys.argv[1], "r") as chimeras_links:
        for link in chimeras_links:
            l = link.strip()
            l = l.split()
            # print chimeras only if their not intra/copies
            if (l[0] not in clustersDict or l[1] not in clustersDict[l[0]]) and \
                    (l[1] not in clustersDict or l[0] not in clustersDict[l[1]]) and \
                    l[0] != l[1]:
                print(l[0], l[1], sep='\t')


createClustersDict()      # create clusters dictionary
updateClustersDictUTRs()  # update the clusters of UTR genes
filterIntraChimeras()     # print chimeras filtered for intra/copies
