# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 21:58:59 2016

@author: LiaHarrington

Ontology_prep.py

Usage:

Used to prepare data for use in enrichment_testing and community_detection

Description:

Prepares ontology data for enrichment and community detection analysis.
First,it creates the list of list of genes in ontology and then it removes
any genes in KEGG but not in IMP. Read in and parse data into dictionary of
pathways to form gene ontology. Saves path names in path names and all genes
in each pathway as a list of list in path_genes. Assumes tab delination and
that pathway name is the first string and the gene symbols start as 3rd string.
PATH_NAMES and PATH_GENES are created global variables.

"""
import pandas as pd

#Initiate GLOBAL VARIABLES

PATH_NAMES = []
PATH_GENES = []

# Note ALL_GENES is also a global variable
# ---------------------------------------------------------------------------

# Choose ontology
REACTOME = 'c2.cp.reactome.v5.1.entrez.gmt.txt'
KEGG = 'c2.cp.kegg.v5.1.entrez.gmt.txt'
PID = 'PID.Entrez.DB.txt'

FILENAME = KEGG

try:
    FHAND = open(FILENAME)

    for line in FHAND:
        line = line.rstrip().split('\t')
        PATH_NAMES.append(line[0])
        PATH_GENES.append(line[2:len(line)])
except:
    print 'File not found.'

FHAND.close()

# ---------------------------------------------------------------------------
# Read in full IMP network and determine set difference between IMP and KEGG
# and set new path_genes to exclude those genes. There are 36 genes in KEGG
# not in the IMP network. Casey said this is ok and to be expected since IMP
# was developed in 2010 and the version may be out of date with newest KEGG.
# ---------------------------------------------------------------------------

try:
    EDGE_LST = pd.read_table('global_average.filtered.txt', header=None)  # read in all IMP genes
    IMP_GENES = set(EDGE_LST[0]).union(set(EDGE_LST[1]))  # set of unique IMP genes

    IMP_GENES = set([str(gene) for gene in IMP_GENES])

    ALL_GENES = set([genes for path in PATH_GENES for genes in path])  # unlists genes

    DIFF = list(ALL_GENES.difference(IMP_GENES))  # finds genes in KEGG not in IMP
    DIFF = [gene for gene in DIFF]  # converts back to string form

    PATH_GENES2 = []  # creates gene ontology wo/o genes in KEGG  but not in IMP
    for path in PATH_GENES:
        path = [gene for gene in path if gene not in DIFF]
        PATH_GENES2.append(path)

    PATH_GENES = PATH_GENES2

    # This version removes any genes in IMP but not in ontology
    ALL_GENES = set([genes for path in PATH_GENES for genes in path])  # unlists genes

except:
    print 'Invalid edege-list.'
