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
    FHAND.close()
except:
    print 'File not found.'

# ---------------------------------------------------------------------------
# Read in full IMP network and determine set difference between IMP and KEGG
# and set new path_genes to exclude those genes.
# ---------------------------------------------------------------------------

try:
    # read in all IMP genes
    EDGE_LST = pd.read_table('global_average.filtered.txt', header=None)

    # set of unique IMP genes
    IMP_GENES = set(EDGE_LST[0]).union(set(EDGE_LST[1]))

    # convert IMP gene from int to string
    IMP_GENES = set([str(gene) for gene in IMP_GENES])

    # unlists genes
    ALL_GENES = set([genes for path in PATH_GENES for genes in path])

    # finds genes in ontology not in IMP
    DIFF = list(ALL_GENES.difference(IMP_GENES))

    # creates gene ontology wo/o genes in ontology  but not in IMP
    PATH_GENES2 = []
    for path in PATH_GENES:
        path = [gene for gene in path if gene not in DIFF]
        PATH_GENES2.append(path)

    PATH_GENES = PATH_GENES2

    # This version removes any genes in IMP but not in ontology
    # unlists genes
    ALL_GENES = set([genes for path in PATH_GENES for genes in path])

except:
    print 'Invalid edege-list.'
