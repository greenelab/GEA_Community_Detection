# -*- coding: utf-8 -*-
"""
Created on Fri May 27 20:56:51 2016

@author: LiaHarrington
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 21:58:59 2016

@author: LiaHarrington
"""
import pandas as pd
from sets import Set
import re

# -------------------------- Code Description -------------------------------
# Prepares ontology data for enrichment and community detection analysis. First,
# it creates the list of list of genes in ontology and then it removes any genes in KEGG
# but not in IMP. 
# ---------------------------------------------------------------------------
# Read in and parse data into dictionary of pathways to form gene ontology 
# Saves path names in path names and all genes in each pathway as a list of list
# in path_genes. Assumes tab delination and that pathway name is the first string
# and the gene symbols start as 3rd string. PATH_NAMES and PATH_GENES are global variables.
# ---------------------------------------------------------------------------
#Initiate GLOBAL VARIABLES 

PATH_NAMES = []
PATH_GENES = []

# Note ALL_GENES is also a global variable 
# ---------------------------------------------------------------------------

# Choose ontology 
reactome = 'c2.cp.reactome.v5.1.entrez.gmt.txt'
kegg = 'c2.cp.kegg.v5.1.entrez.gmt.txt'
pid = 'PID.Entrez.DB.txt'

filename = kegg

try: 
    fhand = open(filename)

    for line in fhand:
        line = line.rstrip().split('\t')
        PATH_NAMES.append(line[0])
        PATH_GENES.append(line[2:len(line)])
except:
    print 'File not found.'

fhand.close()

# ---------------------------------------------------------------------------
# Read in full IMP network and determine set difference between IMP and KEGG
# and set new path_genes to exclude those genes. There are 36 genes in KEGG
# not in the IMP network. Casey said this is ok and to be expected since IMP
# was developed in 2010 and the version may be out of date with newest KEGG.
# ---------------------------------------------------------------------------
    
try: 
    edge_lst = pd.read_table('global_average.filtered.txt', header=None)     # read in all IMP genes
    imp_genes = set(edge_lst[0]).union(set(edge_lst[1]))       # set of unique IMP genes
    
    imp_genes = set([str(gene) for gene in imp_genes])
    
    ALL_GENES = set([genes for path in PATH_GENES for genes in path]) # unlists genes 
    
    diff = list(ALL_GENES.difference(imp_genes))    # finds genes in KEGG not in IMP
    diff = [gene for gene in diff]             # converts back to string form 
        
    path_genes2 = []        # creates gene ontology wo/o genes in KEGG  but not in IMP
    for path in PATH_GENES:
        path = [gene for gene in path if gene not in diff]
        path_genes2.append(path)
        
        
    PATH_GENES = path_genes2
    
    # This version removes any genes in IMP but not in ontology
    ALL_GENES = set([genes for path in PATH_GENES for genes in path]) # unlists genes 

except:
    print 'Invalid edege-list.'




