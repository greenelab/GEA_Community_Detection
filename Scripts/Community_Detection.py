# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 10:10:37 2016

@author: LiaHarrington
"""
execfile('Enrichment_prep.py')

# -------------------------- Code Description -------------------------------
# Creates network from IMP and subgraph using m selected pathways
# ---------------------------------------------------------------------------
# Import packages and set seed 
# ---------------------------------------------------------------------------
import igraph 
from igraph import *
from sets import Set
import random
import pandas 

random.seed(123)
# ---------------------------------------------------------------------------
# Reads in IMP ege-list and creates complete IMP network  
# ---------------------------------------------------------------------------

# read in IMP network 
imp = igraph.Graph.Read_Ncol('global_filtered_short.txt',directed=False)
    
# ---------------------------------------------------------------------------
# Creates subgraph with genes from the m pathways 
# ---------------------------------------------------------------------------
    
def weights(): 
    '''Read in network weights from text file and saves it locally'''    
    m_weights = open('weights_full.txt')
    
    weights = []
    for w in m_weights:
        weight = float(w.rstrip().split()[0])
        weights.append(weight)
    return weights
    
weights = weights()

def community_detection(sub_genes, method, weights): 
    '''Performs community detection using the genes from the m
    pathways to create the relevant subgraphs with only those gene nodes'''
     
    sub_s = imp.subgraph(sub_genes)
    
    if weights != 'NULL': 
        if method == 'fastgreedy': 
            cl = sub_s.community_fastgreedy(weights).as_clustering()
        if method == 'walktrap': 
            cl = sub_s.community_walktrap(weights, steps = 20).as_clustering()
        if method == 'infomap': 
            cl = sub_s.community_infomap(weights)
        if method == 'multilevel': 
            cl = sub_s.community_multilevel(weights)
    elif weights == 'NULL':
        if method == 'fastgreedy': 
            cl = sub_s.community_fastgreedy().as_clustering()
        if method == 'walktrap': 
            cl = sub_s.community_walktrap().as_clustering()
        if method == 'infomap': 
            cl = sub_s.community_infomap()
        if method == 'multilevel': 
            cl = sub_s.community_multilevel()
    
    return cl
    
def index_to_edge_name(cl): 
    '''Transforms the index labeled nodes in each community from cl to thier actual names'''
    cl_named_groups = []    
    for group in range(len(cl)): 
        gp = [cl.graph.vs['name'][cl[group][i]] for i in range(len(cl[group]))]
        cl_named_groups.append(gp)
    return cl_named_groups






















