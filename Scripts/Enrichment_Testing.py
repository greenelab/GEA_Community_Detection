# -*- coding: utf-8 -*-
"""
Created on Fri May 27 20:54:32 2016

@author: LiaHarrington
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 12:51:01 2016

@author: LiaHarrington
"""

execfile('Enrichment_prep.py')
execfile('community_net2.py')
# -------------------------- Code Description -------------------------------
# Completes enrichment analysis of m randonly chosen pathways from the an
# ontology N times with and without community detection 
# ---------------------------------------------------------------------------
# Import packages and set seed 
# ---------------------------------------------------------------------------
import random 
from scipy.stats import hypergeom
import scipy.stats as stats
from sets import Set
import numpy as np
import matplotlib.pyplot as plt
import operator 

# ---------------------------------------------------------------------------
# Perform enrichment analysis of geneList with pathways in ontology and returns
# a  three entry result where the first entry is the p-values for all
# 186 pathways, the second column is True or False depending if 
# meet bonferroni adj p, and entry is number of significant pathways 
# ---------------------------------------------------------------------------
def enrichment(geneList,ontology, alpha, ALL_GENES):
    '''Returns p-values indicating enrichment level of genes in geneList in ontology
    where M = # of distinct genes in ontology, n = # of genes in ith
    pathway, N = # of genes in genelist, and k = # of genes in genelist that are 
    also in pathway. This returns a tuple of the pathway number and associated 
    T/F for significantly enriched in ranked order by most signif p-val.'''

    M = len(ALL_GENES)
    geneList = set(geneList)
    pvals = []
    
    for path in range(len(ontology)): 
        pathway = set(ontology[path])
        n = len(pathway)
        N = len(geneList)
        k = len(pathway.intersection(geneList))
        pvals.append(hypergeom.sf(k-1,M,n, N))
    
    p_tups = [(p, pvals[p]) for p in range(len(pvals))]
    p_tups.sort(key=operator.itemgetter(1))
        
    bonf_alpha = alpha/len(ontology) 
    signif = [pval < bonf_alpha for pval in pvals]
    
    sig_tup = [(p_tups[i][0],p_tups[i][1] < bonf_alpha) for i in range(len(p_tups))]

    return [sig_tup, signif, pvals, 'Number pathways significant: {0}'.format(sum(signif))]

# ---------------------------------------------------------------------------
# Randomly selects n% of m pathways and a% (using floor function) additional random 
# genes and then writes to text file. m is any whole number, n and a are any number
# between 0 and 1. 
# ---------------------------------------------------------------------------

def m_gene_list(m, n, a):
    selected_path_ids = random.sample(range(len(PATH_GENES)), m)
    geneList = []
    for selected_path in selected_path_ids:
        path = PATH_GENES[selected_path]
        path = random.sample(path, int(float(n)*len(path)))
        geneList += path
    geneList += random.sample(ALL_GENES, int(float(a)*len(geneList)))
    return [selected_path_ids, geneList]
      
def write_gene_list(geneList):
    text_file = open('m_genes.txt', 'w')

    for g in geneList:
        text_file.write(g) 
        text_file.write("\n")
    text_file.close()

# ---------------------------------------------------------------------------
# Performs N experiments of draws m pathways from ontology and prints the nuber
# of false positives and power accounting only the most significanlty
# enriched path. 
# ---------------------------------------------------------------------------

# Control arm 
def Control_most_signif(iterations, m, n, a, alpha): 
    ''' Simulation of N iterations of m chosen paths using n% of each path with a%
    additional random genes. Returns the average power and false positives. Power and 
    success here is defined for if the m most signif paths returned are the seeded paths. 
    Com_method can be 'fastgreedy', 'walktrap', 'infomap', or 'multilevel' ''' 
    
    random.seed(123)
    total_false_positive = 0
    successes = 0
    
    for iteration in range(iterations):
        gl = m_gene_list(m,n,a)
        selected_path_ids = gl[0]
        geneList = gl[1]
        
        # results returns [sig_tup, signif (T/F), pvals, # signif]
        results = enrichment(geneList, PATH_GENES,alpha ,ALL_GENES)
        
        # get the m most signif results 
        relevant_results = results[0][0:m]
        
        # extract path number 
        sig_paths = [relevant_results[i][0] for i in range(len(relevant_results))]
    
        
        # check if top m results equal to the m seeded paths 
        if sum([relevant_results[i][1] for i in range(len(relevant_results))]) == m and set(sig_paths)==set(selected_path_ids):
            successes+=1
        
        total_positive = sum(results[1])
        
        # subtract out corectly detected paths from all detected 
        total_false_positive += total_positive - len((set(selected_path_ids).intersection(set(sig_paths))))
    
        
    power = float(successes)/iterations 
    avg_false_pos = float(total_false_positive) / iterations 
    return [power, avg_false_pos]
    

# Experimental arm 
def Experimental_most_signif(com_method, iterations, m, n, a, weights, alpha): 
    ''' Simulation of N iterations of n% of m randomly chosen paths and a% additional
    random genes. Returns the average power and false positives. Weights can either 
    be "weights" using IMP weights or 'NULL' to use no weights. Com_method can be
    'fastgreedy', 'walktrap', 'infomap', or 'multilevel'. '''
    
    random.seed(123)
    total_false_positive = 0
    successes = 0 
    
    for iteration in range(iterations):
        
        gl = m_gene_list(m,n, a)
        selected_path_ids = gl[0]
        #print selected_path_ids
        geneList = gl[1]
        
        # find communities of genes 
        cd_genes = community_detection(geneList,com_method, weights)
        cd_genes_lst = index_to_edge_name(cd_genes)
        
        # keep communities with at least 3 genes
        cd_com_lst = []
        for com in cd_genes_lst: 
            if len(com) >= 3:
                cd_com_lst.append(com)
        
        selected_com_ids = set(selected_path_ids)   # make path ids a set
        success = 0
        detected_com = []
        
        for com in cd_com_lst:                      # for each community, do enrichment
            results = enrichment(com, PATH_GENES, alpha, ALL_GENES)
            total_positive = sum(results[1])
        
            # results returns [sig_tup, signif (T/F), pvals, # signif]
            # sig tup has form (path_id, T/F enriched)
            # results[0][1][1] is most sig enriched tuple getting T/F value 
        
            if results[0][1][1] == True:       # checks if the top result is significant
                most_enriched = results[0][0][0]   # if yes, then set most_enriched to path id
                
            detected = 0
            if most_enriched in selected_path_ids:      # if most_enriched a seeded path, then correct detection 
                detected = 1
                
            total_false_positive += total_positive - detected       # for any given iteration, this will always be 1 or 0 
            
            detected_com.append(most_enriched)
            
        detected_com = set(detected_com)
            
        # number found paths equal to num of seeded pathss
        #print detected_com
        #print selected_com_ids
        #print detected_com == selected_com_ids
        
        if detected_com == selected_com_ids:
            successes += 1
              
    power = float(successes)/iterations 
    avg_false_pos = float(total_false_positive) / iterations 
    return [power, avg_false_pos]
    
# ---------------------------------------------------------------------------
# Performs N experiments of draws m pathways from ontology and prints the nuber
# of false positives and power accounting for all detected paths. 
# ---------------------------------------------------------------------------

 # This is control arm 
def Experiment_partial_pathway(iterations, m, n): 
    ''' Simulation of N iterations of n% of m randomly chosen paths. Returns the 
    average power and false positives'''
    total_false_positive = 0
    successes = 0

    random.seed(123)
    for iteration in range(iterations):
        selected_path_ids = random.sample(range(len(path_genes)), m)
        geneList = []
        
        successfulEnrichment = 0
        for selected_path in selected_path_ids:
            path = path_genes[selected_path]
            path = random.sample(path,len(path)/n)  # choose 1/n% of genes 
            geneList += path
        
        results = enrichment(geneList, path_genes)
            
        for selected_path in selected_path_ids:
            if results[1][selected_path]:
                successfulEnrichment += 1

        if successfulEnrichment == m:
            successes += 1
        
        total_positive = sum(results[1])
        total_false_positive += total_positive - successfulEnrichment
    
    power = float(successes)/iterations 
    avg_false_pos = float(total_false_positive) / iterations 
    return [power, avg_false_pos]
     
# This is experimental arm     
def Experiment_partial_pathway_community(com_method, iterations, m, n, weights): 
    ''' Simulation of N iterations of 1/n% of m randomly chosen paths. Returns the 
    average power and false positives'''
    total_false_positive = 0
    successes = 0
    
    random.seed(123)
    for iteration in range(iterations):
        selected_path_ids = random.sample(range(len(path_genes)), m)
        
        geneList = []
        # create gene list
        for selected_path in selected_path_ids:
            path = path_genes[selected_path]
            path = random.sample(path,len(path)/n)  # choose 1/n% of genes 
            geneList += path
        
        # find communities of genes 
        cd_genes = community_detection(geneList,com_method, weights)
        cd_genes_lst = index_to_edge_name(cd_genes)
        
        # keep communities with at least 3 genes
        cd_com_lst = []
        for com in cd_genes_lst: 
            if len(com) >= 3:
                cd_com_lst.append(com)
        
        selected_com_ids = set(selected_path_ids)   # make path ids a set
        success = 0
        detected_com = set()
        
        for com in cd_com_lst:                      # for each community, do enrichment
            successfulEnrichment = 0
            results = enrichment2(com,path_genes)
            total_positive = sum([results[r][1] for r in range(len(results))])
        
            results2 = []                   # create list of only signif results
            for r in range(len(results)): 
                if results[r][1] == True: 
                    results2.append(results[r])
            
            # path ids of significant paths 
            enrich_path_ids = set([results2[re][0] for re in range(len(results2))])
            
            # set of paths seeded also found
            sig_seeded_paths = selected_com_ids.intersection(enrich_path_ids)
            detected_com = set(detected_com.union(sig_seeded_paths))
            
            total_false_positive += total_positive - len(sig_seeded_paths) 
            
        # number found paths equal to num of seeded pathss
        if detected_com == selected_com_ids:
            successes += 1
              
    power = float(successes)/iterations 
    avg_false_pos = float(total_false_positive) / iterations 
    return [power, avg_false_pos]
    









    
     
    
    

   

















    
    
    

    