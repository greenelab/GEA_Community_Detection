#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 00:21:52 2017

@author: LiaHarrington
"""

execfile('ontology_prep.py')  #  creates ALL_GENES and PATH_GENES

import random
import operator
from scipy.stats import hypergeom
import numpy as np
from community_detection import community_detection, index_to_edge_name
from enrichment_testing import enrichment
 

hgsc = pd.read_csv('entrezid_hgsc.txt', sep = '\t')
hgsc.reset_index(level=0, inplace=True)

hgsc.rename(columns={'index': 'Genes'}, inplace=True)

k4c1 = hgsc.Genes[hgsc.K4Cluster1_Pos==1].append(hgsc.Genes[hgsc.K4Cluster1_Neg==1])
k4c2 = hgsc.Genes[hgsc.K4Cluster2_Pos==1].append(hgsc.Genes[hgsc.K4Cluster2_Neg==1])
k4c3 = hgsc.Genes[hgsc.K4Cluster3_Pos==1].append(hgsc.Genes[hgsc.K4Cluster3_Neg==1])
k4c4 = hgsc.Genes[hgsc.K4Cluster4_Pos==1].append(hgsc.Genes[hgsc.K4Cluster4_Neg==1])

def not_in(gene_list):
    '''finds genes either not in IMP or PID ontology'''
    not_imp = []
    not_all = []
    both = []
    for gene in gene_list: 
        if gene not in IMP_GENES and gene in ALL_GENES: 
            not_imp.append(gene)
        elif gene not in ALL_GENES and gene in IMP_GENES: 
            not_all.append(gene)
        else: 
            both.append(gene)
    return not_imp, not_all, both

def convert_string(gene_list): 
    '''converts the list of integer gene ids to string type'''
    return [str(g) for g in gene_list]

def remove_not_in_IMP(gene_list):
    '''removes all genes in the gene_list that are not in IMP'''
    missing = []
    for gene in gene_list: 
        if gene not in IMP_GENES: 
            missing.append(gene)
    return [gene for gene in gene_list if gene not in missing]
            

all_genes_lst = [k4c1, k4c2, k4c3, k4c4]

def cd_gea_pathways(all_genes_lst, com_method, alpha=.05, min_com_size=3, weights=None): 
    '''loops over all the gene lists using selected community detection method 
        and returns a dictionary with cluster id as key and the community genes
        along with identified significant KEGG pathways'''
    
    path_dict = dict.fromkeys(['k4c1', 'k4c2', 'k4c3', 'k4c4'])
    
    keys = sorted(path_dict.keys())
    
    for i, lst in enumerate(all_genes_lst): 
        gene_list = remove_not_in_IMP(convert_string(lst))
    
        cd_genes = community_detection(gene_list, com_method, weights=None)
        cd_genes_lst = index_to_edge_name(cd_genes)
    
        cd_com_lst = []
        for com in cd_genes_lst:
        # keep communities with at min community size
            if len(com) >= min_com_size:
                cd_com_lst.append(com)
                        
        top_signif_paths = set()
        nonsig_paths = set()
        
        for com in cd_com_lst:
            results = enrichment(com, PATH_GENES, alpha, ALL_GENES)
            top_path = results[0][0][0] # integar id of top path
    
            if results[0][0][1]:  # is the top path significant? (T/F)
                top_signif_paths.add(top_path)
                
            path_names = [PATH_NAMES[j] for j in list(top_signif_paths)]
                
            path_dict[keys[i]] = (path_names, com)
                
    return path_dict
