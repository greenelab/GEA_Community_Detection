# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 12:51:01 2016

@author: LiaHarrington

Usage:
Imported by Enrichment_Main. Requires Ontology_prep and Community_Detection.

Description:
# Completes enrichment analysis of randonly chosen pathways from the an
# ontology N times with and without community detection

"""

import random
import operator
from scipy.stats import hypergeom
import numpy as np
import pickle 
import pandas as pd
from community_detection import community_detection, index_to_edge_name
import os

all_genes_file = os.path.join('Data', 'PID_all_genes.pkl')
ALL_GENES = pickle.load(open(all_genes_file, 'r'))

path_genes_file = os.path.join('Data', 'PID_path_genes.pkl')
PATH_GENES = pickle.load(open(path_genes_file, 'r'))

def enrichment(gene_list, ontology, alpha, all_genes):
    '''
    Description
    Performs gene enrichment analysis

    Arguments
    :param geneList: a list of genes (can be created using m_gene_list())
    :param ontology: ontology created from Ontology_prep
    :param alpha: desired alpha level, usually .05 is chosen
    :param all_genes: unlisted by path list of all genes in the ontology

    Internal Arguments
    :internal: len_all_genes = # of distinct genes in ontology
    :internal: len_path = # of genes in ith pathway
    :internal: len_gene_list = # of genes in genelist
    internal: intersection = # of genes in genelist that are also in pathway.

    Output
    :outpug sig_tup: a sorted list of tuples based on pvalue of the pathway
    number and associated pvale in form (path_id, pval)
    :output signif: list of T/F value of if pval less than set alpha
    :output pvals: pvals of paths. Note, pvals are not sorted here. So
    sig_tup[ix]!= pvals[idx].
    :output number signif: a string output stating number of statistically
    significant pathways detected

    '''

    len_all_genes = len(all_genes)
    gene_list = set(gene_list)
    pvals = []

    for path in range(len(ontology)):
        pathway = set(ontology[path])
        len_path = len(pathway)
        len_gene_list = len(gene_list)
        intersection = len(pathway.intersection(gene_list))
        pvals.append(hypergeom.sf(intersection-1, len_all_genes, len_path,
                                  len_gene_list))

    p_tups = [(p, pvals[p]) for p in range(len(pvals))]
    p_tups.sort(key=operator.itemgetter(1))

    bonf_alpha = alpha/len(ontology)
    signif = [pval < bonf_alpha for pval in pvals]

    sig_tup = [(p_tups[i][0], p_tups[i][1] < bonf_alpha) for i in range(len(p_tups))]

    return [sig_tup, signif, pvals, 'Number pathways significant: {0}'.format(sum(signif))]

# ---------------------------------------------------------------------------
# Randomly selects n% of m pathways and a% (using floor function) additional
# random genes and then writes to text file. m is any whole number, n and a
# are any number between 0 and 1.
# ---------------------------------------------------------------------------

def m_gene_list(num_paths, percent_path, percent_addit):
    '''
    Description
    creates a gene list to input into the Enrichment function

    Arguments
    :param num_paths: integer value of number of desired chosen paths
    :param: percent_path: value [0,1] to indicate proportion of path selected
    :param percent_addit: value [0,1] to proportion of extra genes to add

    Output
    :output selected_path_ids: ids of the paths selected
    :output geneList: the desired gene list with specified charateristics
    '''

    selected_path_ids = random.sample(range(len(PATH_GENES)), num_paths)
    gene_list = []

    for selected_path in selected_path_ids:
        path = PATH_GENES[selected_path]
        path = random.sample(path, int(float(percent_path)*len(path)))
        gene_list += path
    gene_list += random.sample(ALL_GENES,
                               int(float(percent_addit)*len(gene_list)))

    return [selected_path_ids, gene_list]

def write_gene_list(gene_list, text_fh):
    '''
    Description
    Writes the current gene list to the specified file name.

    Arguments
    : param geneList: an input list of gene names
    : param text_fh: the file name to save

    Output
    :Writes a genelist to file
    '''

    with open(text_fh, 'w') as text_file:
        for gene in gene_list:
            text_file.write(str(gene)+'\n')


def gea_performance(iterations, exp_type, num_paths, percent_path, percent_addit,
                    ctr_method=None, com_method=None, weights=None,
                    min_com_size=None, alpha=.05):
    '''
    Description
    Simulation of N iterations of m chosen paths using n% of each path with a%
    additional random genes using community detection

    Arguments
    :param iterations: number of desired iterations of experiment
    :param method: 'ctr_m', 'ctr_all', 'exp'
    :param com_method: can be 'fastgreedy', 'walktrap', 'infomap,' or
    'multilevel', defaults to None if control condition
    :param num_paths: number of paths that should be randomly selected
    :param percent_path: value [0,1] of proportion of genes in each selected
    path taken
    :param percent_addit: value [0,1] of proportion of each pathway of
    random extra genes from the ontology that should be added
    :param weights: Weights can either be IMP weights using WEIGHTS w/o quotes or
    no weights using "NULL". Note n can't be too small for the community detection
    methods because insufficent genes to create network. Defaults to None if control.
    :param min_com_size: minimize number of genes in an acceptable community.
    Defaults to None if control.
    :param alpha: desired alpha level, defaults to 0.05

    Output
    :Returns true positive, false posotive, true negative, and false negative
    for each parameter level combination for desired experiment type and
    community detection method

    '''
    
    results_columns = ['iter_num', 'method', 'num_paths', 'percent_path', 'percent_addit',
                       'true_positive', 'false_positive', 'true_negative',
                       'false_negative']

    summary_results_df = pd.DataFrame(columns=results_columns, index=range(iterations))

    for iteration in range(iterations):
        
        selected_path_ids, gene_list = m_gene_list(num_paths, percent_path, percent_addit)
        set_selected_pathids = set(selected_path_ids)

        if exp_type == 'exp':
            cd_genes = community_detection(gene_list, com_method, weights=None)
            cd_genes_lst = index_to_edge_name(cd_genes)

            cd_com_lst = []
            for com in cd_genes_lst:
            # keep communities with at min community size
                if len(com) >= min_com_size:
                    cd_com_lst.append(com)
        else:
            cd_com_lst = [gene_list]

        if exp_type == 'ctr': 
            top_signif_paths = set()
            nonsig_paths = set()
            for com in cd_com_lst:
                results = enrichment(com, PATH_GENES, alpha, ALL_GENES)
            
    
                if ctr_method == 'ctr_m':
                     # the top m significant paths
                    relevant_results = results[0][0:num_paths]
                    # list of top m paths
                    top_m_paths = [relevant_results[i] for i in range(len(relevant_results))]
    
                    # signif paths in the top m paths
                    top_signif_paths = set([top_m_paths[i][0] for i in range(len(relevant_results))
                                            if top_m_paths[i][1]])
    
                    # paths that are not in the top m and signif
                    non_top_paths = set(range(len(PATH_GENES))).difference(top_signif_paths)
    
    
                elif ctr_method == 'ctr_all':
                    top_signif_paths = set([results[0][i][0] for i in range(len(results[0])) if
                                            results[0][i][1]])
    
                    non_top_paths = set(range(len(PATH_GENES))).difference(top_signif_paths)
        
        elif exp_type == 'exp':  
            top_signif_paths = set()
            nonsig_paths = set()
            for com in cd_com_lst:
                results = enrichment(com, PATH_GENES, alpha, ALL_GENES)
                top_path = results[0][0][0] # integar id of top path

                if results[0][0][1]:  # is the top path significant? (T/F)
                    top_signif_paths.add(top_path)

                ns_paths = set.difference(set(range(len(PATH_GENES))), set([top_path]))

                nonsig_paths = nonsig_paths.union(ns_paths)
            
            non_top_paths = nonsig_paths.difference(top_signif_paths)

        true_pos = set_selected_pathids.intersection(top_signif_paths)
    
        false_pos = top_signif_paths.difference(set_selected_pathids)
    
        false_neg = non_top_paths.intersection(set_selected_pathids)
    
        true_neg = set(range(len(PATH_GENES))).difference(set.union(true_pos,
                                                                    false_pos, false_neg))

        iter_num = np.matrix(range(iterations)).transpose()
    
        summary_results_df.true_positive[iteration] = float(len(true_pos))
        summary_results_df.false_positive[iteration] = float(len(false_pos))
        summary_results_df.false_negative[iteration] = float(len(false_neg))
        summary_results_df.true_negative[iteration] = float(len(true_neg))
    
        all_num_paths = [num_paths for i in range(iterations)]
        all_percent_path = [round(percent_path, 3) for i in range(iterations)]
        all_percent_addit = [round(percent_addit, 3) for i in range(iterations)]
    
        if ctr_method in ['ctr_all', 'ctr_m']:
            summary_results_df = summary_results_df.assign(method = ctr_method)
        else:
            summary_results_df = summary_results_df.assign(method = com_method)
    
    summary_results_df = summary_results_df.assign(num_paths = all_num_paths)
    summary_results_df = summary_results_df.assign(percent_path = all_percent_path)
    summary_results_df = summary_results_df.assign(percent_addit = all_percent_addit)
    summary_results_df = summary_results_df.assign(iter_num = iter_num) 

    return summary_results_df
