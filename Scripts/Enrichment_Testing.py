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

execfile('ontology_prep.py')  #  creates ALL_GENES and PATH_GENES
execfile('community_detection.py')  # community detection functions used in the
                                    # experimental arm

import random
import operator
from scipy.stats import hypergeom

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

# ---------------------------------------------------------------------------
#                       GOAL: MOST SIGNIFICANT PATHWAY
# ---------------------------------------------------------------------------
#                               Control Arm
# Performs N experiments drawing m pathways using n% of the genes in each pathway
# with an additional a% random genes from ontology and returns the
# power, false positive rate, total true positive, total false positive, and
# TP to FP ratio. The output will look something like this:

#['Power: 0.683',
# 'False positive rate: 0.00109.',
# 'Total true positive: 41.0',
# 'Total false positive: 4.0',
# 'TP to FP ratio: 10.25']


#For the control arm, no community detection is performed and one can simply
# fill in the first four paramters as the other paramters concerining
# the experimental arm default to None. When performing the experimental arm,
# be sure to compelte all default paramters except alpha, unless a different
# level is desired.
# ---------------------------------------------------------------------------

def gsea_performance(iterations, num_paths, percent_path, percent_addit,
                     exp_type='control_all', com_method=None, weights=None,
                     min_com_size=None, alpha=.05):
    '''
    Description
    Simulation of N iterations of m chosen paths using n% of each path with a%
    additional random genes using community detection

    Arguments
    :param iterations: number of desired iterations of experiment
    :param control: True/False and indicates if control or experimental condition
    :must set to False if experimental condition desired
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
    :Returns the average power and false positives. Power and
    success here is defined for if the most signif paths is the seeded path.

    '''

    random.seed(123)

    true_pos_nums = []
    false_pos_nums = []
    false_neg_nums = []
    true_neg_nums = []

    for iteration in range(iterations):
        # Randomly select a gene list
        gl_object = m_gene_list(num_paths, percent_path, percent_addit)
        selected_path_ids = gl_object[0]
        set_selected_pathids = set(selected_path_ids)
        gene_list = gl_object[1]

        if exp_type == 'ctr_all': # get all signif paths
            # results returns [sig_tup, signif (T/F), pvals, # signif]
            results = enrichment(gene_list, PATH_GENES, alpha, ALL_GENES)

            top_signif_paths = set([results[0][i][0] for i in range(len(results[0])) if
                           results[0][i][1]])

            non_top_paths = set(range(len(PATH_GENES))).difference(top_signif_paths)


        elif exp_type == 'ctr_m':  # get signif paths in top m

            results = enrichment(gene_list, PATH_GENES, alpha, ALL_GENES)

            # the top m significant paths
            relevant_results = results[0][0:num_paths]

            # list of top m paths
            top_m_paths = [relevant_results[i] for i in range(len(relevant_results))]

            # signif paths in the top m paths
            top_signif_paths = set([top_m_paths[i][0] for i in range(len(relevant_results))
                                if top_m_paths[i][1]])

            # paths that are not in the top m and signif
            non_top_paths = set(range(len(PATH_GENES))).difference(top_signif_paths)

        elif exp_type == 'exp':
            # find communities of genes
            cd_genes = community_detection(gene_list, com_method, weights)
            cd_genes_lst = index_to_edge_name(cd_genes)

            cd_com_lst = []
            for com in cd_genes_lst:
            # keep communities with at min community size
                if len(com) >= min_com_size:
                    cd_com_lst.append(com)

            top_signif_paths = set()
            nonsig_paths = set()
            for com in cd_com_lst:  # for each community, do enrichment
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

        true_pos_nums.append(float(len(true_pos)))
        false_pos_nums.append(float(len(false_pos)))
        false_neg_nums.append(float(len(false_neg)))
        true_neg_nums.append(float(len(true_neg)))

    power = float(sum(true_pos_nums))/(sum(true_pos_nums)+sum(false_neg_nums))
    false_pos_rate = float(sum(false_pos_nums)) / (sum(false_pos_nums)+sum(true_neg_nums))

#    return ['Power: {0}'.format(round(power, 3)),
#            'False positive rate: {0}.'.format(round(false_pos_rate, 5)),
#            'Total true positive: {0}'.format(sum(true_pos_nums)),
#            'Total false positive: {0}'.format(sum(false_pos_nums)),
#            'TP/FP: {0}'.format(round(sum(true_pos_nums)/sum(false_pos_nums), 3))]

    return [round(power, 3),
            round(false_pos_rate, 5),
            sum(true_pos_nums),
            sum(false_pos_nums),
            round(sum(true_pos_nums)/sum(false_pos_nums), 3)]
