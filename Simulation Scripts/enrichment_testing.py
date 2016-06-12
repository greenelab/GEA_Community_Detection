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
    : Performs gene enrichment analysis

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
# Randomly selects n% of m pathways and a% (using floor function) additional random
# genes and then writes to text file. m is any whole number, n and a are any number
# between 0 and 1.
# ---------------------------------------------------------------------------

def m_gene_list(num_paths, percent_path, percent_addit):
    '''
    Description
    :creates a gene list to input into the Enrichment function

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
    :Writes the current gene list to the specified file name.

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
# Performs N experiments of draws m pathways from ontology and prints the nuber
# of false positives and power accounting only the most significanlty
# enriched path.
# ---------------------------------------------------------------------------

# Control arm
def control_most_signif(iterations, num_paths, percent_path,
                        percent_addit, alpha):
    '''
    Description
    :Simulation of N iterations of m chosen paths using n% of each path with a%
    additional random genes.

    Arguments
    :param iterations: number of desired iterations of experiment
    :param num_paths: number of paths that should be randomly selected
    :param percent_path: value [0,1] of proportion of genes in each selected
    path taken
    :param percent_addit: value [0,1] of proportion of each pathway of
    random extra genes from the ontology that should be added
    :param alpha: desired alpha level, usually .05 is chosen

    Output
    :Returns the average power and false positives. Power and
    success here is defined for if the most signif paths is the seeded path.

    '''

    random.seed(123)
    total_false_positive = 0
    successes = 0

    for iteration in range(iterations):
        # Randomly select a gene list
        gl_object = m_gene_list(num_paths, percent_path, percent_addit)
        selected_path_ids = gl_object[0]
        gene_list = gl_object[1]

        # results returns [sig_tup, signif (T/F), pvals, # signif]
        results = enrichment(gene_list, PATH_GENES, alpha, ALL_GENES)

        # get the m most signif results
        relevant_results = results[0][0:num_paths]

        # extract path number
        sig_paths = [relevant_results[i][0] for i in range(len(relevant_results))]

        # check if top m results equal to the m seeded paths
        if sum([relevant_results[i][1] for i in range(len(relevant_results))]) == num_paths \
        and set(sig_paths) == set(selected_path_ids):
            successes += 1

        total_positive = sum(results[1])

        # subtract out corectly detected paths from all detected

        correct_sig_paths = len((set(selected_path_ids).intersection(set(sig_paths))))

        total_false_positive += total_positive - correct_sig_paths

    power = float(successes)/iterations
    avg_false_pos = float(total_false_positive) / iterations
    print [power, avg_false_pos]

# Experimental arm
def experimental_most_signif(iterations, com_method, num_paths, percent_path,
                             percent_addit, weights, alpha, min_com_size):
    '''
    Description
    :Simulation of N iterations of m chosen paths using n% of each path with a%
    additional random genes using community detection

    Arguments
    :param iterations: number of desired iterations of experiment
    :param com_method: can be 'fastgreedy', 'walktrap', 'infomap,' or
    'multilevel'
    :param num_paths: number of paths that should be randomly selected
    :param percent_path: value [0,1] of proportion of genes in each selected
    path taken
    :param percent_addit: value [0,1] of proportion of each pathway of
    random extra genes from the ontology that should be added
    :param weights: Weights can either be IMP weights using WEIGHTS w/o quotes or
    no weights using "NULL". Note n can't be too small for the community detection
    methods because insufficent genes to create network.
    :param alpha: desired alpha level, usually .05 is chosen
    :param min_com_size: minimize number of genes in an acceptable community

    Output
    :Returns the average power and false positives. Power and
    success here is defined for if the most signif paths is the seeded path.

    '''

    random.seed(123)
    total_false_positive = 0
    successes = 0

    for iteration in range(iterations):
        # Randomly select a gene list
        gl_object = m_gene_list(num_paths, percent_path, percent_addit)
        selected_path_ids = gl_object[0]
        gene_list = gl_object[1]

        # find communities of genes
        cd_genes = community_detection(gene_list, com_method, weights)
        cd_genes_lst = index_to_edge_name(cd_genes)

        # keep communities with at min community size
        cd_com_lst = []
        for com in cd_genes_lst:
            if len(com) >= min_com_size:
                cd_com_lst.append(com)

        selected_com_ids = set(selected_path_ids)  # make path ids a set
        detected_com = []

        for com in cd_com_lst:  # for each community, do enrichment
            results = enrichment(com, PATH_GENES, alpha, ALL_GENES)
            total_positive = sum(results[1])  # sum of T/F list

            # results returns [sig_tup, signif (T/F), pvals, # signif]
            # sig tup has form (path_id, T/F enriched)
            # results[0][1][1] is most sig enriched tuple getting T/F value

            if results[0][1][1]:  # checks if the top result is significant
                most_enriched = results[0][0][0]  # if yes, then set most_enriched to path id
            else:
                most_enriched = 'None'

            detected = 0
            # if most_enriched a seeded path, then correct detection
            if most_enriched in selected_path_ids:
                detected = 1
            else:
                detected = 0

            # for any given iteration, this will always be 1 or 0
            total_false_positive += total_positive - detected

            detected_com.append(most_enriched)

        detected_com = set(detected_com)

        # number found paths equal to num of seeded pathss
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

def control_all_signif(iterations, num_paths, percent_path, percent_addit,
                       alpha):
    '''
    Description
    :Simulation of N iterations of m chosen paths using n% of each path with a%
    additional random genes.

    Arguments
    :param iterations: number of desired iterations of experiment
    :param num_paths: number of paths that should be randomly selected
    :param percent_path: value [0,1] of proportion of genes in each selected
    path taken
    :param percent_addit: value [0,1] of proportion of each pathway of
    random extra genes from the ontology that should be added
    :param alpha: desired alpha level, usually .05 is chosen

    Output
    :Returns the average power and false positives. Power and
    success here is defined for if the m signif paths returned are the
    seeded paths.

    '''

    total_false_positive = 0
    successes = 0

    random.seed(123)
    for iteration in range(iterations):
        # Randomly select a gene list
        gl_object = m_gene_list(num_paths, percent_path, percent_addit)
        selected_path_ids = gl_object[0]
        gene_list = gl_object[1]

        successful_enrichment = 0

        results = enrichment(gene_list, PATH_GENES, alpha, ALL_GENES)

        for selected_path in selected_path_ids:
            if results[1][selected_path]:
                successful_enrichment += 1

        if successful_enrichment == num_paths:
            successes += 1

        total_positive = sum(results[1])
        total_false_positive += total_positive - successful_enrichment

    power = float(successes)/iterations
    avg_false_pos = float(total_false_positive) / iterations
    return [power, avg_false_pos]

# This is experimental arm
def experimental_all_signif(iterations, com_method, num_paths, percent_path,
                            percent_addit, weights, alpha, min_com_size):

    '''
    Description
    :Simulation of N iterations of m chosen paths using n% of each path with a%
    additional random genes.

    Arguments
    :param iterations: number of desired iterations of experiment
    :param com_method: can be 'fastgreedy', 'walktrap', 'infomap,' or
    'multilevel'
    :param num_paths: number of paths that should be randomly selected
    :param percent_path: value [0,1] of proportion of genes in each selected
    path taken
    :param percent_addit: value [0,1] of proportion of each pathway of
    random extra genes from the ontology that should be added
    :param weights: Weights can either be IMP weights using WEIGHTS w/o quotes or
    no weights using "NULL". Note n can't be too small for the community detection
    :param alpha: desired alpha level, usually .05 is chosen
    :param min_com_size: minimize number of genes in an acceptable community

    Output
    :Returns the average power and false positives. Power and
    success here is defined for if the m signif paths returned are the
    seeded paths.

    '''

    total_false_positive = 0
    successes = 0

    random.seed(123)
    for iteration in range(iterations):
        # Randomly select a gene list
        gl_object = m_gene_list(num_paths, percent_path, percent_addit)
        selected_path_ids = gl_object[0]
        gene_list = gl_object[1]

        # find communities of genes
        cd_genes = community_detection(gene_list, com_method, weights)
        cd_genes_lst = index_to_edge_name(cd_genes)

        # keep communities with at least min_com_size genes
        cd_com_lst = []
        for com in cd_genes_lst:
            if len(com) >= min_com_size:
                cd_com_lst.append(com)

        selected_com_ids = set(selected_path_ids)   # make path ids a set
        detected_com = set()

        for com in cd_com_lst:                      # for each community, do enrichment
            results = enrichment(com, PATH_GENES, alpha, ALL_GENES)
            total_positive = sum(results[1])

            results2 = []                   # create list of only signif results
            for result in range(len(results)):
                if results[0][result][1]:
                    results2.append(results[0][result])

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
    print [power, avg_false_pos]
