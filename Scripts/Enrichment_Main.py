# -*- coding: utf-8 -*-
"""
Created on Wed May 11 00:19:58 2016

@author: LiaHarrington

Usage:

Used to call relevant functions from Enrichment_Testing to perform
experiment. Note, requires functions and variables created from Ontology_prep
and Community_Detection

Description:

Runs simulations of of specified iterations using a gene list with n% of
m number of path and an additional a% added genes. Control involves usual
enrichment analysis whereas Experimental invokes community detection before
enrichment analyses.

"""

execfile('ontology_prep.py')
execfile('community_detection.py')
execfile('enrichment_esting.py')

# Note: Control_most_signif and Experimental_most_signif from Enrichment_Testing

import numpy as np

def run_control(iterations, percent_addit, alpha):
    '''
    Description
    :Runs paramter sweep of control arm of gene enrichment experiment

    Arguments
    :param iterations: number of desired iterations of experiment
    :param percent_addit: value [0,1] of proportion of each pathway of
    random extra genes from the ontology that should be added
    :param alpha: desired alpha level, usually .05 is chosen

    Internal Arguments:
    :num_paths (implicit and created on line 52)
    :percent paths (implicit and created on line 53)

    Output
    :returns (num_paths, percent_path): [power, avg_false_pos] for each level
    of num_paths and percent_path

    '''

    results = {}
    for num_paths in range(1, 11):
        for percent_path in np.linspace(.3, 1, 5):
            results[(num_paths, round(percent_path, 2))] = control_most_signif(iterations,
                    num_paths, percent_path, percent_addit, alpha)
    return results

def run_experimental(iterations, com_method, percent_addit,
                     weights, alpha, min_com_size):
    '''
    Description
    :Runs paramter sweep of experimental arm of gene enrichment experiment
    using selected community detection method

    Arguments
    :param iterations: number of desired iterations of experiment
    :param com_method: can be 'fastgreedy', 'walktrap', 'infomap,' or
    'multilevel'
    :param percent_addit: value [0,1] of proportion of each pathway of
    random extra genes from the ontology that should be added
    :param weights: Weights can either be IMP WEIGHTS using weights w/o quotes or
    no weights using "NULL". Note n can't be too small for the community detection
    methods because insufficent genes to create network.
    :param alpha: desired alpha level, usually .05 is chosen

    Internal Arguments:
    :num_paths (implicit and created on line 87)
    :percent paths (implicit and created on line 88)

    Output
    :returns (num_paths, percent_path): [power, avg_false_pos] for each level
    of num_paths and percent_path

    '''

    results = {}
    for num_paths in range(1, 11):
        for percent_path in np.linspace(.3, 1, 5):
            results[(num_paths, round(percent_path, 2))] = experimental_most_signif(iterations,
                    com_method, num_paths, percent_path, percent_addit,
                    weights, alpha, min_com_size)
    return results
