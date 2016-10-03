# -*- coding: utf-8 -*-
"""
Created on Sun Sep 18 22:39:51 2016

@author: LiaHarrington
"""

execfile('ontology_prep.py')
execfile('community_detection.py')
execfile('enrichment_testing.py')

# Note: Control_most_signif and Experimental_most_signif from Enrichment_Testing

import os
import numpy as np
import matplotlib.pyplot as plt

def data_generation(iterations, num_paths_min, num_paths_max, percent_min,
                    percent_max, addit_min, addit_max, file_name, exp_type=None, com_method=None,
                    weights=None, min_com_size=None, alpha=.05):
    '''
    Description
    Calls gsea_performance() over desired number of iterations for each paramter
    level combination of num_paths, percet, and percent_addit to create data for
    enrichment simulation

    Arguments
    :param iterations: number of desired iterations of experiment
    :param exp_type: must set to False if experimental condition desired
    :param num_paths: number of paths that should be randomly selected
    :param percent_path: value [0,1] of proportion of genes in each selected
    path taken
    :param percent_addit: value [0,1] of proportion of each pathway of
    random extra genes from the ontology that should be added
    :param file_name: desired name of file
    :param com_method: can be 'fastgreedy', 'walktrap', 'infomap,' or
    'multilevel', defaults to None if control condition
    :param weights: Weights can either be IMP weights using WEIGHTS w/o quotes or
    no weights using "NULL". Note n can't be too small for the community detection
    methods because insufficent genes to create network. Defaults to None if control.
    :param min_com_size: minimize number of genes in an acceptable community.
    Defaults to None if control.
    :param alpha: desired alpha level, defaults to 0.05

    Output
    A tab separated file of the data generated
    '''

    # initializes empty vector
    results = np.empty((0, 8))
    labels = []
    conditions = []

    # loops over methods and creates enrichment data
    methods = ['ctr_all', 'ctr_m', 'fastgreedy', 'walktrap', 'infomap', 'multilevel']
    for method in methods:
        if method in ['ctr_all', 'ctr_m']:
            exp_type = method
            com_method = None
        else:
            exp_type = 'exp'
            com_method = method
        for num_paths in range(num_paths_min, num_paths_max):
            for percent_path in np.linspace(percent_min, percent_max, 2):
                for percent_addit in np.linspace(addit_min, addit_max, 2):
                    res = gsea_performance(iterations, num_paths, percent_path, percent_addit,
                                           exp_type=exp_type, com_method=com_method, weights=None,
                                           min_com_size=None, alpha=.05)
                    results = np.concatenate((results, res), axis=0)
                    conditions.extend([method]*iterations)

    data = np.concatenate((np.matrix(conditions).transpose(), results), axis=1)

    pd.DataFrame(data,
                 columns=['Exp_type', 'Iter_num', '# Path', '% Path', '% Addit', 'TP',
                          'FP', 'FN', 'TN']).to_csv('./results/{0}.csv'.format(file_name), sep='\t')
