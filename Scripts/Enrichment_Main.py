# -*- coding: utf-8 -*-
"""
Created on Sun Sep 18 22:39:51 2016

@author: LiaHarrington

Usage:
Imports functions and constructs from ontology_prep, community_detection
and enrichment_testing

Description:
Creates simulation data over the two controls and four experiemntal conditions
over the various # paths, % path, and % addit gene combinations as a tsv file
"""

# Note: Control_most_signif and Experimental_most_signif from Enrichment_Testing

import numpy as np
import pandas as pd
import sys
from enrichment_testing import gea_performance

def data_generation(iterations, num_paths_min, num_paths_max, percent_min,
                    percent_max, addit_min, addit_max, file_name,
                    weights=None, min_com_size=None, alpha=.05):
    '''
    Description
    Calls gsea_performance() over desired number of iterations for each paramter
    level combination of num_paths, percent, and percent_addit to create data for
    enrichment simulation

    Arguments
    :param iterations: number of desired iterations of experiment
    :param num_paths: number of paths that should be randomly selected
    :param percent_path: value [0,1] of proportion of genes in each selected
    path taken
    :param percent_addit: value [0,1] of proportion of each pathway of
    random extra genes from the ontology that should be added
    :param file_name: desired name of file
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
    results = np.empty((0, 9))

    # loops over methods and creates enrichment data
    methods = ['ctr_all', 'ctr_m', 'fastgreedy', 'walktrap', 'infomap', 'multilevel']
    for method in methods:
        if method in ['ctr_all', 'ctr_m']:
            ctr_method = method
            exp_type = 'ctr'
            com_method = None
        else:
            com_method = method
            exp_type = 'exp' 
            ctr_method = None

        for num_paths in range(num_paths_min, num_paths_max):
            for percent_path in np.linspace(percent_min, percent_max, 5):
                for percent_addit in np.linspace(addit_min, addit_max, 5):
                    res = gea_performance(iterations, exp_type, num_paths, percent_path, percent_addit,
                                          ctr_method=method, com_method=com_method, weights=None,
                                          min_com_size=None, alpha=.05)
                    results = np.concatenate((results, res), axis=0)

    results_columns = ['iter_num', 'method', 'num_paths', 'percent_path', 'percent_addit',
                       'true_positive', 'false_positive', 'true_negative',
                       'false_negative']

    results_df = pd.DataFrame(results, columns=results_columns)

    results_df.to_csv('./Data/{0}'.format(file_name), index=False, sep='\t')
    
if __name__ == '__main__': 
    iterations = sys.argv[1]
    num_paths_min = sys.argv[2]
    num_paths_max = sys.argv[3]
    percent_min = sys.argv[4]
    percent_max = sys.argv[5]
    addit_min = sys.argv[6]
    addit_max = sys.argv[7]
    file_name = sys.argv[8]
    
    data_generation(int(iterations), int(num_paths_min), int(num_paths_max), 
                    float(percent_min), float(percent_max), float(addit_min), 
                         float(addit_max), str(file_name))
    