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

def ctrm_data(iterations, num_paths_min, num_paths_max, percent_min, 
              percent_max, addit_min, addit_max):
                  
#    iterations = 5
#    num_paths_min = 3
#    num_paths_max = 5
#    percent_min = .3
#    percent_max = 1
#    addit_min = 0
#    addit_max = .7 
    
    # initializes empty vector 
    results = np.zeros([1, 5])
    labels = []
    
    for num_paths in range(num_paths_min, num_paths_max):
        for percent_path in np.linspace(percent_min, percent_max, 2):
            for percent_addit in np.linspace(addit_min, addit_max, 2):
                res = gsea_performance(iterations, num_paths, percent_path, percent_addit,
                     exp_type='ctr_m', com_method=None, weights=None,
                     min_com_size=None, alpha=.05)
                results = np.concatenate((results, res), axis = 0)
                labels.append((num_paths, round(percent_path, 3), round(percent_addit, 3)))
        
     # gets rid of inital zero vector row
     results = results[1:, ]
     full_labels = np.matrix([label for label in labels for i in range(iterations)])
        
     data = np.concatenate((full_labels, results), axis = 1) 
     
     pd.DataFrame(data,
                 columns=['# Path', '% Path', '% Addit', 'Iter_num', 'TP', 
                 'FP', 'FN', 'TN']).to_csv('./results/data.csv')

          