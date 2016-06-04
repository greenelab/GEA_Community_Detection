# -*- coding: utf-8 -*-
"""
Created on Wed May 11 00:19:58 2016

@author: LiaHarrington
"""
execfile('Enrichment_prep.py')
execfile('community_net2.py')
execfile('Enrichment.py')

import numpy as np 

def run_control(iterations, n, a, bonf):
    results = {}
    for m in range(1,11): 
        for n in np.linspace(.3,n,5): 
            results[(m,n)]  = Control_most_signif(iterations, m, n, a, bonf)
    return results
        

def run_experimental(iterations,com_method, n, a, weights, bonf): 
    '''Weights can either be IMP weights using wegiths w/o quotes or 
    no weights using "NULL". Note n can't be too small for the community detection
    methods because insufficent genes to create network.'''
    results = {}
    for m in range(1,11): 
        for n in np.linspace(.3,n,5):
            results[(m,n)] = Experimental_most_signif(com_method,iterations,m,n,a, weights, bonf)
    return results
    
    