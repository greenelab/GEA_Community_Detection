#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 11:30:43 2017

@author: LiaHarrington
"""

import pandas as pd
import random
from enrichment_testing import enrichment

gene_nums = [50, 100, 500]
iterations = 1000

def control(iterations, gene_nums): 
    
    columns = ['num_genes', 'num_enriched']
    results_df = pd.DataFrame(columns=columns)
    
    num_enriched = []
    num_genes = []

    for num in gene_nums: 
        
        for iteration in range(iterations): 
            
            gene_list = random.sample(ALL_GENES, num)
            
            # enrichment has return format of 
            # [sig_tup, signif, pvals, 'Number pathways significant: {0}'.format(sum(signif))]
            # thus to get num significant (after adjustment), we want enrich[1]

            enrich = enrichment(gene_list, PATH_GENES, .05, ALL_GENES)
            
            num_enriched.append(sum(enrich[1]))
            num_genes.append(num)
            
    results_df.num_enriched = num_enriched
    results_df.num_genes = num_genes    
           
    return results_df

control_df = control(iterations, gene_nums)
control_df.to_csv('ctr_all_control', sep = ',')
