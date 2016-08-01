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
m number of path and an additional a% added genes. Control involves standard
enrichment analysis whereas Experimental invokes community detection before
enrichment analyses.

"""

execfile('ontology_prep.py')
execfile('community_detection.py')
execfile('enrichment_testing.py')

# Note: Control_most_signif and Experimental_most_signif from Enrichment_Testing

import numpy as np

def run_test(iterations, exp_type, com_method, weights, min_com_size, alpha=.05):
    '''
    Description
    :runs an enrichment test using either control or experimental community
    :detection methods

    Arguments
    :param iterations: number of desired iterations of experiment
    :param exp_type: specifies if should be 'ctr_all' (all significant pathways),
    :'ctr_m' (only top m signif pathways), or 'exp'
    :param com_method: can be 'fastgreedy', 'walktrap', 'infomap,' or
    'multilevel'
    :param weights: Weights can either be IMP WEIGHTS using weights w/o quotes or
    no weights using None.
    :param min_com_size: if using community detection, 3 should be default
    :param alpha: desired alpha level, usually .05 is chosen

    Output
    :returns a list of all the TP/FP rate [[ 0.233], [ 0.18 ], ...] and the
    :combination of number of pathways selected, percentage of pathway chosen,
    :percentage of additional genes in format
    :[(3, 0.3, 0.1), (3, 0.3, 0.25), ...] where the first tuple is 3 paths,
    :using 30% of genes, and 10% additional random genes. The labels are used
    :in the other functions below for formatting output

    '''

    num_rows = len(np.linspace(.3, 1, 5))*len(np.linspace(.1, .7, 5))
    reps = len(range(3, 5))
    tot_rows = num_rows*reps

    results = np.zeros([tot_rows, 1])

    labels = []

    row = 0
    for num_paths in range(3, 5):
        for percent_path in np.linspace(.3, 1, 5):  # from 30% to 100%
            for percent_addit in np.linspace(.1, .7, 5): # from 10% to 70%
                results[row][0] = gsea_performance(iterations,
                                                   num_paths, percent_path,
                                                   percent_addit,
                                                   exp_type, com_method, weights,
                                                   min_com_size, alpha=.05)[4]
                labels.append((num_paths, round(percent_path, 3),
                               round(percent_addit, 3)))
                row += 1

    return [results, labels]

def param_sweep(iterations):
    '''
    Description
    :performs paramter sweep over control and experimental methods selecting

    Arguments
    :param iterations: number of iterations to run each case


    Output
    :matrix of every num_paths, percent_path, and percent_addit combination
    :for every method and the associated
    :(num_paths, percent_path, percent_addit) tuple labels

    '''

    # len of all combinations of num_paths, percent_path, percent_addit
    num_rows = len(np.linspace(.3, 1, 5))*len(range(3, 5))

    # 6 methods (2 control, 4 community)
    num_cols = 6

    comparison_results = np.zeros((num_rows, num_cols))

    ctr_all = run_test(iterations, 'ctr_all', None, None, .05, None)
    ctr_m = run_test(iterations, 'ctr_m', None, None, .05, None)
    fg = run_test(iterations, 'exp', 'fastgreedy', None, .05, 3)
    wt = run_test(iterations, 'exp', 'walktrap', None, .05, 3)
    im = run_test(iterations, 'exp', 'infomap', None, .05, 3)
    ml = run_test(iterations, 'exp', 'multilevel', None, .05, 3)

    labels = ctr_all[1]  # grab labels from first, same for all methods

    # the first row is labels and the other
    param_results = np.concatenate((ctr_all[0], ctr_m[0], fg[0],
                                    wt[0], im[0], ml[0]), axis=1)

    return [param_results, labels]

def param_summary(sweep_object):
    '''
    Description
    :provides summary of which method is top performing over all methods
    :and parameter variations and also provides which method is best at each
    :level of paramter/method combination

    Arguments
    :object created from using param_sweep()

    Output
    :prints dictionary of count of each method and how many times it was
    :best performing method (highest TP/FP ratio), prints the single best
    :method and count, returns for each (num_paths, percent_path, percent_addit)
    :tuple, the associated best performing method

    '''

    param_nums = sweep_object[0]

    param_labels = sweep_object[1]

    location_dict = dict.fromkeys(range(0, 6), 0)

    for row in range(len(param_nums)):
        location = list(param_nums[row]).index(max(param_nums[row]))
        if location in location_dict:
            location_dict[location] += 1
        else:
            location_dict[location] = 0

    # output for location_dict looks like
    # {0: 0, 1: 4, 2: 33, 3: 0, 4: 5, 5: 8}

    method_labels = dict.fromkeys(['ctr_all', 'ctr_m', 'fg', 'wt', 'im', 'ml'])

    # converts location_dict numbered methods to abreviations
    # output looks like: {'ctr_all': 0, 'ctr_m': 4, 'fg': 33, 'im': 0, 'ml': 5, 'wt': 8}
    labeled_locations = {sorted(method_labels.keys())[i]:location_dict[i]
                         for i in range(len(location_dict))}

    maxx = max(location_dict.values())
    max_methods = [(x, y) for x, y in labeled_locations.items() if y == maxx]

    print '-------- Dictionary count of best performing methods --------'

    print labeled_locations

    print '-------- Best performing methods --------'

    print max_methods  # looks like [('fg', 33)]

    param_dict = dict.fromkeys(param_labels)
    combinations = {}
    for row in range(len(param_labels)):
        location = list(param_nums[row]).index(max(param_nums[row]))
        combinations[param_dict.keys()[row]] = location

    # converts all the numbers (0-5) which code for each of the 6 methods
    # to more readible abreviations
    for item in combinations:
        if combinations[item] == 0:
            combinations[item] = 'ctr_all'
        if combinations[item] == 1:
            combinations[item] = 'ctr_m'
        if combinations[item] == 2:
            combinations[item] = 'fg'
        if combinations[item] == 3:
            combinations[item] = 'wt'
        if combinations[item] == 4:
            combinations[item] = 'im'
        if combinations[item] == 5:
            combinations[item] = 'ml'

    print '-------- Best method for each path number and path percentage combo --------'

# sample output for combinations:
#
#    {(3, 0.3, 0.1): 'ml',
#     (3, 0.3, 0.25): 'fg',
#     (3, 0.3, 0.4): 'fg'}
#
# best method for (3, .3, .1) using 3 pathways, 30% of genes, and 10% additional
# genes is 'ml' which is multilevel etc

    return combinations  # returns used so that each paramater level combo
                         # gets own line

def times_better(sweep_object):
    '''
    Description
    :from inital runs, fg appears to be best method so this finds the ratio of
    :TP/FP of the 'fg' to 'ctr_all' methods

    Arguments
    :object created from using param_sweep()

    Output
    :average ratio of the TP/FP for the 'fg' and 'ctr_all' methods over all the
    :iterations and paramter combinations

    '''

    param_nums = sweep_object[0]

    ratios = []

    for combination in param_nums:
        # 'ctr_all' is the first method and 'fg' is the third method
        ratios.append(float(round(combination[2]/combination[0], 2)))

    mean_ratio = np.mean(ratios)

    print 'On average the fastgreedy method is {0} times better than the \
    control using all signif paths.'.format(mean_ratio)
