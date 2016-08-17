# -*- coding: utf-8 -*-
"""
Created on Wed May 11 00:19:58 2016

@author: LiaHarrington

Usage:

Used to call relevant functions from Enrichment_Testing to perform
experiment. Note, requires functions and variables created from Ontology_prep
and Community_Detection

Description:

Runs simulations of specified iterations using a gene list with n% of
m number of path and an additional a% added genes. Control involves standard
enrichment analysis whereas Experimental invokes community detection before
enrichment analyses.

"""

execfile('./Scripts/ontology_prep.py')
execfile('./Scripts/community_detection.py')
execfile('./Scripts/enrichment_testing.py')

# Note: Control_most_signif and Experimental_most_signif from Enrichment_Testing

import os
import numpy as np
import matplotlib.pyplot as plt

def param_sweep(iterations, num_paths_min, num_paths_max, percent_min,
                percent_max, addit_min, addit_max):
    '''
    Description
    performs paramter sweep over control and experimental methods

    Arguments
    :param iterations: number of iterations to run each case
    :param num_paths_min: min number of paths, cannot be zero
    :param num_paths_max: max number of paths, cannot be greater than length
    :of ontology
    :param percent_min: min percent of genes in each pathway taken, lowest
    :reccomended is .3 for community detection construction
    ::param percent_max: max percent of genes in each pathway taken, can be up
    :to 1
    :param addit_min: min percent of additional random genes added
    :param addit_max: max percent of additional random genes added, use intution,
    :but more than .7 may become less realistic

    Output
    :matrix of every num_paths, percent_path, and percent_addit combination
    :for every method and the associated
    :(num_paths, percent_path, percent_addit) tuple labels

    '''

    # 5 chosen because unlikely want more than 5 possible percentage levels,
    # but this can be changed
    num_rows = len(np.linspace(percent_min,
                               percent_max, 5))*len(np.linspace(addit_min, addit_max, 5))

    reps = len(range(num_paths_min, num_paths_max))

    tot_rows = num_rows*reps

    ctr_all = np.zeros([tot_rows, 1])
    ctr_m = np.zeros([tot_rows, 1])
    fg = np.zeros([tot_rows, 1])
    wt = np.zeros([tot_rows, 1])
    im = np.zeros([tot_rows, 1])
    ml = np.zeros([tot_rows, 1])

    labels = []
    row = 0
    for num_paths in range(num_paths_min, num_paths_max):
        for percent_path in np.linspace(percent_min, percent_max, 5):
            for percent_addit in np.linspace(addit_min, addit_max, 5):
                ctr_all[row] = gsea_performance(iterations, num_paths,
                                                percent_path, percent_addit,
                                                'ctr_all', None, None, None, .05)[4]
                ctr_m[row] = gsea_performance(iterations, num_paths, percent_path,
                                              percent_addit, 'ctr_m', None, None, None, .05)[4]
                fg[row] = gsea_performance(iterations, num_paths, percent_path,
                                           percent_addit, 'exp', 'fastgreedy', None, 3, .05)[4]
                wt[row] = gsea_performance(iterations, num_paths, percent_path,
                                           percent_addit, 'exp', 'walktrap', None, 3, .05)[4]
                im[row] = gsea_performance(iterations, num_paths, percent_path,
                                           percent_addit, 'exp', 'infomap', None, .05)[4]
                ml[row] = gsea_performance(iterations, num_paths, percent_path,
                                           percent_addit, 'exp', 'multilevel', None, 3, .05)[4]
                labels.append((num_paths, round(percent_path, 3),
                               round(percent_addit, 3)))
                row += 1

    # the first row is labels and the other
    param_results = np.concatenate((ctr_all, ctr_m, fg,
                                    wt, im, ml), axis=1)

    return [param_results, labels]

def param_summary(sweep_object):
    '''
    Description
    provides summary of which method is top performing over all methods
    and parameter variations and also provides which method is best at each
    level of paramter/method combination

    Arguments
    :param sweep_object: object created from using param_sweep()

    Output
    :prints dictionary of count of each method and how many times it was
    :best performing method (highest TP/FP ratio), prints the single best
    :method and count, returns for each (num_paths, percent_path, percent_addit)
    :tuple, the associated best performing method. Additionally, plots bar graph
    :of count of when each method had highest TP/FP ratio and exports dictionary
    :data as a tab deliminated csv

    '''

    param_nums, param_labels = sweep_object

    location_dict = dict.fromkeys(range(0, 6), 0)

    for row in range(len(param_nums)):
        location = list(param_nums[row]).index(max(param_nums[row]))
        if location in location_dict:
            location_dict[location] += 1
        else:
            location_dict[location] = 0

    # output for location_dict looks like
    # {0: 0, 1: 4, 2: 33, 3: 0, 4: 5, 5: 8}

    method_labels = ['ctr_all', 'ctr_m', 'fg', 'wt', 'im', 'ml']

    # converts location_dict numbered methods to abreviations
    # output looks like: {'ctr_all': 0, 'ctr_m': 4, 'fg': 33, 'im': 0, 'ml': 5, 'wt': 8}
    labeled_locations = {sorted(dict.fromkeys(method_labels).keys())[i]:location_dict[i]
                         for i in range(len(location_dict))}

    maxx = max(location_dict.values())
    max_methods = [(x, y) for x, y in labeled_locations.items() if y == maxx]

    print '-------- Dictionary count of best performing methods --------'

    print labeled_locations

    # save dictionary as dataframe and output as csv
    try:
        os.mkdir('results')
    except:
        print 'Results folder already exists'

    pd.DataFrame(labeled_locations.items(),
                 columns=['Method', 'Frequency']).to_csv('./results/labeled_locations.csv')

    # make barplot of methods and TP/FP ratio count
    plt.bar(range(len(labeled_locations)), labeled_locations.values(),
            align='center', color='w')
    plt.xticks(range(len(labeled_locations)), labeled_locations.keys())
    plt.ylabel('Frequency best TP/FP ratio')
    plt.xlabel('Method')
    plt.show()

    print '-------- Best performing methods --------'

    print max_methods  # looks like [('fg', 33)]

    combinations = []
    for  row in range(len(param_labels)):
        # find where max TP/FP ratio is
        location = list(param_nums[row]).index(max(param_nums[row]))
        # what method is location associated with?
        identifier = method_labels[location]
        # set method to paramter level combo
        combinations.append((param_labels[row], identifier))

    # save dictionary as dataframe and output as csv
    pd.DataFrame(combinations,
                 columns=['Combination', 'Method']).to_csv('./results/combinations.csv')

    print '-------- Best method for each path number and path percentage combo --------'

#[((3, 0.3, 0.1), 'ml'),
# ((3, 0.3, 0.25), 'im'),
# ((3, 0.3, 0.4), 'fg'),...]
# best method for (3, .3, .1) using 3 pathways, 30% of genes, and 10% additional
# genes is 'ml,' which is multilevel etc

    return combinations  # returns used so that each paramater level combo
                         # gets own line

def rate_ratio(sweep_object, reference_method):
    '''
    Description
    from chosen reference_method, finds the ratio TP/FP ratio compared to
    all other methods including itself

    Arguments
    :param sweep_object: object created from using param_sweep()
    :param reference_method: method that is the reference for the other methods
    :to be compared against. can be 'ctr_all', 'ctr_m', 'fg', 'wt', 'im', 'ml'

    Output
    :average ratio of the TP/FP for the 'fg' and 'ctr_all' methods over all the
    :iterations and paramter combinations

    '''

    param_nums = sweep_object[0]

    names = ['ctr_all', 'ctr_m', 'fg', 'wt', 'im', 'ml']

    reference_location = names.index(reference_method)

    methods = [[] for _ in range(len(names))]

    for combination in param_nums:
        methods[0].append(combination[0]/combination[reference_location])
        methods[1].append(combination[1]/combination[reference_location])
        methods[2].append(combination[2]/combination[reference_location])
        methods[3].append(combination[3]/combination[reference_location])
        methods[4].append(combination[4]/combination[reference_location])
        methods[5].append(combination[5]/combination[reference_location])

    tb_dict = dict()
    for m in range(len(methods)):
        tb_dict[names[m]] = np.mean(methods[m])

    # tb_dict looks something like:
    # {'ctr_all': 0.021869912051460418,
    # 'ctr_m': 0.52751028000383549,
    # 'fg': 1.0,
    # 'im': 0.71620368516520383,
    # 'ml': 0.86533500528478879,
    # 'wt': 0.67474526297791027}

    # save dictionary as dataframe and output as csv
    pd.DataFrame(tb_dict.items(),
                 columns=['Method', 'Ratio']).to_csv('./results/tb_dict.csv')

    # plots bar graph of ratio of TP/FP of each method to chosen reference
    # method
    plt.bar(range(len(tb_dict)), tb_dict.values(), align='center', color='w')
    plt.xticks(range(len(tb_dict)), tb_dict.keys())
    plt.ylabel('Ratio of TP/FP to reference method')
    plt.xlabel('Method')
    plt.show()
