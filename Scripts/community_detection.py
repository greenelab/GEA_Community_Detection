# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 10:10:37 2016

@author: LiaHarrington

Usage:
Imported by Enrichment_Testing.py and Enrichment_Main.py

Description:
Creates network from IMP and subgraph using m selected pathways

"""

# ---------------------------------------------------------------------------
# Import packages and set seed
# ---------------------------------------------------------------------------
import random
import igraph
import os

random.seed(123)
# ---------------------------------------------------------------------------
# Reads in IMP ege-list and creates complete IMP network
# ---------------------------------------------------------------------------

# read in IMP network

imp_file = os.path.join('Data', 'global_average.filtered.txt')
IMP_NETWORK = igraph.Graph.Read_Ncol(imp_file, directed=False)

# ---------------------------------------------------------------------------
# Creates subgraph with genes from the m pathways
# ---------------------------------------------------------------------------

def weights():
    '''
    Description
    :Creates the weights for the IMP network construction using the bayesian
    posterior weights from the IMP network. Weights are pulled from the
    'weights_full.txt' or can be gotten from the third column of the
    IMP dataset.

    Arguments
    :none

    Output
    :returns a list of relevant IMP edgelist weights
    '''

    m_weights = open('weights_full.txt')

    imp_weights = []
    for weight in m_weights:
        weight = float(weight.rstrip().split()[0])
        imp_weights.append(weight)
    return imp_weights

def community_detection(sub_genes, method, weights):
    '''
    Description
    :Performs community detection using the genes from the m
    pathways to create the relevant subgraphs with only those gene nodes

    Arguments
    :param sub_genes: the gene list created using m_genes_list
    :param method: 'fastgreedy', 'walktrap', 'infomap' or 'multilevel'
    :param weights: either IMP weights using 'weights' or no weights using
    'NULL'

    Output
    :returns the community detection object which contains all the communities
    and corresponding genes as a list of list.

    '''

    sub_s = IMP_NETWORK.subgraph(sub_genes)

    if weights != 'NULL':
        if method == 'fastgreedy':
            sub_com = sub_s.community_fastgreedy(weights).as_clustering()
        if method == 'walktrap':
            sub_com = sub_s.community_walktrap(weights).as_clustering()
        if method == 'infomap':
            sub_com = sub_s.community_infomap(weights)
        if method == 'multilevel':
            sub_com = sub_s.community_multilevel(weights)
    elif weights == 'NULL':
        if method == 'fastgreedy':
            sub_com = sub_s.community_fastgreedy().as_clustering()
        if method == 'walktrap':
            sub_com = sub_s.community_walktrap().as_clustering()
        if method == 'infomap':
            sub_com = sub_s.community_infomap()
        if method == 'multilevel':
            sub_com = sub_s.community_multilevel()

    return sub_com

def index_to_edge_name(sub_com):
    '''
    Description
    :Transforms the index labeled nodes in each community from sub_com to thier
    actual names

    Arguments
    :param sub_com: the object returned by community_detection()

    Output
    :a list of list of name labled nodes by community instead of index labled
    nodes

    '''

    subcom_named_groups = []
    for group in range(len(sub_com)):
        com_group = [sub_com.graph.vs['name'][sub_com[group][i]] for i in
                     range(len(sub_com[group]))]
        subcom_named_groups.append(com_group)
    return subcom_named_groups
