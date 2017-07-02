#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 00:23:54 2017

@author: LiaHarrington
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 00:21:52 2017

@author: LiaHarrington
"""

from community_detection import community_detection, index_to_edge_name
from enrichment_testing import enrichment
import pickle
import pandas as pd
import os

all_genes_file = os.path.join('Data', 'PID_all_genes.pkl')
ALL_GENES = pickle.load(open(all_genes_file, 'r'))

path_genes_file = os.join.path('Data', 'PID_path_genes.pkl')
PATH_GENES = pickle.load(open(path_genes_file, 'r'))

path_names_file = os.join.path('Data', 'PID_path_names.pkl')
PATH_NAMES = pickle.load(open(path_names_file, 'r'))

imp_file = os.join.path('Data', 'IMP_genes.pkl')
IMP_GENES = pickle.load(open(imp_file, 'r'))

hgsc_file = os.join.path('Data', 'entrezid_hgsc.txt')
hgsc = pd.read_csv(hgsc_file, sep='\t')

hgsc.reset_index(level=0, inplace=True)

hgsc.rename(columns={'index': 'Genes'}, inplace=True)

# create genelist for each cluster from hgsc data
k4c1 = hgsc.Genes[hgsc.K4Cluster1_Pos == 1].append(hgsc.Genes[hgsc.K4Cluster1_Neg == 1])
k4c2 = hgsc.Genes[hgsc.K4Cluster2_Pos == 1].append(hgsc.Genes[hgsc.K4Cluster2_Neg == 1])
k4c3 = hgsc.Genes[hgsc.K4Cluster3_Pos == 1].append(hgsc.Genes[hgsc.K4Cluster3_Neg == 1])
k4c4 = hgsc.Genes[hgsc.K4Cluster4_Pos == 1].append(hgsc.Genes[hgsc.K4Cluster4_Neg == 1])

k3c1 = hgsc.Genes[hgsc.K3Cluster1_Pos == 1].append(hgsc.Genes[hgsc.K3Cluster1_Neg == 1])
k3c2 = hgsc.Genes[hgsc.K3Cluster2_Pos == 1].append(hgsc.Genes[hgsc.K3Cluster2_Neg == 1])
k3c3 = hgsc.Genes[hgsc.K3Cluster3_Pos == 1].append(hgsc.Genes[hgsc.K3Cluster3_Neg == 1])

k2c1 = hgsc.Genes[hgsc.K2Cluster1_Pos == 1].append(hgsc.Genes[hgsc.K2Cluster1_Neg == 1])
k2c2 = hgsc.Genes[hgsc.K2Cluster2_Pos == 1].append(hgsc.Genes[hgsc.K2Cluster2_Neg == 1])

# create combined genelist over the cluster for each k
all_genes_lst_k4 = [k4c1, k4c2, k4c3, k4c4]
all_genes_lst_k3 = [k3c1, k3c2, k3c3]
all_genes_lst_k2 = [k2c1, k2c2]

def not_in(gene_list):
    '''
    Description
    Finds genes either not in IMP or PID ontology to help
    identify which genes are missing from which source

    Arguments
    :gene_lst: the gene_lst of interest

    Output
    :not_imp: genes in gene_lst not in IMP
    :not_all: genes in gene_lst not in ontology (ALL_GENES)
    :both: genes neither in IMP or ontology
    '''

    not_imp = []
    not_all = []
    both = []

    for gene in gene_list:
        if gene not in IMP_GENES and gene in ALL_GENES:
            not_imp.append(gene)
        elif gene not in ALL_GENES and gene in IMP_GENES:
            not_all.append(gene)
        else:
            both.append(gene)
    return not_imp, not_all, both

def convert_string(gene_list):
    '''
    Description
    converts the list of integer gene ids to string type

    Arguments
    :gene_list: the list of genes of interest

    Output
    :A list of genes in string form instead of integer. This is for compatibility
    with PATH_GENES and ALL_GENES
    '''

    return [str(g) for g in gene_list]

def remove_not_in_IMP(gene_list):
    '''
    Description
    Removes all genes in the gene_list that are not in IMP

    Arguments
    :gene_list: list of genes of interest

    Output
    :Returns list of genes that have genes not in IMP removed. This is to
    help ensure that the sub_graph from the gene list can be created from IMP
    network for community detection methods.

    '''
    missing = []
    for gene in gene_list:
        if gene not in IMP_GENES:
            missing.append(gene)
    return [gene for gene in gene_list if gene not in missing]

def cd_gea_pathways(all_genes_lst, all_names_lst, com_method, alpha=.05,
                    min_com_size=3, weights=None):
    '''
    Description
    Loops over all the gene lists using selected community detection method
    and returns a dataframe with significant pathway_name, pvalue, cluster id
    and community number

    Arguments
    :all_gene_lst: list of list of genes in each cluster for specified hgsc k
    :all_names_lst: list of cluster name, for example ['k2c1', 'k2c2']
    :com_method: can be 'fastgreedy', 'walktrap', 'infomap', or 'multilevel'
    :alpha: preset to .05
    :min_com_size: minimum size of community for it to be considered for further analysis
    :weights: weights used in IMP network, default is None

    Output
    :dataframe with pathway_name, pval, cluster id and community id for each significant
    pathway
    '''

    path_dict = dict.fromkeys(all_names_lst, [])

    keys = sorted(path_dict.keys())

    cluster_df = pd.DataFrame()
    for i, lst in enumerate(all_genes_lst):
        gene_list = remove_not_in_IMP(convert_string(lst))

        cd_genes = community_detection(gene_list, com_method, weights=None)
        cd_genes_lst = index_to_edge_name(cd_genes)

        cd_com_lst = []
        for com in cd_genes_lst:
        # keep communities with at min community size
            if len(com) >= min_com_size:
                cd_com_lst.append(com)

        for c, com in enumerate(cd_com_lst):
            cd_com_dict = dict()
            results = enrichment(com, PATH_GENES, alpha, ALL_GENES)
            # check if top path signif
            if results[0][0][1]:
                cd_com_dict[PATH_NAMES[results[0][0][0]]] = [results[2][results[0][0][0]]]


            community_df = pd.DataFrame(cd_com_dict)
            community_melt_df = pd.melt(community_df)
            community_melt_df = community_melt_df.assign(cluster='{0}'.format(keys[i]))
            community_melt_df = community_melt_df.assign(community=c)

            # check non-empty dataframe
            if len(community_melt_df) > 0:

                hgsc_name_dict = {'k2c1': ['K2Cluster1_Pos','K2Cluster1_Neg'],
                                  'k2c2': ['K2Cluster2_Pos','K2Cluster2_Neg'],
                                  'k3c1': ['K3Cluster1_Pos','K3Cluster1_Neg'],
                                  'k3c2': ['K3Cluster2_Pos','K3Cluster2_Neg'],
                                  'k3c3': ['K3Cluster3_Pos','K3Cluster3_Neg'],
                                  'k4c1': ['K4Cluster1_Pos','K4Cluster1_Neg'],
                                  'k4c2': ['K4Cluster2_Pos','K4Cluster2_Neg'],
                                  'k4c3': ['K4Cluster3_Pos','K4Cluster3_Neg'],
                                  'k4c4': ['K4Cluster4_Pos','K4Cluster4_Neg']}

                # must do zero-index to actually just get the cluster string object
                cluster_pos, cluster_neg = hgsc_name_dict[community_melt_df.cluster[0]]

                cluster_pos_genes = convert_string(hgsc.Genes[hgsc[cluster_pos] == 1])
                cluster_neg_genes = convert_string(hgsc.Genes[hgsc[cluster_neg] == 1])

                relv_path_id = PATH_NAMES.index(community_melt_df.variable[0])
                relv_path_genes = PATH_GENES[relv_path_id]

                path_pos = set(cluster_pos_genes).intersection(relv_path_genes)
                path_neg = set(cluster_neg_genes).intersection(relv_path_genes)

                path_com_intersect = set(relv_path_genes).intersection(set(com))

                path_com_intersect_pos = set(cluster_pos_genes).intersection(
                    path_com_intersect)
                path_com_intersect_neg = set(cluster_neg_genes).intersection(
                    path_com_intersect)

                community_melt_df = community_melt_df.assign(
                    num_pos_genes_pathway=len(path_pos))
                community_melt_df = community_melt_df.assign(
                    num_neg_genes_pathway=len(path_neg))
                community_melt_df = community_melt_df.assign(
                    pathway_com_overlap_pos=str(list(path_com_intersect_pos))).astype(object)
                community_melt_df = community_melt_df.assign(
                    pathway_com_overlap_neg=str(list(path_com_intersect_neg))).astype(object)

            cluster_df = cluster_df.append(community_melt_df, ignore_index=True)

        cluster_df2 = cluster_df.rename(columns = {'variable': 'pathway_name', 'value': 'pval'})

    return cluster_df2

master_genes_lst = [all_genes_lst_k4, all_genes_lst_k3, all_genes_lst_k2]
master_namelst = [['k4c1', 'k4c2', 'k4c3', 'k4c4'], ['k3c1', 'k3c2', 'k3c3'], ['k2c1', 'k2c2']]

def run_all_cd(master_genes_lst, master_namelst):
    '''
    Description
    Creates master dataframe over all community detection methods and k=2,3,4

    Arguments
    :master_genes_lst: list of the list of list of hgsc genes for specified k
    :master_name_lst: list of list of names for each k and respective cluster

    Output
    :dataframe of pathway_name, pval, cluster, community for all community detection
    methods and k values

    '''

    methods = ['fastgreedy', 'walktrap', 'infomap', 'multilevel']
    cols = ['pathway_name', 'pval', 'cluster', 'community', 'method', 'num_pos_genes_pathway',
            'num_neg_genes_pathway', 'pathway_com_overlap_pos', 'pathway_com_overlap_neg']

    master_df = pd.DataFrame(columns=cols)

    for idx, gl in enumerate(master_genes_lst):
        for m in methods:
            cd_df = cd_gea_pathways(gl,master_namelst[idx], m, alpha=.05, 
                                    min_com_size=3, weights=None)
            if cd_df.shape[0]>0:
                cd_df = cd_df.assign(method = m)
                master_df = master_df.append(cd_df, ignore_index=True)

    master_df.to_csv('./Data/master_cd_runs.csv', sep =',', index=False)

run_all_cd(master_genes_lst, master_namelst)
    
if __name__ == '__main__': 
    run_all_cd(master_genes_lst = master_genes_lst, master_namelst = master_namelst)
