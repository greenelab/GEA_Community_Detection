#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 00:21:52 2017

@author: LiaHarrington
"""

from community_detection import community_detection, index_to_edge_name
from enrichment_testing import enrichment

ALL_GENES = pickle.load(open('pid_all_genes.pkl', 'r'))
PATH_GENES = pickle.load(open('pid_path_genes.pkl', 'r'))
PATH_NAMES = pickle.load(open('pid_path_names.pkl', 'r'))


hgsc = pd.read_csv('entrezid_hgsc.txt', sep = '\t')
hgsc.reset_index(level=0, inplace=True)

hgsc.rename(columns={'index': 'Genes'}, inplace=True)

k4c1 = hgsc.Genes[hgsc.K4Cluster1_Pos==1].append(hgsc.Genes[hgsc.K4Cluster1_Neg==1])
k4c2 = hgsc.Genes[hgsc.K4Cluster2_Pos==1].append(hgsc.Genes[hgsc.K4Cluster2_Neg==1])
k4c3 = hgsc.Genes[hgsc.K4Cluster3_Pos==1].append(hgsc.Genes[hgsc.K4Cluster3_Neg==1])
k4c4 = hgsc.Genes[hgsc.K4Cluster4_Pos==1].append(hgsc.Genes[hgsc.K4Cluster4_Neg==1])

k3c1 = hgsc.Genes[hgsc.K3Cluster1_Pos==1].append(hgsc.Genes[hgsc.K3Cluster1_Neg==1])
k3c2 = hgsc.Genes[hgsc.K3Cluster2_Pos==1].append(hgsc.Genes[hgsc.K3Cluster2_Neg==1])
k3c3 = hgsc.Genes[hgsc.K3Cluster3_Pos==1].append(hgsc.Genes[hgsc.K3Cluster3_Neg==1])

k2c1 = hgsc.Genes[hgsc.K2Cluster1_Pos==1].append(hgsc.Genes[hgsc.K2Cluster1_Neg==1])
k2c2 = hgsc.Genes[hgsc.K2Cluster2_Pos==1].append(hgsc.Genes[hgsc.K2Cluster2_Neg==1])

all_genes_lst_k4 = [k4c1, k4c2, k4c3, k4c4]
all_genes_lst_k3 = [k3c1, k3c2, k3c3]
all_genes_lst_k2 = [k2c1, k2c2]


def not_in(gene_list):
    '''finds genes either not in IMP or PID ontology'''
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
    '''converts the list of integer gene ids to string type'''
    return [str(g) for g in gene_list]

def remove_not_in_IMP(gene_list):
    '''removes all genes in the gene_list that are not in IMP'''
    missing = []
    for gene in gene_list: 
        if gene not in IMP_GENES: 
            missing.append(gene)
    return [gene for gene in gene_list if gene not in missing]

all_genes_lst = [k4c1, k4c2, k4c3, k4c4]
all_names_lst = ['k4c1', 'k4c2', 'k4c3', 'k4c4']

def cd_gea_pathways(all_genes_lst, all_names_lst, com_method, alpha=.05, min_com_size=3, weights=None): 
    '''loops over all the gene lists using selected community detection method 
        and returns a dictionary with cluster id as key and the community genes
        along with identified significant pathways'''
        
    #all_genes_lst = all_genes_lst_k3
    #all_names_lst = ['k3c1', 'k3c2', 'k3c3']
#    com_method = 'walktrap'
#    min_com_size = 3
#    weights = None
#    alpha = .05
#    
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
                
         
            com_df = pd.DataFrame(cd_com_dict)
            cd_com2 = pd.melt(com_df)
            cd_com2 = cd_com2.assign(cluster = '{0}'.format(keys[i]))
            cd_com2 = cd_com2.assign(community = c)
            
            cluster_df = cluster_df.append(cd_com2, ignore_index=True)
    
    cluster_df.columns = ['pathway_name', 'pval', 'cluster', 'community' ]
    
    return cluster_df
            
master_genes_lst = [all_genes_lst_k4, all_genes_lst_k3, all_genes_lst_k2]
master_namelst = [['k4c1', 'k4c2', 'k4c3', 'k4c4'], ['k3c1', 'k3c2', 'k3c3'], ['k2c1', 'k2c2']]
    
def run_all_cd(master_genes_lst, master_namelst): 
    methods = ['fastgreedy', 'walktrap', 'infomap', 'multilevel']
    cols = ['pathway_name', 'pval', 'cluster', 'community', 'method']
    
    master_df = pd.DataFrame(columns=cols)
    
    for idx, gl in enumerate(master_genes_lst): 
        for m in methods: 
            cd_df = cd_gea_pathways(gl,master_namelst[idx], m, alpha=.05, min_com_size=3, weights=None)
            if cd_df.shape[0]>0: 
                cd_df = cd_df.assign(method = m) 
                master_df = master_df.append(cd_df, ignore_index=True)
        
    #return master_df
    master_df.to_csv('./Data/master_cd_runs.csv', sep =',', index=False)
        