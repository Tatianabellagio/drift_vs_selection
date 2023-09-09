# once i have wholegenome_offset.trees
import pandas as pd
import tskit
import allel
import random
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import tsinfer
import pyslim
import json
import os
from collections import defaultdict

def overlap_neutral_mut (ts_new, ts, mapper_realid_metadataid):
    ## extract surviving ndoes and comapre them to our old ndoes to place mtuations in the right place
    surviving_nodes = []
    for i in ts_new.tables.nodes:
        surviving_nodes.append(i.metadata['slim_id'])
    ## new nodes id and the ids i gave them in the past
    new_mapper = pd.DataFrame({'new_ids': range(0, len(ts_new.tables.nodes)), 'my_ids_metadata':surviving_nodes})
    ## map old nodes with new nodes
    mapper_lost_nodes = new_mapper.merge(mapper_realid_metadataid, left_on = 'my_ids_metadata', right_on = 'my_ids_metadata', how= 'right')

    ## create a mask to only keep from the old nodes the ones that survived the simulation
    mask = mapper_lost_nodes['new_ids'].notna()

    tables_og = ts.dump_tables()

    ## now filter old tables only based on surviving nodes 
    tables_og.nodes.replace_with(tables_og.nodes[mask])

    ## now filter mutation table based on the surviving nodes, for that, extract the nodes 
    old_nodes = tables_og.mutations.node

    old_nodes = pd.Series(old_nodes)

    old_nodes.name = 'old_nodes'

    ## create a dataframe relating the new and old nodes
    replace_oldbynew_nodes = pd.merge(old_nodes, mapper_lost_nodes, left_on ='old_nodes', right_on = 'real_id', how= 'left')

    ## create a mask to filter out all the mutations than has been lost 
    mask_mutations_lost = replace_oldbynew_nodes['new_ids'].notna()

    ## filter out mutations that has been lost 
    table_mutations = tables_og.mutations[mask_mutations_lost]

    ## now replace the old nodes ids by the new nodes ids with the mapper
    ids_to_replace = replace_oldbynew_nodes.dropna()['new_ids']
    table_mutations.node = np.array(ids_to_replace.astype('int32'))

    ## and jsut set the sites from 0 to the length of mutation table 
    table_mutations.site = np.array(range(0, len(table_mutations))).astype('int32')

    ## apply the same filter from the mutations table to the sites table 
    table_sites = tables_og.sites[mask_mutations_lost]  

    ## now replace all this filter old tables in the new tree seq! 
    new_tables = ts_new.dump_tables()

    new_tables.mutations.replace_with(table_mutations)

    new_tables.sites.replace_with(table_sites)

    ## make sure to compure mutations parents
    new_tables.compute_mutation_parents()

    ## create tree seq based on tables
    tree_nm = new_tables.tree_sequence()

    return tree_nm.simplify()

def convert_tree_to_vcf (tree,name_vcf):
    # create a vcf file from the treeseq 
    with open(name_vcf, 'w') as file:
        # Pass the file object as the output parameter
        tree.write_vcf(output=file)
        
#import the old tree
#ts_old = tskit.load("../treeseq/wholegenome_offset_baselinetree.trees")
#import mapper old nodes to new nodes
#mapper_realid_metadataid = pd.read_csv('../treeseq/mapper_realid_metadataid_wholegenome.csv')

def generate_nonhet_pos():
    vcf = allel.read_vcf('../treeseq/wholegenome_offset.vcf')
    ## extract the genotype from the vcf file
    geno = vcf["calldata/GT"]
    ## calcualte the genotype counts
    geno = allel.GenotypeArray(vcf["calldata/GT"])
    het_sites = geno.count_het(axis=1)
    mask_non_het_sites = het_sites==0 
    nonhet_pos = vcf['variants/POS'][mask_non_het_sites]
    pd.Series(nonhet_pos).to_csv('nonhet_pos_wholegenome_offset.csv')
    return

def filtering_pos (nonhet_pos, pos_new, geno_og, geno_new):
    pos_to_keep = np.intersect1d(nonhet_pos, pos_new)
    mask_pos_ogvcf = pd.Series(pos_og).isin(pos_to_keep)
    geno_og_rpos  = geno_og[mask_pos_ogvcf]
    mask_pos_newvcf = pd.Series(pos_new).isin(pos_to_keep)
    geno_new_rpos  = geno_new[mask_pos_newvcf]
    return geno_og_rpos, geno_new_rpos
    

def get_ecotype_geno_mapper(geno_og_rpos):
    geno_og_rpos = np.swapaxes(geno_og_rpos, 0, 1)
    ecotype_geno_mapper = {}
    for i,j in zip(geno_og_rpos, samples):
        geno = i.tobytes()
        ecotype_geno_mapper[geno] = j
    return ecotype_geno_mapper

def get_ecotype_counts(geno_new_rpos, pop_name):
    geno_new_rpos = np.swapaxes(geno_new_rpos, 0, 1)
    # Initialize a defaultdict to store the genotype counts
    ecotype_counts = defaultdict(int)
    # Count genotypes in geno_drift_resh
    for i in geno_new_rpos:
        sample = i.tobytes()
        ecotype = ecotype_geno_mapper.get(sample, 'other')
        ecotype_counts[ecotype] += 1
    name = 'count'+ pop_name
    ecotype_countsdf = pd.DataFrame(list(ecotype_counts.items()), columns=['ecotype', name])
    ecotype_countsdf['ecotype'] = ecotype_countsdf['ecotype'].str.split('_').str[0]
    return ecotype_countsdf

nonhet_pos = np.array(pd.read_csv('nonhet_pos_wholegenome_offset.csv',usecols=[1])['0'])

vcf_og = allel.read_vcf('../treeseq/wholegenome_offset.vcf', fields=['calldata/GT', 'variants/POS', 'samples'])
geno_og = vcf_og['calldata/GT']
pos_og = vcf_og['variants/POS']
samples = vcf_og['samples']


variances = pd.read_csv('variances.txt',header=None)[0]
optimas = pd.read_csv('optimas.txt',header=None)[0]
print('?')
## import each of the selection trees, add mutations and save their vcf file 
#for i in variances:
#    for j in optimas:
#        file_name = f'output_selection/var{i}_optima{j}_result.trees'
#        if os.path.exists(file_name) and os.path.getsize(file_name) < 0:
#            print('empty_tree')
#            with open(f'output_selection/result_selection_var{i}_optima{j}.vcf', "w"):
#                pass  # Create an empty vcf file 
#        elif os.path.exists(file_name) and os.path.getsize(file_name) > 0:
#            print('tree with individuals')
#            ts_new = tskit.load(file_name)
#            ts_nm = overlap_neutral_mut(ts_new, ts_old, mapper_realid_metadataid)
#            convert_tree_to_vcf(ts_nm, f'output_selection/result_selection_var{i}_optima{j}.vcf')

ecotypes_grenenet = pd.read_csv('ecotypes_grenenet.txt',header=None, dtype=object)
ecotypes_grenenet.columns= ['ecotype']
ecotypes_grenenet = pd.concat([ecotypes_grenenet, pd.DataFrame(data = {'ecotype': ['other']}, index=[231])],axis=0)
print('hola')

for i in variances:
    print(i)
    for j in optimas:
        vcf_filename = f'output_selection/result_selection_var{i}_optima{j}.vcf'
        print(vcf_filename)
        # Check if the VCF file exists and is not empty (all individuals died)
        if os.path.exists(vcf_filename) and os.path.getsize(vcf_filename) > 0:
        
            ## import each of teh vcf realted to the drift simulation
            vcf_new = allel.read_vcf(vcf_filename, fields = ['calldata/GT', 'variants/POS'])
            ##extract the posisions and the geno array
            pos_new = vcf_new['variants/POS']
            geno_new = vcf_new['calldata/GT']
            ## for each of them create the ecotype geno mapper, depending on the positions that made it 
            geno_og_rpos, geno_new_rpos = filtering_pos (nonhet_pos, pos_new, geno_og, geno_new)
            ecotype_geno_mapper = get_ecotype_geno_mapper(geno_og_rpos)
            ecotype_countsdf = get_ecotype_counts(geno_new_rpos, f'_selection_var{i}_optima{j}')
            print(ecotype_countsdf)
            ## merge with previous 
            ecotypes_grenenet = ecotypes_grenenet.merge(ecotype_countsdf, how='left', on ='ecotype')
            print(ecotypes_grenenet)
        elif os.path.exists(vcf_filename) and os.path.getsize(vcf_filename) == 0:
            empty_col_name = f'counts_selection_var{i}_optima{j}'
            print('added empty column in ' + empty_col_name)
            ecotypes_grenenet[empty_col_name] = 0
            
ecotypes_grenenet.to_csv('ecotype_counts_selection.csv')