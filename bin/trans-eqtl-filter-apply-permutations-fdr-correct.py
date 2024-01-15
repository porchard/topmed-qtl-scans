#!/usr/bin/env python3.8
# coding: utf-8


import pandas as pd
import numpy as np
import sys
import argparse
import pickle
from tensorqtl import trans

parser = argparse.ArgumentParser()
parser.add_argument('--pairs', required=True)
parser.add_argument('--permutations', required=True)
# parser.add_argument('--biotypes', default=['protein_coding', 'lincRNA'], nargs='+')
args = parser.parse_args()

#PAIRS_DF = 'Monocyte.trans_qtl_pairs.with_mappability.txt.gz'
#PERMUTATIONS_PICKLE = 'Monocyte.permutations.pickle'
PAIRS_DF = args.pairs
PERMUTATIONS_PICKLE = args.permutations
BIOTYPES = ['protein_coding', 'lincRNA']

def padjust_bh(p):
    """
    Benjamini-Hochberg adjusted p-values
    Replicates p.adjust(p, method="BH") from R
    From pyqtl: https://github.com/broadinstitute/pyqtl/blob/master/qtl/stats.py
    """
    n = len(p)
    i = np.arange(n,0,-1)
    o = np.argsort(p)[::-1]
    ro = np.argsort(o)
    return np.minimum(1, np.minimum.accumulate(float(n)/i * np.array(p)[o]))[ro]


with open(PERMUTATIONS_PICKLE, 'rb') as f:
    permutations_df = pickle.load(f)
    permutations_df.index = ('chr' + permutations_df.index).to_list()

pairs_df = pd.read_csv(PAIRS_DF, sep='\t')


# mappability and biotype filters
pairs_df = pairs_df[(pairs_df.gene_mappability>=0.8) & (~pairs_df.gene_crossmaps_to_gene_near_variant) & (pairs_df.biotype.isin(BIOTYPES))]

# filter to inter-chromosomal only
pairs_df['variant_chr'] = pairs_df.variant_id.str.split('_', expand=True)[0]
pairs_df = pairs_df[pairs_df.variant_chr != pairs_df.phenotype_chr]

# keep top p-value per gene
top_df = pairs_df.groupby('phenotype_id').apply(lambda x: x.loc[x['pval'].idxmin()])
assert(pairs_df.phenotype_id.nunique() == len(top_df))

# apply permutations
trans.apply_permutations(permutations_df, top_df)

# fdr correct
top_df['pval_beta_no_zero'] = np.maximum(top_df.pval_beta, 1e-300)
top_df['qvalue'] = padjust_bh(top_df['pval_beta'])

# output
top_df.to_csv(sys.stdout, sep='\t', index=False)
