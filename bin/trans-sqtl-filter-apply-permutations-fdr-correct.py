#!/usr/bin/env python3.8
# coding: utf-8


import pandas as pd
import numpy as np
import sys
import argparse
import pickle
from tensorqtl import trans
import scipy.stats as stats
import re


parser = argparse.ArgumentParser()
parser.add_argument('--pairs', required=True)
parser.add_argument('--permutations', required=True)
parser.add_argument('--phenotype-groups', required=True)
# parser.add_argument('--biotypes', default=['protein_coding', 'lincRNA'], nargs='+')
args = parser.parse_args()

#PAIRS_DF = 'Monocyte.trans_qtl_pairs.with_mappability.txt.gz'
#PERMUTATIONS_PICKLE = 'Monocyte.permutations.pickle'
PAIRS_DF = args.pairs
PERMUTATIONS_PICKLE = args.permutations
BIOTYPES = ['protein_coding', 'lincRNA']


def phenotype_id_to_gene_id(x):
    # first, try to match with a version
    ENSEMBL_RE_WITH_VERSION = 'ENSG\d+\.\d+'
    ENSEMBL_RE_WITHOUT_VERSION = 'ENSG\d+'
    with_version = re.search(ENSEMBL_RE_WITH_VERSION, x)
    without_version = re.search(ENSEMBL_RE_WITHOUT_VERSION, x)
    if with_version:
        return with_version.group(0)
    elif without_version:
        return without_version.group(0)
    else:
        raise ValueError(f'Not able to infer gene ID from {x}')


def correct_across_phenotypes(beta_shape1, beta_shape2, pval_true_df, number_phenotypes):
    return 1 - ((1 - stats.beta.cdf(pval_true_df, beta_shape1, beta_shape2)) ** number_phenotypes)


def correct_across_phenotypes_log_space(beta_shape1, beta_shape2, pval_true_df, number_phenotypes):
    # same as the above function, but with additional precision
    x = np.log1p(-1*stats.beta.cdf(pval_true_df, beta_shape1, beta_shape2))
    return -1*np.expm1(x * number_phenotypes)



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

pairs_df['gene_id'] = pairs_df.phenotype_id.map(phenotype_id_to_gene_id)
# get phenotypes per gene
phenotype_groups = pd.read_csv(args.phenotype_groups, sep='\t', header=None, names=['phenotype_id', 'gene_id'])
phenotypes_per_gene = phenotype_groups.groupby('gene_id').phenotype_id.nunique().to_dict()

# keep top p-value per gene
top_df = pairs_df.groupby('gene_id').apply(lambda x: x.loc[x['pval'].idxmin()])
assert(pairs_df.gene_id.nunique() == len(top_df))

# apply permutations
trans.apply_permutations(permutations_df, top_df)
top_df['pval_beta_no_zero'] = np.maximum(top_df.pval_beta, 1e-300)

# correct across phenotypes
top_df['phenotypes_tested_for_gene'] = top_df.gene_id.map(phenotypes_per_gene)
top_df['pval_beta_corrected_across_phenotypes'] = [correct_across_phenotypes_log_space(beta_shape1, beta_shape2, pval_true_df, phenotypes_per_gene[gene_id]) for beta_shape1, beta_shape2, pval_true_df, gene_id in zip(top_df.beta_shape1, top_df.beta_shape2, top_df.pval_true_df, top_df.gene_id)]
top_df['pval_beta_corrected_across_phenotypes_no_zero'] = np.maximum(top_df.pval_beta_corrected_across_phenotypes, 1e-300)

# if have only one phenotype, pval_beta should be same as pval_beta_corrected_across_phenotypes
# approx, to allow for precision issues
assert((np.log10(top_df[top_df.phenotypes_tested_for_gene==1].pval_beta) - np.log10(top_df[top_df.phenotypes_tested_for_gene==1].pval_beta_corrected_across_phenotypes)).abs().max() < 1e-10)


# fdr correct
top_df['qvalue'] = padjust_bh(top_df['pval_beta_corrected_across_phenotypes'])

# output
top_df.to_csv(sys.stdout, sep='\t', index=False)