#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--prefix', default='pcs.', help='Prefix for output files')
parser.add_argument('--fdr', type=float, default=0.05, help='FDR threshold. default = 0.05')
parser.add_argument('--permutation-files', nargs='+', dest='permutation_files', help='Permutation files. Filenames must match {tissue}.{phenotype_pcs}.cis_qtl.txt.gz; files must contain columns qval, pval_perm, pval_beta, tss_distance')
args = parser.parse_args()

PREFIX = args.prefix
PERMUTATION_FILES = args.permutation_files

FDR_COLUMN = 'qval'

def load_permutation_file(f):
    tissue, pcs = os.path.basename(f).split('.')[:2]
    tmp = pd.read_csv(f, sep='\t')
    tmp['tissue'] = tissue
    tmp['pcs'] = int(pcs)
    return tmp


results = pd.concat([load_permutation_file(f) for f in PERMUTATION_FILES])

# check that fitted p-values are well-calibrated
g = sns.FacetGrid(results, col='pcs', col_wrap=4)
g.map(plt.scatter, 'pval_perm', 'pval_beta')
for ax in g.axes.flatten():
    ax.plot([0, 1], [0, 1], c='red', linestyle='--')
g.fig.savefig(f'{PREFIX}fitted-vs-empirical-p.png')
g.fig.clf()


# plot distance between gene and top variants
g = sns.FacetGrid(results[results[FDR_COLUMN]<=args.fdr].assign(dist_trans=lambda df: df.tss_distance/1e3), col='pcs', col_wrap=4)
g.map(plt.hist, 'dist_trans', bins=30)
g.set_xlabels('Distance from eGene to top variant (kb)')
g.fig.tight_layout()
g.fig.savefig(f'{PREFIX}distance-from-gene-to-top-variant.png')
g.fig.clf()


# plot number of eGenes as a function of PCs
number_pcs_vs_egenes = results.groupby('pcs')[FDR_COLUMN].apply(lambda x: sum(x<=args.fdr)).rename('egenes').reset_index()
if results.groupby('pcs').phenotype_id.nunique().nunique() != 1:
    raise ValueError('Different numbers of genes were tested with different numbers of PCs')
GENES_TESTED = results.groupby('pcs').phenotype_id.nunique().unique()[0]

fig, ax = plt.subplots()
sns.scatterplot(x='pcs', y='egenes', data=number_pcs_vs_egenes, ax=ax)
sns.lineplot('pcs', 'egenes', data=number_pcs_vs_egenes, ax=ax)
ax.axhline(GENES_TESTED, color='black', linestyle='--')
ax.set_xlabel('PCs')
ax.set_ylabel('Number eGenes')
ax.grid(True)
fig.tight_layout()
fig.savefig(f'{PREFIX}number-pcs-vs-number-egenes.png')
fig.clf()