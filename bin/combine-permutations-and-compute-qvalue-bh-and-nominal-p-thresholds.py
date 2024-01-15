#!/usr/bin/env python

import sys
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import logging

PREFIX = sys.argv[1]
PERMUTATION_FILES = sys.argv[2:]

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')


def calculate_qvalues(res_df, fdr=0.05):
    """Annotate permutation results with q-values, p-value threshold"""

    logging.info('Computing q-values')
    logging.info(f'  * Number of phenotypes tested: {res_df.shape[0]}')
    r = stats.pearsonr(res_df['pval_perm'], res_df['pval_beta'])[0]
    logging.info(f'  * Correlation between Beta-approximated and empirical p-values: {r:.4f}')
    
    qval = multipletests(res_df.pval_beta, method='fdr_bh')[1]
    res_df['qval'] = qval
    logging.info(f"  * QTL phenotypes @ FDR {fdr:.2f}: {(res_df['qval'] <= fdr).sum()}")

    # determine global min(p) significance threshold and calculate nominal p-value threshold for each gene
    lb = res_df.loc[res_df['qval'] <= fdr, 'pval_beta'].sort_values()
    ub = res_df.loc[res_df['qval'] > fdr, 'pval_beta'].sort_values()

    if lb.shape[0] > 0:  # significant phenotypes
        lb = lb[-1]
        if ub.shape[0] > 0:
            ub = ub[0]
            pthreshold = (lb+ub)/2
        else:
            pthreshold = lb
        logging.info(f'  * min p-value threshold @ FDR {fdr}: {pthreshold:.6g}')
        res_df['pval_nominal_threshold'] = stats.beta.ppf(pthreshold, res_df['beta_shape1'], res_df['beta_shape2'])

permutations = pd.concat([pd.read_csv(f, sep='\t') for f in PERMUTATION_FILES]).set_index('phenotype_id')
for i in ['qval', 'pval_nominal_threshold']:
    if i in permutations.columns.to_list():
        permutations = permutations.drop(columns=[i])
calculate_qvalues(permutations)
out_file = f'{PREFIX}.cis_qtl.txt.gz'
permutations.to_csv(out_file, sep='\t', float_format='%.6g')
