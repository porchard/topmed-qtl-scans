#!/usr/bin/env python

import numpy as np
import pandas as pd
import scipy.stats as stats
import sys
import os
import glob
from datetime import datetime
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--permutations', required=True, help='Path to single permutation file')
parser.add_argument('--groups', required=True, help='Path to permutation groups file')
parser.add_argument('--nominal-files', nargs='+', required=True, help='Paths of nominal files')
parser.add_argument('--out', required=True, help='Path to output file (parquet format)')
args = parser.parse_args()


def get_significant_pairs(res_df, nominal_files, group_s=None, fdr=0.05):
    """Significant variant-phenotype pairs based on nominal p-value threshold for each phenotype"""
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] tensorQTL: parsing all significant variant-phenotype pairs', flush=True)
    assert 'qval' in res_df

    # significant phenotypes (apply FDR threshold)
    if group_s is not None:
        df = res_df.loc[res_df['qval'] <= fdr, ['pval_nominal_threshold', 'pval_nominal', 'pval_beta', 'group_id']].copy()
        df.set_index('group_id', inplace=True)
    else:
        df = res_df.loc[res_df['qval'] <= fdr, ['pval_nominal_threshold', 'pval_nominal', 'pval_beta']].copy()
    df.rename(columns={'pval_nominal': 'min_pval_nominal'}, inplace=True)
    signif_phenotype_ids = set(df.index)
    threshold_dict = df['pval_nominal_threshold'].to_dict()
    
    nominal_files = {os.path.basename(i).split('.')[-2]:i for i in nominal_files}

    chroms = sorted(nominal_files.keys(), key=lambda x: int(x.replace('chr', '').replace('X', '23')))
    signif_df = []
    for k,c in enumerate(chroms, 1):
        print(f'  * processing chr. {k}/{len(chroms)}', end='\r', flush=True)
        nominal_df = pd.read_parquet(nominal_files[c])
        nominal_df = nominal_df[nominal_df['pval_nominal'] <= df['pval_nominal_threshold'].max()]
        if group_s is not None:
            nominal_df.insert(1, 'group_id', nominal_df['phenotype_id'].map(group_s))
            nominal_df = nominal_df[nominal_df['group_id'].isin(signif_phenotype_ids)]
            m = nominal_df['pval_nominal'] < nominal_df['group_id'].apply(lambda x: threshold_dict[x])
        else:
            nominal_df = nominal_df[nominal_df['phenotype_id'].isin(signif_phenotype_ids)]
            m = nominal_df['pval_nominal'] < nominal_df['phenotype_id'].apply(lambda x: threshold_dict[x])
        signif_df.append(nominal_df[m])
    print()
    signif_df = pd.concat(signif_df, axis=0)
    if group_s is not None:
        signif_df = signif_df.merge(df, left_on='group_id', right_index=True)
    else:
        signif_df = signif_df.merge(df, left_on='phenotype_id', right_index=True)
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] done', flush=True)
    return signif_df.reset_index(drop=True)


res_df = pd.read_csv(args.permutations, sep='\t')
group_s = pd.read_csv(args.groups, sep='\t', index_col=0, header=None).squeeze('columns')
df = get_significant_pairs(res_df, args.nominal_files, group_s)
df.to_parquet(args.out)